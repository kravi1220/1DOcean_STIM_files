# Block B: Oceananigans LES of a KG6 window (KG6 design 2026-10-05).
#
# Numerics validated against PALM case c02 of Garanaik et al. (2024) (Section 3.5 of the paper):
# WENO(order=9) advection, no explicit subgrid closure (implicit LES), RK3, CFL 0.5.
# KG6 specifics (inputs les_kg6_<B>.nc, read from the working directory):
#   initial T (conservative), S (absolute), u, v from the GOTM k-eps 2 m run at window start
#   hourly surface fluxes exactly as received by that GOTM ocean: non-solar heat Q, shortwave
#   I0 absorbed with depth by Jerlov I (A 0.58, g1 0.35 m, g2 23 m, as GOTM), kinematic wind
#   stress (taux, tauy), upward salt flux ssf (precipitation, ice melt/growth, brine);
#   piecewise-linear in time between the hourly values
#   TEOS-10 equation of state (rho0 1027 kg/m3, cp 3991.87 J/kg/K as GOTM), f at 67.58 N
#   domain 256 m x 256 m x 160 m deep, 1 m grid; no-flux, free-slip bottom
#   turbulence seeded with 1e-3 m/s velocity and 1e-4 K noise in the upper 40 m
#
# Usage: julia --project=. les_kg6_window.jl <B1|B2> [hours] [dx]

using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: on_architecture
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using CUDA, Random, Printf, NCDatasets

win   = get(ARGS, 1, "B1")
hours = parse(Float64, get(ARGS, 2, "120"))
Δ     = parse(Float64, get(ARGS, 3, "1.0"))

ds = NCDataset("les_kg6_$(win).nc")
zp, Tp, Sp, up, vp = (Array{Float64}(ds[k][:]) for k in ("z", "T", "S", "u", "v"))
ts = Array{Float64}(ds["time"][:])
Q, I0, τx, τy, ssf = (Array{Float64}(ds[k][:]) for k in ("Q", "I0", "taux", "tauy", "ssf"))
close(ds)

const ρ₀ = 1027.0
const cₚ = 3991.86795711963
const A_sw, g1_sw, g2_sw = 0.58, 0.35, 23.0
f = 2 * 7.2921e-5 * sind(67.58005324675324)

Lx = Ly = 256.0; Lz = 160.0
Nx, Ny, Nz = round(Int, Lx/Δ), round(Int, Ly/Δ), round(Int, Lz/Δ)
arch = GPU()
grid = RectilinearGrid(arch, size=(Nx, Ny, Nz), x=(0, Lx), y=(0, Ly), z=(-Lz, 0),
                       topology=(Periodic, Periodic, Bounded), halo=(5, 5, 5))

# hourly series -> device arrays; linear interpolation in time inside the kernels
dev(a) = on_architecture(arch, a)
series = (ts=dev(ts), Q=dev(Q), I0=dev(I0), τx=dev(τx), τy=dev(τy), ssf=dev(ssf), Δt=ts[2]-ts[1], n=length(ts))
@inline function interp(t, ts, v, Δt, n)
    k = clamp(unsafe_trunc(Int, t / Δt) + 1, 1, n - 1)
    w = clamp((t - ts[k]) / Δt, 0.0, 1.0)
    return @inbounds (1 - w) * v[k] + w * v[k+1]
end
@inline Jᵀ(x, y, t, p) = -interp(t, p.ts, p.Q,  p.Δt, p.n) / (ρ₀ * cₚ)     # upward T flux
@inline Jˢ(x, y, t, p) =  interp(t, p.ts, p.ssf, p.Δt, p.n)                # upward S flux
@inline Jᵘ(x, y, t, p) = -interp(t, p.ts, p.τx, p.Δt, p.n)                 # upward u flux = -τx/ρ0
@inline Jᵛ(x, y, t, p) = -interp(t, p.ts, p.τy, p.Δt, p.n)
# Shortwave: two-band Jerlov-I absorption integrated over each grid cell (SW fix 2026-10-15).
# Sampling A/g1 exp(z/g1) at cell centres lost the 0.35 m band: only 61 % (2 m) / 44 % (4 m) of I0 was absorbed.
# Cell k gets [Fr(z_top) - Fr(z_bot)] / dz, the bottom cell also the small remainder below -Lz: sum = I0 exactly.
Fr(z) = A_sw * exp(z / g1_sw) + (1 - A_sw) * exp(z / g2_sw)
zf_sw = collect(range(-Lz, 0, length=Nz + 1)); dz_sw = Lz / Nz
w_sw = [(Fr(zf_sw[k+1]) - Fr(zf_sw[k])) / dz_sw for k in 1:Nz]; w_sw[1] += Fr(zf_sw[1]) / dz_sw
@assert abs(sum(w_sw) * dz_sw - 1) < 1e-12
wsw = on_architecture(arch, w_sw)
@inline sw_heating(i, j, k, grid, clock, fields, p) = interp(clock.time, p.ts, p.I0, p.Δt, p.n) / (ρ₀ * cₚ) * @inbounds p.wsw[k]

bcs = (T=FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵀ, parameters=series)),
       S=FieldBoundaryConditions(top=FluxBoundaryCondition(Jˢ, parameters=series)),
       u=FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵘ, parameters=series)),
       v=FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵛ, parameters=series)))
forcing = (T=Forcing(sw_heating, discrete_form=true, parameters=merge(series, (; wsw))),)

buoyancy = SeawaterBuoyancy(equation_of_state=TEOS10EquationOfState(reference_density=ρ₀))
model = NonhydrostaticModel(grid; advection=WENO(order=9), timestepper=:RungeKutta3,
                            tracers=(:T, :S), buoyancy, coriolis=FPlane(f=f),
                            closure=nothing, boundary_conditions=bcs, forcing)

lin(z, zs, vs) = z <= zs[1] ? vs[1] : z >= zs[end] ? vs[end] :
                 (k = searchsortedlast(zs, z); w = (z - zs[k]) / (zs[k+1] - zs[k]); (1 - w) * vs[k] + w * vs[k+1])
Random.seed!(2026)
seed(z) = z > -40 ? 1.0 : 0.0
set!(model, T=(x, y, z) -> lin(z, zp, Tp) + 1e-4 * randn() * seed(z),
            S=(x, y, z) -> lin(z, zp, Sp),
            u=(x, y, z) -> lin(z, zp, up) + 1e-3 * randn() * seed(z),
            v=(x, y, z) -> lin(z, zp, vp) + 1e-3 * randn() * seed(z),
            w=(x, y, z) -> 1e-3 * randn() * seed(z))

simulation = Simulation(model, Δt=0.5, stop_time=hours * hour)
conjure_time_step_wizard!(simulation, cfl=0.5, max_Δt=10.0)

wall0 = Ref(time_ns())
progress(sim) = (@printf("t = %7.2f h  Δt = %5.2f s  max|w| = %.2e  wall = %.1f min\n", time(sim)/hour, sim.Δt,
                         maximum(abs, sim.model.velocities.w), (time_ns()-wall0[])*1e-9/60); flush(stdout))
add_callback!(simulation, progress, TimeInterval(3hour))

u, v, w = model.velocities
T, S = model.tracers.T, model.tracers.S
avg(φ) = Average(φ, dims=(1, 2))
profiles = (; T=avg(T), S=avg(S), u=avg(u), v=avg(v), w=avg(w),
            wT=avg(w * T), wS=avg(w * S), ww=avg(w * w), uu=avg(u * u), vv=avg(v * v),
            wu=avg(w * u), wv=avg(w * v))
prefix = "les_kg6_$(win)_dx$(Δ)"
simulation.output_writers[:profiles] =
    NetCDFWriter(model, profiles, filename="$(prefix)_profiles.nc",
                 schedule=AveragedTimeInterval(1hour, window=1hour), overwrite_files=true)
simulation.output_writers[:slices] =
    NetCDFWriter(model, (; T, w), filename="$(prefix)_xz.nc", indices=(:, 1, :),
                 schedule=TimeInterval(6hour), overwrite_files=true)
simulation.output_writers[:checkpointer] = Checkpointer(model, schedule=TimeInterval(12hour), prefix="$(prefix)_chk", cleanup=true)   # keep only the latest checkpoint (quota)

@info "KG6 LES" win hours Δ Nx Ny Nz f
run!(simulation)
@printf("done: %.1f h simulated in %.1f min wall\n", time(simulation)/hour, (time_ns()-wall0[])*1e-9/60)
