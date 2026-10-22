# Block B: ERA5-driven (interactive-flux) Oceananigans LES at KG6.
#
# Same numerics/setup as les_kg6_window.jl (validated against PALM c02: WENO-9 implicit
# LES, RK3, CFL 0.5; TEOS-10, f at 67.58 N, Jerlov-I shortwave, 256 x 256 x 160 m), but the
# surface fluxes are recomputed every 60 s of model time from ERA5 (KG6 meteo.dat) and the LES's own
# horizontally averaged top-cell temperature, with GOTM's Kondo bulk formulae
# (GotmAirSea.jl; matches the GOTM-validated Python port to 3e-5 W/m2):
#   heat (non-solar) and shortwave (Payne albedo) into the ocean, wind stress, and the
#   precipitation salt flux S_top*P (GOTM: evap = 0). SST for the bulk formulae = in-situ
#   temperature at p = 0 from the mean top-cell CT and SA (as GOTM: gsw_t_from_ct).
# No sea ice: the run stops when the mean top-cell CT reaches the TEOS-10 freezing point.
#
# Usage: julia --project=. les_kg6_era5.jl <input.nc> <hours> <dx> <tag> [pickup]
#   [nofreezestop]: keep running past the freezing point (ocean-only, no ice)
#   pickup: continue from the latest checkpoint <prefix>_chk_iteration*.jld2 (outputs and flux log appended)
#   input.nc: profile z, T (CT), S (SA), u, v and global attribute 'start' (ISO date-time)

using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: on_architecture
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using CUDA, Random, Printf, NCDatasets, GibbsSeaWater
using Dates: DateTime, Millisecond, Second      # not `using Dates`: its `hour` clashes with Oceananigans.Units.hour
include("GotmAirSea.jl"); using .GotmAirSea

infile = ARGS[1]
hours  = parse(Float64, ARGS[2])
Δ      = parse(Float64, ARGS[3])
tag    = ARGS[4]
pickup = length(ARGS) >= 5 && ARGS[5] == "pickup"
freezestop = !(length(ARGS) >= 6 && ARGS[6] == "nofreezestop")   # no-ice run past freezing (as the Block A no-ice 1D runs)

ds = NCDataset(infile)
zp, Tp, Sp, up, vp = (Array{Float64}(ds[k][:]) for k in ("z", "T", "S", "u", "v"))
t_start = DateTime(ds.attrib["start"])
close(ds)
met = read_meteo(get(ENV, "LES_METEO", joinpath(@__DIR__, "..", "..", "experiments", "gotm", "meteo.dat")))   # ERA5 forcing (default KG6)
const LAT, LON = parse(Float64, get(ENV, "LES_LAT", "67.58005324675324")), parse(Float64, get(ENV, "LES_LON", "-23.92983766233766"))
const ρ₀ = 1027.0
const cₚ = 3991.86795711963
const A_sw, g1_sw, g2_sw = 0.58, 0.35, 23.0

Lx = Ly = 256.0; Lz = parse(Float64, get(ENV, "LES_LZ", "160"))   # domain depth (KG6 160 m; LB9 256 m)
Nx, Ny, Nz = round(Int, Lx/Δ), round(Int, Ly/Δ), round(Int, Lz/Δ)
arch = GPU()
grid = RectilinearGrid(arch, size=(Nx, Ny, Nz), x=(0, Lx), y=(0, Ly), z=(-Lz, 0),
                       topology=(Periodic, Periodic, Bounded), halo=(5, 5, 5))

# current surface fluxes, updated every time step on the CPU and copied to the device:
# F = [Jᵀ (upward, K m/s), I0 (W/m2 shortwave entering), Jᵘ, Jᵛ (upward, m2/s2), Jˢ (upward, g/kg m/s)]
F = on_architecture(arch, zeros(5))
p = (; F)
@inline Jᵀ(x, y, t, p) = @inbounds p.F[1]
@inline Jᵘ(x, y, t, p) = @inbounds p.F[3]
@inline Jᵛ(x, y, t, p) = @inbounds p.F[4]
@inline Jˢ(x, y, t, p) = @inbounds p.F[5]
# Shortwave: two-band Jerlov-I absorption integrated over each grid cell (SW fix 2026-10-15).
# Sampling A/g1 exp(z/g1) at cell centres lost the 0.35 m band: only 61 % (2 m) / 44 % (4 m) of I0 was absorbed.
# Cell k gets [Fr(z_top) - Fr(z_bot)] / dz, the bottom cell also the small remainder below -Lz: sum = I0 exactly.
Fr(z) = A_sw * exp(z / g1_sw) + (1 - A_sw) * exp(z / g2_sw)
zf_sw = collect(range(-Lz, 0, length=Nz + 1)); dz_sw = Lz / Nz
w_sw = [(Fr(zf_sw[k+1]) - Fr(zf_sw[k])) / dz_sw for k in 1:Nz]; w_sw[1] += Fr(zf_sw[1]) / dz_sw
@assert abs(sum(w_sw) * dz_sw - 1) < 1e-12
wsw = on_architecture(arch, w_sw)
@inline sw_heating(i, j, k, grid, clock, fields, p) = @inbounds p.F[2] / (ρ₀ * cₚ) * p.wsw[k]

bcs = (T=FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵀ, parameters=p)),
       S=FieldBoundaryConditions(top=FluxBoundaryCondition(Jˢ, parameters=p)),
       u=FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵘ, parameters=p)),
       v=FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵛ, parameters=p)))
forcing = (T=Forcing(sw_heating, discrete_form=true, parameters=(; F, wsw)),)

buoyancy = SeawaterBuoyancy(equation_of_state=TEOS10EquationOfState(reference_density=ρ₀))
model = NonhydrostaticModel(grid; advection=WENO(order=9), timestepper=:RungeKutta3,
                            tracers=(:T, :S), buoyancy, coriolis=FPlane(f=2*7.2921e-5*sind(LAT)),
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

T, S = model.tracers.T, model.tracers.S
Ttop = Field(Average(view(T, :, :, Nz), dims=(1, 2)))
Stop = Field(Average(view(S, :, :, Nz), dims=(1, 2)))
flog = open("les_kg6_era5_$(tag)_fluxes.csv", pickup ? "a" : "w")
pickup || println(flog, "time,CT_top,SA_top,sst_insitu,Q_nonsolar,I0,taux,tauy,precip")
state = Dict(:next_log => 0.0, :stop => false)

function update_fluxes!(sim)
    compute!(Ttop); compute!(Stop)
    ct = Array(interior(Ttop))[1]; sa = Array(interior(Stop))[1]   # averaged top cell (index Nz)
    sst = gsw_t_from_ct(sa, ct, 0.0)
    tnow = t_start + Millisecond(round(Int, 1000 * time(sim)))
    r = bulk(tnow, sst, meteo_at(met, tnow), LON, LAT)
    copyto!(F, [-r.heat / (ρ₀ * cₚ), r.swr_net, -r.taux / ρ₀, -r.tauy / ρ₀, sa * r.precip])
    state[:next_log] = max(state[:next_log], 3600.0 * floor(time(sim) / 3600.0))   # after a pickup
    if time(sim) >= state[:next_log]
        @printf(flog, "%s,%.5f,%.5f,%.5f,%.4f,%.4f,%.6f,%.6f,%.4e\n", tnow, ct, sa, sst, r.heat, r.swr_net, r.taux, r.tauy, r.precip)
        flush(flog); state[:next_log] += 3600.0
    end
    if freezestop && ct <= gsw_ct_freezing(sa, 0.0, 0.0)
        @printf("surface reached the freezing point at %s (CT %.4f, SA %.3f): stopping (no ice model)\n", tnow, ct, sa)
        sim.running = false
    end
    return nothing
end
update_fluxes!(Simulation(model, Δt=1.0, stop_iteration=0))   # fluxes for the first step

simulation = Simulation(model, Δt=0.5, stop_time=hours * hour)
conjure_time_step_wizard!(simulation, cfl=0.5, max_Δt=10.0)
add_callback!(simulation, update_fluxes!, TimeInterval(60.0))   # every 60 s (GOTM updates every 112.5 s)

wall0 = Ref(time_ns())
progress(sim) = (@printf("t = %8.2f h (%s)  Δt = %5.2f s  max|w| = %.2e  wall = %.1f min\n", time(sim)/hour,
                         t_start + Second(round(Int, time(sim))), sim.Δt, maximum(abs, sim.model.velocities.w),
                         (time_ns()-wall0[])*1e-9/60); flush(stdout))
add_callback!(simulation, progress, TimeInterval(24hour))

u, v, w = model.velocities
avg(φ) = Average(φ, dims=(1, 2))
profiles = (; T=avg(T), S=avg(S), u=avg(u), v=avg(v), w=avg(w),
            wT=avg(w * T), wS=avg(w * S), ww=avg(w * w), uu=avg(u * u), vv=avg(v * v), wu=avg(w * u), wv=avg(w * v))
prefix = "les_kg6_era5_$(tag)"
simulation.output_writers[:profiles] =
    NetCDFWriter(model, profiles, filename="$(prefix)_profiles.nc",
                 schedule=AveragedTimeInterval(1hour, window=1hour), overwrite_files=!pickup,
                 global_attributes=Dict("start" => string(t_start), "input" => infile))
simulation.output_writers[:slices] =
    NetCDFWriter(model, (; T, w), filename="$(prefix)_xz.nc", indices=(:, 1, :),
                 schedule=TimeInterval(1day), overwrite_files=!pickup)
simulation.output_writers[:checkpointer] = Checkpointer(model, schedule=TimeInterval(1day), prefix="$(prefix)_chk", cleanup=true)   # keep only the latest checkpoint (quota)

@info "KG6 ERA5-driven LES" infile t_start hours Δ Nx Ny Nz pickup
run!(simulation; pickup)
close(flog)
@printf("done: %.1f h simulated in %.1f min wall\n", time(simulation)/hour, (time_ns()-wall0[])*1e-9/60)
