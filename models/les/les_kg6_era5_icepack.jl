# Oceananigans LES x Icepack (mushy) at KG6: ERA5-driven, interactive fluxes, hourly ice-ocean coupling.
#
# Same numerics, domain, forcing and initial state as les_kg6_era5.jl (ocean-only LES), but coupled from the start
# (checkpoint argument 'none') or continued from its OWN checkpoint, and couples a persistent Icepack column (ice_server.py -> icepack_pipe, the same column
# model and settings as the GOTM and MPAS-ADC Icepack runs). Coupling, as in the MPAS-ADC driver:
#   every hour: hour-start mean top-cell T (in situ), S (practical) and the top-cell thickness are
#     sent to Icepack, which returns the ice state and hour-mean ice->ocean fluxes for that hour;
#   every 60 s: open-water bulk fluxes (GOTM Kondo, LES SST) weighted by (1 - aice_mean), plus
#     heat  fhocn + frazil heat, shortwave fswthru, freshwater fresh and salt fsalt from Icepack;
#     wind stress (1 - aice_mean) * tau (motionless ice); rain/snow on ice go to Icepack.
#   upward salt flux = SA_top ((1 - a) P + fresh/1000) - 1000 fsalt / rho0   (as in the window input files)
# No freezing stop: Icepack handles freeze-up (frazil from the supercooled top cell).
#
# Usage: julia --project=. les_kg6_era5_icepack.jl <input.nc> <stop_hours> <dx> <tag> <checkpoint.jld2|none> <ice_dir>
#   ice_dir: Icepack run folder whose icepack_in start (year_init, istep0) equals the model start time
#   stop_hours: model time since the input 'start' at which to stop

using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: on_architecture
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using CUDA, Random, Printf, NCDatasets, GibbsSeaWater, JLD2
using Dates: DateTime, Millisecond, Second, Hour, format
include("GotmAirSea.jl"); using .GotmAirSea

infile, hours, Δ, tag, chk, icedir = ARGS[1], parse(Float64, ARGS[2]), parse(Float64, ARGS[3]), ARGS[4], ARGS[5], ARGS[6]
ds = NCDataset(infile); t_start = DateTime(ds.attrib["start"]); close(ds)
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
F = on_architecture(arch, zeros(5)); p = (; F)
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
buoyancy = SeawaterBuoyancy(equation_of_state=TEOS10EquationOfState(reference_density=ρ₀))
model = NonhydrostaticModel(grid; advection=WENO(order=9), timestepper=:RungeKutta3,
                            tracers=(:T, :S), buoyancy, coriolis=FPlane(f=2*7.2921e-5*sind(LAT)),
                            closure=nothing, boundary_conditions=bcs, forcing=(T=Forcing(sw_heating, discrete_form=true, parameters=(; F, wsw)),))
if chk == "none"      # fresh start from the input profile, exactly as les_kg6_era5.jl (same seed and noise)
    ds = NCDataset(infile); zp, Tp, Sp, up, vp = (Array{Float64}(ds[k][:]) for k in ("z", "T", "S", "u", "v")); close(ds)
    lin(z, zs, vs) = z <= zs[1] ? vs[1] : z >= zs[end] ? vs[end] :
                     (k = searchsortedlast(zs, z); w = (z - zs[k]) / (zs[k+1] - zs[k]); (1 - w) * vs[k] + w * vs[k+1])
    Random.seed!(2026)
    seed(z) = z > -40 ? 1.0 : 0.0
    set!(model, T=(x, y, z) -> lin(z, zp, Tp) + 1e-4 * randn() * seed(z), S=(x, y, z) -> lin(z, zp, Sp),
                u=(x, y, z) -> lin(z, zp, up) + 1e-3 * randn() * seed(z), v=(x, y, z) -> lin(z, zp, vp) + 1e-3 * randn() * seed(z),
                w=(x, y, z) -> 1e-3 * randn() * seed(z))
else                 # continue this coupled run from its own checkpoint (fields stored with halos; RK3 keeps no history)
    jldopen(chk) do f
        g = "simulation/model/"
        for (name, fld) in (("velocities/u", model.velocities.u), ("velocities/v", model.velocities.v), ("velocities/w", model.velocities.w),
                            ("tracers/T", model.tracers.T), ("tracers/S", model.tracers.S))
            copyto!(parent(fld), f[g * name * "/data"])
        end
        c = f[g * "clock"]
        model.clock.time = getproperty(c, :time); model.clock.iteration = getproperty(c, :iteration)
    end
end
@info "initial state" chk time=t_start + Second(round(Int, model.clock.time))

T, S = model.tracers.T, model.tracers.S
Ttop = Field(Average(view(T, :, :, Nz), dims=(1, 2))); Stop = Field(Average(view(S, :, :, Nz), dims=(1, 2)))
function top_state()
    compute!(Ttop); compute!(Stop)
    return Array(interior(Ttop))[1], Array(interior(Stop))[1]
end

# ---- Icepack column (separate process; its own environment) ----
# ICE_ENV: optional shell file sourced before starting the Icepack server (compiler runtime, Python with numpy/netCDF4)
env = get(ENV, "ICE_ENV", "/dev/null")
ice_server = joinpath(@__DIR__, "..", "..", "coupling", "ice_server.py")
tnow0 = t_start + Second(round(Int, model.clock.time))
ice = open(`bash -c "source $env >/dev/null 2>&1; exec python3 -u $ice_server $icedir $(format(tnow0, "yyyy-mm-ddTHH:MM:SS"))"`, "r+")
while true
    l = readline(ice); println("ice_server: ", l); flush(stdout)
    l == "READY" && break
    eof(ice) && error("ice_server ended during start-up")
end
const ICEKEYS = ("aice", "vice", "vsno", "Tsfc", "Sice", "fhocn", "fswthru", "fresh", "fsalt", "frazil", "aice_mean")
icestate = Dict(k => 0.0 for k in ICEKEYS)
ilog = open("les_kg6_icepack_$(tag)_ice.csv", "w"); println(ilog, "time,CT_top,SA_top," * join(ICEKEYS, ","))
flog = open("les_kg6_icepack_$(tag)_fluxes.csv", "w"); println(flog, "time,CT_top,SA_top,sst_insitu,Q_ocean,I0_ocean,taux,tauy,ssf,aice_mean")

function ice_hour!(sim)
    ct, sa = top_state()
    tnow = t_start + Millisecond(round(Int, 1000 * time(sim)))
    sst = gsw_t_from_ct(sa, ct, 0.0); sp = gsw_sp_from_sa(sa, 0.0, LON, LAT)
    println(ice, format(tnow, "yyyy-mm-ddTHH:MM:SS"), " ", sst, " ", sp, " ", Δ); flush(ice)
    v = parse.(Float64, split(readline(ice)))
    length(v) == length(ICEKEYS) || error("ice_server: bad reply at $tnow")
    for (k, x) in zip(ICEKEYS, v); icestate[k] = x; end
    @printf(ilog, "%s,%.5f,%.5f,%s\n", tnow, ct, sa, join((@sprintf("%.6e", x) for x in v), ",")); flush(ilog)
    return nothing
end

state = Dict(:next_log => 0.0)
function update_fluxes!(sim)
    ct, sa = top_state()
    sst = gsw_t_from_ct(sa, ct, 0.0)
    tnow = t_start + Millisecond(round(Int, 1000 * time(sim)))
    r = bulk(tnow, sst, meteo_at(met, tnow), LON, LAT)
    a = icestate["aice_mean"]; ow = 1 - a
    Q = ow * r.heat + icestate["fhocn"] + icestate["frazil"]
    I0 = ow * r.swr_net + icestate["fswthru"]
    ssf = sa * (ow * r.precip + icestate["fresh"] / 1000) - 1000 * icestate["fsalt"] / ρ₀
    copyto!(F, [-Q / (ρ₀ * cₚ), I0, -ow * r.taux / ρ₀, -ow * r.tauy / ρ₀, ssf])
    state[:next_log] = max(state[:next_log], 3600.0 * floor(time(sim) / 3600.0))
    if time(sim) >= state[:next_log]
        @printf(flog, "%s,%.5f,%.5f,%.5f,%.4f,%.4f,%.6f,%.6f,%.4e,%.4f\n", tnow, ct, sa, sst, Q, I0, ow * r.taux, ow * r.tauy, ssf, a)
        flush(flog); state[:next_log] += 3600.0
    end
    return nothing
end

simulation = Simulation(model, Δt=0.5, stop_time=hours * hour)
ice_hour!(simulation); update_fluxes!(simulation)                # fluxes for the first hour
conjure_time_step_wizard!(simulation, cfl=0.5, max_Δt=10.0)
add_callback!(simulation, ice_hour!, TimeInterval(3600.0))         # registered first: ice before fluxes
add_callback!(simulation, update_fluxes!, TimeInterval(60.0))

wall0 = Ref(time_ns())
progress(sim) = (@printf("t = %8.2f h (%s)  Δt = %5.2f s  max|w| = %.2e  aice = %.3f  hi = %.3f m  wall = %.1f min\n",
                         time(sim)/hour, t_start + Second(round(Int, time(sim))), sim.Δt, maximum(abs, sim.model.velocities.w),
                         icestate["aice"], icestate["aice"] > 0 ? icestate["vice"] / icestate["aice"] : 0.0,
                         (time_ns()-wall0[])*1e-9/60); flush(stdout))
add_callback!(simulation, progress, TimeInterval(24hour))

u, v, w = model.velocities
avg(φ) = Average(φ, dims=(1, 2))
profiles = (; T=avg(T), S=avg(S), u=avg(u), v=avg(v), w=avg(w),
            wT=avg(w * T), wS=avg(w * S), ww=avg(w * w), uu=avg(u * u), vv=avg(v * v), wu=avg(w * u), wv=avg(w * v))
prefix = "les_kg6_icepack_$(tag)"
simulation.output_writers[:profiles] =
    NetCDFWriter(model, profiles, filename="$(prefix)_profiles.nc", schedule=AveragedTimeInterval(1hour, window=1hour),
                 overwrite_files=true, global_attributes=Dict("start" => string(t_start), "restart_from" => chk))
simulation.output_writers[:slices] =
    NetCDFWriter(model, (; T, w), filename="$(prefix)_xz.nc", indices=(:, 1, :), schedule=TimeInterval(1day), overwrite_files=true)
simulation.output_writers[:checkpointer] = Checkpointer(model, schedule=TimeInterval(1day), prefix="$(prefix)_chk", cleanup=true)   # keep only the latest checkpoint (quota)

@info "KG6 LES x Icepack" infile t_start hours Δ Nx Ny Nz icedir
run!(simulation)
println(ice, "END"); flush(ice); close(ilog); close(flog)
@printf("done: %.1f h simulated (to %s) in %.1f min wall\n", time(simulation)/hour, t_start + Second(round(Int, time(simulation))),
        (time_ns()-wall0[])*1e-9/60)
