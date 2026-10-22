# Reproduction package: GOTM with sea-ice thermodynamics at station KG6

Kumar Ravi Prakash and Rohit Shukla. Code, model set-up and inputs for the manuscript
*"GOTM with Sea-Ice Thermodynamics as a Foundation for Ocean–Ice Process Studies: Testing a
Nonlocal Mass-Flux Closure, Large-Eddy Simulation and Mushy-Layer Sea Ice at a Subpolar Station"*
(submitted to the *Journal of Advances in Modeling Earth Systems*).

With this package you can build the four models, rerun every experiment in the paper and
regenerate every figure and the numbers quoted in the text. It contains **no model output**:
all of it is produced by the run scripts below.

## Contents

| Folder | Contents |
|---|---|
| `models/gotm/` | GOTM v6.0.0-211 (`84f8d8e9`) with STIM (`c54fc32`), CVMix, GSW and flexout; Icepack (CICE-Consortium `a348e4c`, unmodified) added to STIM as the ice model `icepack`, and the Lebedev fix. |
| `models/icepack_pipe/` | Stand-alone Icepack column (same STIM–Icepack wrapper) driven over stdin/stdout; used by the MPAS-ADC and LES couplings. |
| `models/mpas_adc/` | MPAS-Ocean with the ADC closure of Garanaik et al. (2024) (fork `amrapallig/MPAS-Model`, commit `7b8141e3`), ocean core only, with the two salt-conservation fixes and two build fixes. |
| `models/les/` | Oceananigans LES scripts (ocean only, coupled to Icepack, fixed-flux windows), the GOTM bulk-flux port `GotmAirSea.jl`, and the exact Julia environment (`Project.toml`, `Manifest.toml`). |
| `coupling/` | Python: GOTM Kondo bulk fluxes (`gotm_airsea.py`), the hourly MPAS-Ocean driver (`mpas_driver.py`), the Icepack column interface (`icepack_pipe.py`) and the Icepack server for the LES (`ice_server.py`). |
| `patches/` | Every change we made to third-party code, as unified diffs against the upstream versions above (already applied in `models/`). |
| `experiments/gotm/` | All GOTM configurations (`gotm.yaml`, `icepack_in`) by experiment block, with the ERA5 forcing (`meteo.dat`, `meteo_snow.dat`) and initial profiles. |
| `experiments/mpas_adc/` | All MPAS-ADC cases (`init.nc`, `graph.info`, `namelist.ocean`, `streams.ocean`, `case.json`, `icepack_in` or `window_fluxes.nc`). |
| `experiments/les/` | LES initial-state/forcing files (`les_kg6_*.nc`) and the Icepack set-up of the four coupled-run segments. |
| `analysis/` | Scripts that compute the tables and the numbers quoted in the paper and draw Figures 1–11. |
| `data/` | OSI SAF ice concentration at the KG6 grid cell (`osisaf_kg6.csv`) and daily regional fields 1 Dec 2018–31 May 2019 (`osisaf_region_*.npz`). The CTD validation casts are not included (see below). |
| `scripts/` | Build, run and analysis scripts. |

The code that created the input files (initial profiles, ERA5 forcing, model configurations,
LES input files) and the code that downloaded ERA5 and OSI SAF data are not included; the
input files themselves are.

## Experiments

Blocks as in Table 1 of the paper. All runs: station KG6 (67.58°N, 23.93°W), ERA5 forcing
through GOTM's Kondo bulk formulae, 6 Aug 2018 12:00 – 10 Aug 2019 12:00 unless stated.

| Block | Where | Runs |
|---|---|---|
| R – reference | `experiments/gotm/blockR/` | GOTM k-ε / KPP with STIM-Winton (R1, R2) and without ice (R3, R4) at time steps 1800–28.125 s; k-ε with STIM-MyLake and the patched STIM-Lebedev (`R1_lebedev_fix`; the unpatched behaviour, no ice at any time step, is obtained by building with `patches/gotm_stim_lebedev_fix.patch` reversed); 5-member initial-condition ensemble (`R5_dt112.5`). |
| A – closures | `experiments/gotm/blockA/`, `experiments/mpas_adc/blockA/`, `blockA_dt/` | GOTM k-ε, KPP and MPAS-ADC on uniform 2, 5, 10 m grids, no ice; ADC outer-step ladder 30 min – 18.75 s (1700 h, against the 9.375 s run). |
| B – LES | `experiments/les/` | Windows B1 (28 Nov – 3 Dec 2018) and B2 (8–13 Feb 2019) at 1/2/4 m; run from the 15 Oct 2018 cast to 16 Feb 2019; full year ocean only and coupled to Icepack, 2 m. |
| W – 1D vs LES | `experiments/gotm/blockW/`, `experiments/mpas_adc/blockW/` | k-ε, KPP, ADC at 2, 5, 10 m over B1 and B2 with the LES initial state and prescribed fluxes. |
| Oct-15 | `experiments/gotm/oct15/`, `experiments/mpas_adc/oct15/` | k-ε, KPP, ADC at 2, 5, 10 m from the 15 Oct 2018 cast (the ADC 2 m run stops in its first hour, Section 4.2.3). |
| C – sea ice | `experiments/gotm/blockC/` | GOTM k-ε with Icepack BL99 and mushy-layer thermodynamics. |
| D – interaction | `experiments/gotm/blockD/`, `experiments/mpas_adc/blockD/` | k-ε (2 and 5 m) and ADC (5 m) coupled to Icepack mushy. |

## Requirements

- Fortran and C compilers with MPI (`mpif90`, `mpicc`), CMake ≥ 3.10, NetCDF-C and NetCDF-Fortran.
  Used: GCC 12.3, Open MPI 4.1, netCDF-Fortran 4.6.1, CMake 3.26 (RHEL 8).
- ParallelIO (PIO) 2.6.2 (`NCAR/ParallelIO`, tag `pio2_6_2`) built with CMake and `-DWITH_PNETCDF=OFF`, for MPAS-Ocean.
- Python 3 with `numpy`, `netCDF4`, `gsw`, `matplotlib`, and `cartopy` for Figure 1
  (used: Python 3.11, numpy 2.4, netCDF4 1.7, gsw 3.6, matplotlib 3.11, cartopy 0.26; cartopy
  downloads Natural Earth coastlines on first use).
- Julia 1.10 and an NVIDIA GPU (CUDA) for the LES; `julia --project=models/les -e 'using Pkg; Pkg.instantiate()'`
  installs the exact package versions of `Manifest.toml` (Oceananigans 0.113.5, CUDA.jl 6.4.2).

## 1. Build

```bash
scripts/build_gotm.sh                       # -> build/gotm/gotm
scripts/build_icepack_pipe.sh               # -> models/icepack_pipe/icepack_pipe (needs the GOTM build)
NETCDF=<netcdf-fortran prefix> PIO=<pio prefix> scripts/build_mpas_ocean.sh   # -> models/mpas_adc/ocean_model
julia --project=models/les -e 'using Pkg; Pkg.instantiate()'
```

## 2. Run

```bash
scripts/run_gotm.sh build/gotm/gotm         # all 65 GOTM runs (minutes each)
scripts/run_mpas.sh                         # all 20 MPAS-ADC cases (full-year cases about 1.5-2 days each, serial)
scripts/run_les.sh windows                  # LES, one GPU each; see the script for wall times and segments
scripts/run_les.sh oct15
scripts/run_les.sh alone
ICE_ENV=<shell file with python3 and the Fortran runtime> scripts/run_les.sh coupled
```

`run_mpas.sh` also takes single case folders, so the cases can be spread over a cluster.
The coupled LES runs in four segments that begin and end without ice, because the Icepack
column is restarted empty at each segment start (`experiments/les/ice_cpl2m_<seg>/icepack_in`).

## 3. Analyse

```bash
scripts/run_analysis.sh     # tables -> results/tables/, Figures 1-11 -> results/figures/
```

| Paper | File (results/figures) | Script |
|---|---|---|
| Fig. 1 | `K0_KG6_osisaf_map.png` | `analysis/fig_map.py` |
| Figs. 2–6 | `K1`–`K5` | `analysis/fig_main.py` |
| Figs. 7–11 | `K7`–`K11` | `analysis/fig_summary.py` |
| Sections 4.1–4.7, Appendix A | `results/tables/*.txt`, `K_numbers.json`, `K_summary_numbers.txt` | the other scripts in `analysis/` |

## Reproducibility notes

- GOTM: rebuilding from `models/gotm` with the compiler and libraries above reproduces the
  output of all 65 GOTM runs bit for bit; MPAS-ADC and the Icepack column reproduce the original
  runs bit for bit over the period tested; the analysis reproduces the paper's tables and figures
  (details in `VERIFICATION.md`). Other compilers or library versions give small round-off differences.
- The Block W GOTM runs start from the window state in `restart_initial.nc` (LES initial profiles,
  quiescent turbulence), which `run_gotm.sh` copies to `restart.nc` before each run.
- LES: GPU runs are deterministic for a given GPU, CUDA and package version, but turbulence is
  chaotic: a restart at a different time, or different hardware, gives a different but
  statistically equivalent realisation (in our runs, restarts changed the profiles by at most
  0.04°C before ice formed). The published segments are listed in `scripts/run_les.sh`.
- The B1 1 m LES (used only for the grid-convergence statement, Section 3.5) was run with the
  earlier shortwave discretisation in the archived results; this package runs it with the
  cell-integrated shortwave used in all other runs. Shortwave in B1 is 1–7 W m⁻², so the
  difference is negligible.

## Changes to third-party code (`patches/`)

- `gotm_stim_icepack_coupling.patch`: Icepack as a STIM ice model and its two-way coupling in `gotm.F90`.
- `gotm_stim_lebedev_fix.patch`: STIM-Lebedev freezing-degree-day accumulation (Section 4.4.1).
- `mpas_adc_salt_conservation.patch`: the two salt-flux errors in the published ADC code (Section 3.4).
- `mpas_build_fixes.patch`: PIO library detection, and CVMix/BGC checks that do not need git metadata.

## Licence

Our code (`coupling/`, `analysis/`, `scripts/`, `models/les/`, `models/icepack_pipe/`, the
patches and configuration files) is released under the MIT licence (`LICENSE`). Bundled
third-party code keeps its own licence: GOTM (GPL v2, `models/gotm/COPYING`), STIM
(`models/gotm/extern/stim/LICENSE`), Icepack (`models/gotm/extern/stim/src/models/icepack/LICENSE.pdf`),
CVMix and GSW (their folders), MPAS (`models/mpas_adc/LICENSE`).

## CTD observations (not included)

The KG6 CTD casts are held by the Marine and Freshwater Research Institute (MFRI), Iceland, and
are not distributed with this package; they can be requested from MFRI. The model runs do not
need them. The scoring against the casts and the CTD points in the figures do: place the casts in
`data/ctd_profiles_KG6.nc` (NetCDF, dimensions `time` and `depth`; variables `time` (with CF units),
`depth` (m, positive down), `temperature` (in-situ, °C) and `salinity` (practical), NaN where not
measured), with the casts of 15 Oct 2018, 15 Feb, 18 May and 10 Aug 2019. The initial-condition
files of the experiments contain the 6 Aug and 15 Oct 2018 casts as model input.

## Data sources

- CTD: Marine and Freshwater Research Institute (MFRI), Iceland, repeat-hydrography station KG6 (not included).
- ERA5 (Hersbach et al., 2020), Copernicus Climate Change Service; contains modified Copernicus information.
- OSI SAF Global Sea Ice Concentration Climate Data Record OSI-450-a and its operational continuation, EUMETSAT.

## Use of AI tools

This package was assembled with the help of Claude (Anthropic): the model couplings, run and
analysis scripts and this documentation were written with AI assistance and checked against
the model output by rerunning (see `VERIFICATION.md`). The authors are responsible for the content.
