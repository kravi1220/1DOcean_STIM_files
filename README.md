#  "Beyond the No-Physics Control: A Verification Framework for 1D Ocean-Ice Modeling"

This repo is for the manuscript submitted to the
*Journal of Advances in Modeling Earth Systems* (JAMES). 

This repo lets you **build GOTM+STIM, run every model configuration used
in the paper, and regenerate every figure and table from that output**. It
does **not** include the code that built the initial-condition, forcing, or
GOTM-config files in the first place, nor any raw model output. 
| Folder | Contents |
|---|---|
| `gotm_source/` | GOTM v6.0.0 source (+ STIM ice module, CVMix, gsw, flexout), with git history and build directories stripped. See `gotm_source/VERSION.txt` for exact component commits and the build configuration used. |
| `run_gotm/run_gotm_case.sh` | Runs the `gotm` binary over every config in a run directory. |
| `model_inputs/case_LB9/`, `model_inputs/case_KG6/` | Every `gotm.yaml` config used in the paper, plus the initial-condition (`t_prof_*.dat`/`s_prof_*.dat`) and surface-forcing (`meteo*.dat`) files each one reads. No model output. |
| `station_screening/screen_stations.py` | Screens all 87 network stations for the coldest-relative-to-freezing occupation (Section 2.1) — confirms LB9 ranks 1st, KG6 2nd. |
| `scoring/` | Scores model output against real CTD checkpoints and a persistence baseline; computes the significance/noise-floor checks and the sweep convergence correction. Writes `tables_output/`. |
| `figures/` | Builds every figure from model output. Writes `figures_output/`. Imports the shared `seaice1d/` plotting-style package. |
| `data/` | CTD profile data for all 87 network stations (`ctd_profiles_<STATION>.nc`) and cached OSI-SAF sea-ice concentration extracts used by the figures. |
| `tables_output/`, `figures_output/` | Pre-computed reference outputs, regenerated fresh from this exact package and verified byte-identical to the paper's own working tables (image files may differ by a few pixels from run to run due to font/anti-aliasing rendering, not data). |


## Software requirements

- **A Fortran compiler, CMake, and a NetCDF-Fortran library** to build GOTM.
- **Python 3** with `numpy`, `netCDF4`, `gsw` (TEOS-10), `matplotlib`, `scipy`.
- No internet access is required for any step below.

## 1. Build GOTM

```bash
cd gotm_source
cmake -B build -DCMAKE_BUILD_TYPE=Release \
      -DGOTM_USE_STIM=ON -DGOTM_USE_CVMIX=ON -DGOTM_USE_FABM=OFF \
      -DSTIM_WINTON=ON -DSTIM_LEBEDEV=ON -DSTIM_MYLAKE=ON -DSTIM_BASAL_MELT=ON
cmake --build build
# produces gotm_source/build/gotm
```

## 2. Run the model

`run_gotm_case.sh` runs every `gotm.yaml`-containing subdirectory of a given
run directory (this GOTM build exits with code 1 even on full success —
the script judges success by the presence of `all.nc`, not the exit code):

```bash
GOTM=$(pwd)/gotm_source/build/gotm

# Main-event factorial + no-ice controls + IC ensembles (LB9)
./run_gotm/run_gotm_case.sh model_inputs/case_LB9/runs "$GOTM"

# Generalization tests
./run_gotm/run_gotm_case.sh model_inputs/case_LB9/runs_event2 "$GOTM"
./run_gotm/run_gotm_case.sh model_inputs/case_LB9/runs_fairall "$GOTM"
./run_gotm/run_gotm_case.sh model_inputs/case_LB9/runs_restart_oct2018 "$GOTM"
./run_gotm/run_gotm_case.sh model_inputs/case_LB9/runs_restart_oct2018_full "$GOTM"
./run_gotm/run_gotm_case.sh model_inputs/case_LB9/runs_restart2019 "$GOTM"
for off in m3 m2 m1 p1 p2 p3; do
  ./run_gotm/run_gotm_case.sh model_inputs/case_LB9/runs_sweep_$off "$GOTM"
done

# Numerical-convergence checks (Sections 4.1.2, 4.4.2, 4.4.3; dt halved down
# to 14 s for KG6). runs_sensitivity/ itself holds several single-config
# runs directly (dt900, dt450, dt225, keps_dt450, kpp_dt450, ...) plus two
# nested containers that each need their own call:
./run_gotm/run_gotm_case.sh model_inputs/case_LB9/runs_sensitivity "$GOTM"
./run_gotm/run_gotm_case.sh model_inputs/case_LB9/runs_sensitivity/correctedpos_dt450 "$GOTM"
./run_gotm/run_gotm_case.sh model_inputs/case_LB9/runs_sensitivity/correctedpos_fairall_dt450 "$GOTM"
./run_gotm/run_gotm_case.sh model_inputs/case_LB9/runs_dt450 "$GOTM"
./run_gotm/run_gotm_case.sh model_inputs/case_LB9/runs_restart2019_dt450 "$GOTM"
./run_gotm/run_gotm_case.sh model_inputs/case_LB9/runs_restart_oct2018_dt450 "$GOTM"
for off in m1 m2 m3 p1 p2 p3; do
  ./run_gotm/run_gotm_case.sh model_inputs/case_LB9/runs_sweep_${off}_dt450 "$GOTM"
done

# KG6 (non-marginal comparison site)
./run_gotm/run_gotm_case.sh model_inputs/case_KG6/runs "$GOTM"
./run_gotm/run_gotm_case.sh model_inputs/case_KG6/runs_dt450 "$GOTM"
./run_gotm/run_gotm_case.sh model_inputs/case_KG6/runs_sensitivity "$GOTM"
```

Each `runs*` directory writes its own `all.nc` alongside its `gotm.yaml`.

## 3. Score every run

Run from the package root (each script resolves its own paths from its
location, so this works regardless of your current shell's directory, but
the examples below assume you're at the package root):

```bash
python station_screening/screen_stations.py   # confirms LB9/KG6 ranking

python scoring/s16_LB9_analysis.py             # -> T18
python scoring/s19_LB9_noice_comparison.py     # -> T19  (Table 1)
python scoring/s20_LB9_table_significant_differences.py  # -> T20 (Table 2)
python scoring/s22_LB9_significance_check.py   # -> T21  (Table 3)
python scoring/s31_LB9_upper_ocean_rmse.py     # -> T33/T34 (Table 4)
python scoring/s33_LB9_ic_noise_floor_depthvarying.py     # -> T35 (Table 2 noise floor)
python scoring/s23_LB9_event2_analysis.py      # -> T23/T24
python scoring/s25_LB9_fairall_analysis.py     # -> T26/T27 (Table 5)
python scoring/s26_LB9_sweep_analysis.py       # -> T28, uncorrected (dt=1800s)
python scoring/s34_LB9_sweep_convergence_fix.py           # -> T28, corrected (Table 6)
python scoring/s28_LB9_kg6_analysis.py         # -> T29/T30 (Table 7)
```

**Run `s26` before `s34`** — `s34` reconverges specific offsets in `s26`'s
own output (the -3/-2 degC sweep offsets originally looked engaged at the
default 1800 s time step; rerun at 450 s, ice thickness is 0.0 m at every
offset — see `s34`'s own docstring for the full account).

## 4. Build every figure

```bash
python figures/s21_LB9_figures_v2.py     # Figures 1, 2, 3 (F24 map, F23 all-experiments, F28 closure-vs-ice)
python figures/s32_LB9_deep_profile_figure.py  # Figure 4 (F35 deep-profile comparison)
python figures/s30_LB9_fairall_figure.py       # Figure 5 (F34 Fairall comparison)
python figures/s27_LB9_sweep_figure.py         # Figure 6 (F32 sweep, run after s34 above)
python figures/s29_LB9_kg6_figure.py           # Figure 7 (F33 KG6 comparison)
python figures/s24_LB9_event2_figure.py        # supplementary: second-event engagement (not embedded in the current manuscript)
```

`s21` and `s24` also write a couple of intermediate figures (`F19`, `F25`,
`F31`) from earlier in the pipeline's development that are no longer
embedded in the current manuscript; they're harmless to regenerate and are
left in for completeness of the working pipeline.

## Figure/table cross-reference (current manuscript)

| Manuscript | File | Built by |
|---|---|---|
| Figure 1 | `F24_LB9_map.png` | `figures/s21_LB9_figures_v2.py` |
| Figure 2 | `F23_LB9_all_experiments.png` | `figures/s21_LB9_figures_v2.py` |
| Figure 3 | `F28_LB9_closure_vs_ice_effects.png` | `figures/s21_LB9_figures_v2.py` |
| Figure 4 | `F35_LB9_deep_profile_comparison.png` | `figures/s32_LB9_deep_profile_figure.py` |
| Figure 5 | `F34_LB9_fairall_experiments.png` | `figures/s30_LB9_fairall_figure.py` |
| Figure 6 | `F32_LB9_sweep_dose_response.png` | `figures/s27_LB9_sweep_figure.py` |
| Figure 7 | `F33_LB9_kg6_experiments.png` | `figures/s29_LB9_kg6_figure.py` |
| Table 1 | `T19_LB9_noice_comparison.csv` | `scoring/s19_LB9_noice_comparison.py` |
| Table 2 | `T20_LB9_significant_differences.csv`, `T35_LB9_ic_noise_floor_depthvarying.csv` | `scoring/s20_LB9_table_significant_differences.py`, `scoring/s33_LB9_ic_noise_floor_depthvarying.py` |
| Table 3 | `T21_LB9_per_date_rmse.csv` | `scoring/s22_LB9_significance_check.py` |
| Table 4 | `T33_LB9_top50_per_date_rmse.csv` | `scoring/s31_LB9_upper_ocean_rmse.py` |
| Table 5 | `T26_LB9_fairall_skill.csv`, `T27_LB9_fairall_noice_comparison.csv` | `scoring/s25_LB9_fairall_analysis.py` |
| Table 6 | `T28_LB9_sweep_summary.csv` | `scoring/s26_LB9_sweep_analysis.py` then `scoring/s34_LB9_sweep_convergence_fix.py` |
| Table 7 | `T30_LB9_kg6_noice_comparison.csv` | `scoring/s28_LB9_kg6_analysis.py` |
| Table A1 | (static — no script; the experiment design itself) | — |


## License

Code in this package is released under the MIT License (`LICENSE`). GOTM,
STIM, and CVMix are third-party and separately licensed — see
`gotm_source/VERSION.txt` and https://github.com/gotm-model/code.

## Citation

If you use this package, please cite the associated manuscript (full
citation to be added on acceptance).
