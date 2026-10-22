# Verification of this package

Checks made on 22 October 2026 on the Elja cluster (RHEL 8, GCC 12.3, Open MPI 4.1, netCDF-Fortran 4.6.1),
comparing this package against the model output and figures of the study.

## Models

| Check | Result |
|---|---|
| GOTM built from `models/gotm` (`scripts/build_gotm.sh`); all 65 configurations in `experiments/gotm` rerun | every variable of every `all.nc` bit-identical to the original runs; `icepack_column.dat` of the four Icepack runs identical |
| `icepack_pipe` built with `scripts/build_icepack_pipe.sh`, driven with 48 h of freezing forcing | output identical to the original executable (ice grows to full cover, 0.21 m) |
| MPAS-Ocean built from `models/mpas_adc` (`scripts/build_mpas_ocean.sh`, PIO 2.6.2 without the empty `libpio.a` the original build needed); `coupling/mpas_driver.py` run for 24 h on `blockA/KG6_dz5`, `blockD/KG6_dz5_icepack`, `blockW/B1_dz2`, `oct15/KG6oct_dz5` | hourly profiles, ADC moments, surface fluxes and the Icepack log bit-identical to the first 24 h of the original runs |
| LES | not rerun (about 6 GPU-days); the scripts differ from those of the published runs only in file locations |

## Analysis

`scripts/run_analysis.sh` was run on the original model output arranged in the package layout.

| Output | Result |
|---|---|
| `blockA_score.txt`, `compare_fullyear.txt`, `blockD_table.txt`, `adc_dt_ladder.txt`, `K_summary_numbers.txt` | identical to the study's tables |
| `compare_oct15.txt` | identical, except that the two rows of runs with the original (not salt-corrected) ADC code are not produced |
| `icepack_salt_budget.txt` | residuals 0.16 % (mushy) and 0.67 % (BL99) of the peak ice salt, as in Appendix A |
| Figures 2–6, 8–11 (`K1`–`K5`, `K8`–`K11`) | pixel-identical to the manuscript figures |
| Figure 1 (`K0`) | same data; longitude/latitude labels drawn horizontally by the cartopy version used here |
| Figure 7 (`K7`) | now drawn from `blockW_compare.txt`; the manuscript version used window scores typed in before the LES shortwave correction (differences up to 0.003°C in temperature and RMS) |
| `les_convergence.txt` | computed from the corrected 2 and 4 m window LES (the original table used the earlier runs; differences in the third decimal) |
