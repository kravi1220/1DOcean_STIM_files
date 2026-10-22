#!/bin/bash
# Run MPAS-Ocean + ADC cases with the hourly flux exchange of coupling/mpas_driver.py.
# Each case runs serially (one MPI rank); a full year at the 9.375 s outer step takes about 1.5-2 days.
# Needs: models/mpas_adc/ocean_model (scripts/build_mpas_ocean.sh), mpirun, Python 3 with numpy, netCDF4;
#        models/icepack_pipe/icepack_pipe for blockD (scripts/build_icepack_pipe.sh).
# Usage: scripts/run_mpas.sh [case_dir ...]   (default: every case below experiments/mpas_adc)
# Output per case: hourly/<site>_adc.nc (profiles), hourly/<site>_fluxes.csv, and for Icepack icepack_column.dat.
# The oct15/KG6oct_dz2 case stops in its first hour by design of the published code (paper Section 4.2.3).
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cases=("$@"); [ ${#cases[@]} -eq 0 ] && cases=($(dirname $(find "$ROOT/experiments/mpas_adc" -name case.json | sort)))
for c in "${cases[@]}"; do
  echo "== $c"; python3 "$ROOT/coupling/mpas_driver.py" "$c" || echo "case $c stopped (see message above)"
done
