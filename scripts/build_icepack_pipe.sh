#!/bin/bash
# Build the stand-alone Icepack column (models/icepack_pipe/icepack_pipe) used by the MPAS-ADC and LES couplings.
# It links the STIM-Icepack wrapper and Icepack libraries of the GOTM build, so run scripts/build_gotm.sh first.
# Usage: scripts/build_icepack_pipe.sh [gotm_build_dir]   (default: build/gotm); FC defaults to gfortran
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd); B=${1:-$ROOT/build/gotm}; FC=${FC:-gfortran}
cd "$ROOT/models/icepack_pipe"
$FC -O3 -I"$B/extern/stim/modules" -I"$B/extern/stim/icepack_modules" icepack_pipe.F90 \
    "$B/extern/stim/libstim.a" "$B/extern/stim/libstim_icepack_lib.a" -o icepack_pipe
rm -f *.mod
ls -l icepack_pipe
