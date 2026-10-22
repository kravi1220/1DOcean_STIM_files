#!/bin/bash
# Build GOTM v6 with STIM (Winton, MyLake, Lebedev, Icepack) and CVMix -> build/gotm.
# Needs: a Fortran compiler (gfortran 12.3 was used), CMake >= 3.10, NetCDF-Fortran (nf-config on PATH).
# Usage: scripts/build_gotm.sh [build_dir]     (default: build/gotm next to the package root)
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd); B=${1:-$ROOT/build/gotm}
cmake -S "$ROOT/models/gotm" -B "$B" -DCMAKE_BUILD_TYPE=Release \
      -DGOTM_USE_STIM=ON -DGOTM_USE_CVMIX=ON -DGOTM_USE_FABM=OFF \
      -DSTIM_WINTON=ON -DSTIM_LEBEDEV=ON -DSTIM_MYLAKE=ON -DSTIM_BASAL_MELT=ON -DSTIM_ICEPACK=ON
cmake --build "$B" -j 8
ls -l "$B/gotm"
