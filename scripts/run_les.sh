#!/bin/bash
# Oceananigans LES runs of the paper (one NVIDIA GPU each; Julia environment in models/les).
# Output (profiles, slices, flux and ice logs, checkpoints) is written to experiments/les/.
# Wall times on one A100: windows 1-3 h each, Oct-Feb run about 20 h, each full-year run about 2.5 days.
# Jobs longer than a queue limit continue from their latest checkpoint ('pickup' / checkpoint argument);
# the segments of the published runs are listed below.
#
# Usage: scripts/run_les.sh <step>     steps: windows | oct15 | alone | coupled
#   Set ICE_ENV to a shell file that makes python3 (numpy, netCDF4) and the compiler runtime of icepack_pipe available
#   to the Icepack server of the coupled run.
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd); M="$ROOT/models/les"
cd "$ROOT/experiments/les"
J="julia --project=$M"
case "$1" in
  windows)   # Block B/W windows: B1 at 1 (24 h), 2 and 4 m; B2 at 2 and 4 m (120 h)
    $J $M/les_kg6_window.jl B1 24 1.0
    for w in B1 B2; do for dx in 2.0 4.0; do $J $M/les_kg6_window.jl $w 120 $dx; done; done ;;
  oct15)     # 15 Oct 2018 cast to 16 Feb 2019 (2958 h), 2 m, stops at the freezing point
    $J $M/les_kg6_era5.jl les_kg6_oct15.nc 2958 2.0 oct15_2m ;;
  alone)     # full year, ocean only, 2 m: 6 Aug 2018 12:00 -> 10 Aug 2019 12:00 (8856 h), past the freezing point.
             # The published run was interrupted at 3600 h and continued from that checkpoint with
             #   $J $M/les_kg6_era5.jl les_kg6_aug06.nc 8856 2.0 aug06_2m pickup nofreezestop
    $J $M/les_kg6_era5.jl les_kg6_aug06.nc 8856 2.0 aug06_2m fresh nofreezestop ;;
  coupled)   # full year, LES x Icepack (mushy), 2 m, in four segments that start and end ice-free
             # (Icepack restarts without ice; ice_cpl2m_<seg>/icepack_in starts at the segment start):
             #   seg1  0    -> 3408 h (26 Dec 2018 12:00)   seg1b 3408 -> 4633 h (15 Feb 2019 13:00)
             #   seg2  4633 -> 7153 h (31 May 2019 13:00)   seg3  7153 -> 8856 h (10 Aug 2019 12:00)
    last() { ls les_kg6_icepack_$1_chk_iteration*.jld2 | sed 's/.*iteration\([0-9]*\).jld2/\1 &/' | sort -n | tail -1 | cut -d' ' -f2; }
    $J $M/les_kg6_era5_icepack.jl les_kg6_aug06.nc 3408 2.0 cpl2m_seg1  none                $PWD/ice_cpl2m_seg1
    $J $M/les_kg6_era5_icepack.jl les_kg6_aug06.nc 4633 2.0 cpl2m_seg1b "$(last cpl2m_seg1)"  $PWD/ice_cpl2m_seg1b
    $J $M/les_kg6_era5_icepack.jl les_kg6_aug06.nc 7153 2.0 cpl2m_seg2  "$(last cpl2m_seg1b)" $PWD/ice_cpl2m_seg2
    $J $M/les_kg6_era5_icepack.jl les_kg6_aug06.nc 8856 2.0 cpl2m_seg3  "$(last cpl2m_seg2)"  $PWD/ice_cpl2m_seg3 ;;
  *) echo "usage: $0 windows|oct15|alone|coupled"; exit 1 ;;
esac
