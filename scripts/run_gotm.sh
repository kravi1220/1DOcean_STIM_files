#!/bin/bash
# Run GOTM in every experiment directory (one gotm.yaml each) below a folder.
# Usage: scripts/run_gotm.sh <path/to/gotm> [folder]   (folder defaults to experiments/gotm)
# Output: all.nc (and, for Icepack runs, icepack_column.dat) in each run directory.
# This GOTM build can exit with code 1 after a complete run, so success is judged by all.nc.
set -u
GOTM=$(realpath "$1"); ROOT=${2:-$(dirname "$0")/../experiments/gotm}
fail=0
for y in $(find "$ROOT" -name gotm.yaml | sort); do
  d=$(dirname "$y")
  [ -f "$d/icepack_in" ] && mkdir -p "$d/restart"
  # Block W windows start from a given state (LES window profiles, quiescent turbulence); GOTM overwrites restart.nc
  [ -f "$d/restart_initial.nc" ] && cp "$d/restart_initial.nc" "$d/restart.nc"
  ( cd "$d" && "$GOTM" > gotm.log 2>&1 )
  if [ -s "$d/all.nc" ]; then echo "ok   $d"; else echo "FAIL $d (see $d/gotm.log)"; fail=1; fi
done
exit $fail
