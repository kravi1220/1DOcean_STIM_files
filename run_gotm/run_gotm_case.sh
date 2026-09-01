#!/bin/bash
# Run every generated GOTM config in a case's run directory: the
# main-event factorial, both generalization events, the Fairall-flux
# test, each forcing-sweep offset, every numerical-convergence check, and
# KG6 all use this same pattern (one gotm.yaml per immediate subdirectory).
#
# NOTE: this GOTM build exits with code 1 even on a fully successful run
# (confirmed by direct inspection: all.nc is written correctly despite
# the nonzero exit code) -- do not use `set -e`; success is judged by
# the presence of all.nc, not the exit code.
#
# Usage:
#   ./run_gotm_case.sh <runs_dir> [gotm_binary]
#
# Examples:
#   ./run_gotm_case.sh ../model_inputs/case_LB9/runs
#   ./run_gotm_case.sh ../model_inputs/case_LB9/runs_event2
#   ./run_gotm_case.sh ../model_inputs/case_LB9/runs_fairall
#   ./run_gotm_case.sh ../model_inputs/case_LB9/runs_sweep_m1
#   ./run_gotm_case.sh ../model_inputs/case_KG6/runs

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS="$1"
GOTM="${2:-$SCRIPT_DIR/../gotm_source/build/gotm}"

if [ -z "$RUNS" ] || [ ! -d "$RUNS" ]; then
  echo "Usage: $0 <runs_dir> [gotm_binary]"
  exit 1
fi

for d in "$RUNS"/*/; do
  d="${d%/}"
  name=$(basename "$d")
  [ -f "$d/gotm.yaml" ] || continue
  echo "=== running $name ==="
  (cd "$d" && "$GOTM" > gotm.log 2>&1)
  if [ -f "$d/all.nc" ]; then
    echo "  OK: $d/all.nc"
  else
    echo "  FAILED: check $d/gotm.log"
    tail -20 "$d/gotm.log"
  fi
done
echo "ALL RUNS IN $RUNS DONE"
