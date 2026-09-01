"""
Writes T20_LB9_significant_differences.csv: full-depth RMSE_T/RMSE_S for
each turbulence closure (representative cell, since no ice scheme ever
differs from any other within a closure at the corrected LB9 position --
Hice = 0 in every one of the 9 cells) against the persistence forecast,
with a significance flag based on the 4-member initial-condition
ensemble's noise floor (5x the ensemble spread). Feeds Table 1 of the
manuscript.

Computed live from T18_LB9_skill.csv and T18_LB9_persistence.csv (both
written by s16_LB9_analysis.py from the current model output) rather than
from hardcoded values, so this table is reproducible from a rerun.
"""

from pathlib import Path
_PKG_ROOT = Path(__file__).resolve().parent.parent  # package root (this script lives one level below it)
import csv
import statistics

TAB_DIR = str(_PKG_ROOT / "tables_output")

with open(f"{TAB_DIR}/T18_LB9_skill.csv") as fh:
    skill = list(csv.DictReader(fh))
with open(f"{TAB_DIR}/T18_LB9_persistence.csv") as fh:
    persist_rows = list(csv.DictReader(fh))

persist_rmse_T = statistics.mean(float(r["rmse_T"]) for r in persist_rows)
persist_rmse_S = statistics.mean(float(r["rmse_S"]) for r in persist_rows)

# IC-ensemble noise floor: spread (max-min) of the 4-member ensemble,
# layered on the k-eps+Winton cell (Section 3.2)
ic_rows = [r for r in skill if r["run"].startswith("ic_ens_")]
ic_rmse_T = [float(r["rmse_T"]) for r in ic_rows]
ic_noise_floor = max(ic_rmse_T) - min(ic_rmse_T)
sig_threshold = 5 * ic_noise_floor

# one representative cell per closure: every ice scheme within a closure
# gives bit-identical RMSE_T/RMSE_S here (Hice = 0 throughout), confirmed
# directly against T18's own rows before picking baseline/turb_komega/turb_kpp
CLOSURE_REP = {"k-eps": "baseline", "k-omega": "turb_komega", "kpp": "turb_kpp"}
by_run = {r["run"]: r for r in skill}

rows = []
for closure, rep in CLOSURE_REP.items():
    r = by_run[rep]
    rmse_T, rmse_S = float(r["rmse_T"]), float(r["rmse_S"])
    delta = rmse_T - persist_rmse_T
    sig = "significant" if abs(delta) > sig_threshold else "not significant (within IC noise floor)"
    rows.append((f"{closure} (all 3 ice schemes; bit-identical; Hice=0 throughout)",
                 rmse_T, rmse_S, delta, f"vs. persistence: {sig}"))
rows.append(("Persistence (null baseline)", persist_rmse_T, persist_rmse_S, 0.0, "reference"))

with open(f"{TAB_DIR}/T20_LB9_significant_differences.csv", "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["configuration", "rmse_T", "rmse_S", "delta_vs_persistence", "note"])
    for r in rows:
        w.writerow(r)

print("T20 written")
for r in rows:
    print(r)

print(f"\nIC noise floor (ensemble spread): {ic_noise_floor:.6f} degC")
print(f"Significance threshold used: {sig_threshold:.6f} degC (5x IC noise floor)")
keps_kpp_diff = abs(by_run["baseline"]["rmse_T"] and float(by_run["baseline"]["rmse_T"]) - float(by_run["turb_kpp"]["rmse_T"]))
print(f"k-eps vs KPP closure-only difference: {keps_kpp_diff:.6f} -> "
      f"{'significant' if keps_kpp_diff > sig_threshold else 'not significant (below threshold)'}")
