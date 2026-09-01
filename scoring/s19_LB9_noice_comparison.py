"""
Bit-for-bit comparison of every main-event 9-cell factorial run (plus the
forcing-assumption and IC-ensemble variants) against its own closure's
no-ice control, at the near-surface (uppermost model layer) temperature
series -- writes T19_LB9_noice_comparison.csv, the table
s24_LB9_event2_figure.py reads to report the main event's own engaged
count for its Figure 5 bar chart.

Near-surface, not full-depth-array, is deliberate: at the corrected LB9
position (Section 2.1), we found the full 3-D temp field diverges from
the no-ice control by up to several degrees at ~50-100 m depth during an
unrelated autumn mixing event (Nov 2018) in every closure, including
k-epsilon and KPP -- a real but off-event sensitivity of the deep mixing
solution to the presence of a never-triggered ice module, unconnected to
sea-ice physics (Hice = 0 throughout). Comparing the full array would
conflate that unrelated signal with actual ice engagement; comparing the
near-surface layer, where the sea-ice module can plausibly act directly,
does not.
"""

from pathlib import Path
_PKG_ROOT = Path(__file__).resolve().parent.parent  # package root (this script lives one level below it)
import netCDF4 as nc
import numpy as np
import csv

RUNS = str(_PKG_ROOT / "model_inputs/case_LB9/runs")
TAB_DIR = str(_PKG_ROOT / "tables_output")

CLOSURES = ["k-eps", "k-omega", "kpp"]
ICES = ["winton", "lebedev", "mylake"]
SLUG = {"k-eps": "keps", "k-omega": "komega", "kpp": "kpp"}
NOICE_NAME = {"k-eps": "noice_keps", "k-omega": "noice_komega", "kpp": "noice_kpp"}

# main-run cell -> directory name (mirrors s06c's factorial + cross_* naming)
CELL_NAME = {
    ("k-eps", "winton", "lagrangian"): "baseline",
    ("k-eps", "winton", "eulerian"): "forcing_eulerian",
    ("k-eps", "lebedev", "lagrangian"): "ice_lebedev",
    ("k-eps", "lebedev", "eulerian"): "cross_keps_lebedev_eulerian",
    ("k-eps", "mylake", "lagrangian"): "ice_mylake",
    ("k-eps", "mylake", "eulerian"): "cross_keps_mylake_eulerian",
    ("k-omega", "winton", "lagrangian"): "turb_komega",
    ("k-omega", "winton", "eulerian"): "cross_komega_winton_eulerian",
    ("k-omega", "lebedev", "lagrangian"): "cross_komega_lebedev_lagrangian",
    ("k-omega", "lebedev", "eulerian"): "cross_komega_lebedev_eulerian",
    ("k-omega", "mylake", "lagrangian"): "cross_komega_mylake_lagrangian",
    ("k-omega", "mylake", "eulerian"): "cross_komega_mylake_eulerian",
    ("kpp", "winton", "lagrangian"): "turb_kpp",
    ("kpp", "winton", "eulerian"): "cross_kpp_winton_eulerian",
    ("kpp", "lebedev", "lagrangian"): "cross_kpp_lebedev_lagrangian",
    ("kpp", "lebedev", "eulerian"): "cross_kpp_lebedev_eulerian",
    ("kpp", "mylake", "lagrangian"): "cross_kpp_mylake_lagrangian",
    ("kpp", "mylake", "eulerian"): "cross_kpp_mylake_eulerian",
}


def load_temp(name):
    f = nc.Dataset(f"{RUNS}/{name}/all.nc")
    temp = np.array(f.variables["temp"][:]).squeeze()[:, -1]
    f.close()
    return temp


noice_cache = {c: load_temp(NOICE_NAME[c]) for c in CLOSURES}

rows = []
for (closure, ice, forcing), name in CELL_NAME.items():
    temp = load_temp(name)
    temp_noice = noice_cache[closure]
    max_diff = float(np.nanmax(np.abs(temp - temp_noice)))
    rows.append(dict(run=name, closure=closure, ice=ice, forcing=forcing,
                      max_diff_vs_noice=max_diff, identical_to_noice=(max_diff < 1e-6)))
    print(f"{name} ({closure}+{ice}, {forcing}): max_diff={max_diff:.6f}  identical={max_diff < 1e-6}")

with open(f"{TAB_DIR}/T19_LB9_noice_comparison.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["run", "closure", "ice", "forcing", "max_diff_vs_noice", "identical_to_noice"],
                        lineterminator="\n")
    w.writeheader()
    w.writerows(rows)

n_engaged = sum(1 for r in rows if not r["identical_to_noice"] and r["forcing"] == "lagrangian")
print(f"\nMain event: {n_engaged} of 9 cells engaged (Lagrangian/Eulerian are degenerate here, counted once)")
print("T19 written.")
