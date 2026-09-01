"""
LB9 synthetic forcing-offset sweep figure: coldest near-surface temperature
reached by each closure's own no-ice trajectory, as a function of the
synthetic offset, with the real local freezing point marked.

Replaces an earlier version of this figure (engaged-cell count vs. k-omega's
own unconverged coldest temperature) after s34_LB9_sweep_convergence_fix.py
found that version's central signal - "2 of 9 cells engage, real ice forms,
at the two coldest offsets" - was itself the same k-omega time-step
non-convergence artifact already diagnosed in Section 4.1.2 for the main
event: Hice = 0.0 m at both offsets once rerun at the converged dt=450s
time step, and k-omega's own coldest excursion moves far warmer
(-2.963degC -> -0.349degC at -3degC). Once every offset's k-omega
trajectory is examined at dt=450s, all three closures respond smoothly
and consistently to the offset (no volatile outlier), and none ever
reaches the local freezing point anywhere in the tested -3 to +3degC
range - a cleaner, fully-converged version of the same "gap is several
degrees" finding, not a different one.

Produces Figure 6 (F32_LB9_sweep_dose_response.png).
"""

from pathlib import Path
_PKG_ROOT = Path(__file__).resolve().parent.parent  # package root (this script lives one level below it)
import sys
sys.path.insert(0, str(_PKG_ROOT))

import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from seaice1d.plotstyle import apply_style, bold_ticks, bold_legend_text, DARK
apply_style()
plt.rcParams.update({
    "font.size": 15, "axes.labelsize": 17, "axes.labelweight": "bold",
    "axes.titlesize": 16, "axes.titleweight": "bold",
    "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 13,
    "lines.linewidth": 3.0, "axes.linewidth": 1.6,
    "xtick.major.width": 1.5, "ytick.major.width": 1.5,
})

TAB_DIR = str(_PKG_ROOT / "tables_output")
FIG_DIR = str(_PKG_ROOT / "figures_output")
FREEZING_PT = -1.8598  # local freezing point at LB9, TEOS-10, in-situ scale at true sample pressure (Section 2.1)

with open(f"{TAB_DIR}/T28_LB9_sweep_summary.csv") as f:
    rows = list(csv.DictReader(f))

offsets = [int(r["offset"]) for r in rows]
coldest_keps = [float(r["coldest_keps_noice"]) for r in rows]
coldest_komega = [float(r["coldest_komega_noice"]) for r in rows]
# KPP's own coldest values are not in T28 (never re-run for this figure
# script); pull them directly from Table 6's own established, already-cited
# numbers (Section 4.4.3), which use the default dt=1800s throughout - KPP,
# like k-epsilon, is already established elsewhere in this paper (Section
# 4.1.2) to be well-converged at the default time step, so no re-run is
# needed here.
coldest_kpp_by_offset = {-3: 0.012, -2: 0.412, -1: 0.820, 0: 1.214, 1: 1.674, 2: 2.144, 3: 2.649}
coldest_kpp = [coldest_kpp_by_offset[o] for o in offsets]

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(offsets, coldest_keps, "o-", color="0.35", label="k-epsilon", markersize=9)
ax.plot(offsets, coldest_komega, "s--", color="#d62728", label="k-omega (converged, dt=450s)", markersize=9)
ax.plot(offsets, coldest_kpp, "^:", color="#1f77b4", label="KPP", markersize=9)

ax.axhline(FREEZING_PT, color="0.3", ls="--", lw=2.2, zorder=2)
ax.text(offsets[0] - 0.15, FREEZING_PT - 0.15, "local freezing point", fontsize=12.5, ha="left",
        color="0.25", style="italic", va="top")

ax.set_xlabel("Synthetic air-temperature/dewpoint offset (°C)")
ax.set_ylabel("Coldest near-surface temperature reached\nby each closure's own no-ice trajectory (°C)")
ax.set_xticks(offsets)
ax.set_xticklabels([f"{o:+d}" if o != 0 else "0\n(real event)" for o in offsets])
ax.set_title("No closure's own trajectory reaches the local freezing point\nat any tested offset, once numerically converged", fontsize=16)
bold_ticks(ax)

leg = ax.legend(fontsize=13, loc="upper left")
bold_legend_text(leg)

fig.tight_layout()
fig.savefig(f"{FIG_DIR}/F32_LB9_sweep_dose_response.png", dpi=200, metadata={"Software": None})
print("Figure 6 (F32) saved")
print("Offsets:", offsets)
print("k-eps:", coldest_keps)
print("k-omega (converged):", coldest_komega)
print("KPP:", coldest_kpp)
