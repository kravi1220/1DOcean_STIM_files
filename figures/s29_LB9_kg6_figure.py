"""
KG6 simulated-results figure (Figure 8): near-surface temperature time
series for representative KG6 configurations, mirroring Figure 2's
design for the main LB9 site, plus the real KG6 CTD checkpoints. This is
the trajectory-level view of the KG6 negative-control test that Table 4's
summary counts alone do not show.
"""

from pathlib import Path
_PKG_ROOT = Path(__file__).resolve().parent.parent  # package root (this script lives one level below it)
import sys
sys.path.insert(0, str(_PKG_ROOT))

import numpy as np
import netCDF4 as nc
import csv
import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from seaice1d.plotstyle import apply_style, bold_ticks, bold_legend_text, DARK
apply_style()
plt.rcParams.update({
    "font.size": 15, "axes.labelsize": 17, "axes.labelweight": "bold",
    "axes.titlesize": 16, "axes.titleweight": "bold",
    "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 13,
    "lines.linewidth": 2.6, "axes.linewidth": 1.6,
    "xtick.major.width": 1.5, "ytick.major.width": 1.5,
})

RUNS = str(_PKG_ROOT / "model_inputs/case_KG6/runs")
CTD_FILE = str(_PKG_ROOT / "data/ctd_profiles_KG6.nc")
FIG_DIR = str(_PKG_ROOT / "figures_output")
VALIDATION_IDX = [4, 5, 6, 7]

REPR_RUNS = {
    "k-ε+Lebedev": ("kg6_keps_lebedev", "0.3", "--"),
    "k-ω+Lebedev": ("kg6_komega_lebedev", "0.55", ":"),
    "k-ε+Winton": ("kg6_keps_winton", "#d62728", "-"),
    "k-ω+MyLake": ("kg6_komega_mylake", "#1f77b4", "-"),
}


def load_run(name, varnames):
    f = nc.Dataset(f"{RUNS}/{name}/all.nc")
    t = f.variables["time"]
    times = nc.num2date(t[:], t.units, calendar=t.calendar,
                         only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    out = {"time": np.array(times)}
    for v in varnames:
        out[v] = np.array(f.variables[v][:]).squeeze()
    f.close()
    return out


ds = nc.Dataset(CTD_FILE)
depth_ctd = ds.variables["depth"][:]
temp_ctd = ds.variables["temperature"][:]
dates_ctd = nc.num2date(ds.variables["time"][:], ds.variables["time"].units,
                         only_use_cftime_datetimes=False, only_use_python_datetimes=True)
ds.close()

obs_dates, obs_sst = [], []
for idx in VALIDATION_IDX:
    t_ = temp_ctd[idx]
    valid = np.isfinite(t_) & np.isfinite(depth_ctd)
    order = np.argsort(depth_ctd[valid])
    obs_dates.append(dates_ctd[idx])
    obs_sst.append(t_[valid][order][0])

fig, axes = plt.subplots(1, 2, figsize=(18, 6.8))

# panel (a): near-surface temperature time series
ax = axes[0]
for label, (run, color, ls) in REPR_RUNS.items():
    d = load_run(run, ["temp"])
    near_surface = d["temp"][:, -1]
    ax.plot(d["time"], near_surface, ls, color=color, lw=2.6, label=label)
ax.scatter(obs_dates, obs_sst, s=220, marker="*", color="black", edgecolor="white",
           linewidth=1.0, zorder=6, label="Real CTD checkpoints")
ax.axhline(-1.7290, color="0.5", ls=":", lw=1.6, zorder=1)
ax.text(obs_dates[0], -1.7290, "  local freezing point", fontsize=11, color="0.4", va="bottom", style="italic")
ax.set_ylabel("Near-surface temperature (°C)")
ax.set_xlabel("Date")
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
ax.set_title("(a) Simulated temperature at KG6")
bold_ticks(ax)
leg1 = ax.legend(fontsize=12.5, loc="upper left", framealpha=0.95, edgecolor=DARK,
                  handlelength=2.4, labelspacing=0.5)
bold_legend_text(leg1)

# panel (b): ice thickness for the two engaged configs
ax2 = axes[1]
for label, (run, color, ls) in [("k-ε+Winton", ("kg6_keps_winton", "#d62728", "-")),
                                  ("k-ω+MyLake", ("kg6_komega_mylake", "#1f77b4", "-"))]:
    d = load_run(run, ["Hice"])
    ax2.plot(d["time"], d["Hice"], ls, color=color, lw=2.6, label=label)
ax2.set_ylabel("Simulated sea-ice thickness (m)")
ax2.set_xlabel("Date")
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
ax2.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
ax2.set_title("(b) Simulated ice thickness at KG6")
bold_ticks(ax2)

# -- OSI-SAF concentration overlay at KG6's own grid point, same masking
#    rule as Figure 3/5b (s04c_fetch_osisaf_KG6.py) --
osisaf_dates, osisaf_conc, osisaf_status = [], [], []
with open(str(_PKG_ROOT / "data/osisaf_kg6.csv")) as fh:
    for r in csv.DictReader(fh):
        d = datetime.date.fromisoformat(r["date"])
        osisaf_dates.append(datetime.datetime.combine(d, datetime.time(12, 0)))
        osisaf_conc.append(float(r["concentration_pct"]) if r["concentration_pct"] else None)
        osisaf_status.append(int(r["status_flag"]) if r["status_flag"] else None)
osisaf_conc_masked = np.array([c if (c is not None and st in (0, 4)) else np.nan
                                for c, st in zip(osisaf_conc, osisaf_status)])

ax2b = ax2.twinx()
ax2b.scatter(osisaf_dates, osisaf_conc_masked, color="k", marker="o", s=45,
             label="OSI-SAF concentration", zorder=5)
ax2b.set_ylabel("OSI-SAF ice concentration (%)")
ax2b.set_ylim(0, 100)
bold_ticks(ax2b)

lines2, labels2 = ax2.get_legend_handles_labels()
lines2b, labels2b = ax2b.get_legend_handles_labels()
leg2 = ax2.legend(lines2 + lines2b, labels2 + labels2b, fontsize=12.5, loc="upper right")
bold_legend_text(leg2)

fig.tight_layout()
fig.savefig(f"{FIG_DIR}/F33_LB9_kg6_experiments.png", dpi=200, metadata={"Software": None})
print("Figure 8 (F33) saved")
