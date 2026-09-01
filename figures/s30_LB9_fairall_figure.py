"""
Fairall-flux robustness figure (Figure 6): actual simulated near-surface
temperature and ice-thickness trajectories under Kondo (1975) vs. Fairall
et al. (1996) bulk flux formulas, for the same configurations, at LB9.
This is the trajectory-level view of the Table 3 summary counts.
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

RUNS_KONDO = str(_PKG_ROOT / "model_inputs/case_LB9/runs")
RUNS_FAIRALL = str(_PKG_ROOT / "model_inputs/case_LB9/runs_fairall")
CTD_FILE = str(_PKG_ROOT / "data/ctd_profiles_LB9.nc")
FIG_DIR = str(_PKG_ROOT / "figures_output")
VALIDATION_IDX = [4, 5, 6, 7]

# (label, kondo_run, fairall_run, color)
CELLS = [
    ("k-ω+Lebedev", "cross_komega_lebedev_lagrangian", "fa_komega_lebedev", "0.35"),
    ("k-ω+Winton", "turb_komega", "fa_komega_winton", "#d62728"),
    ("k-ω+MyLake", "cross_komega_mylake_lagrangian", "fa_komega_mylake", "#1f77b4"),
]


def load_run(rundir, name, varnames):
    f = nc.Dataset(f"{rundir}/{name}/all.nc")
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

ax = axes[0]
for label, kondo_run, fairall_run, color in CELLS:
    dk = load_run(RUNS_KONDO, kondo_run, ["temp"])
    df = load_run(RUNS_FAIRALL, fairall_run, ["temp"])
    ax.plot(dk["time"], dk["temp"][:, -1], "-", color=color, lw=2.6,
             label=f"{label}, Kondo" if color != "0.35" else label)
    ax.plot(df["time"], df["temp"][:, -1], "--", color=color, lw=2.2,
             label=f"{label}, Fairall" if color != "0.35" else None)
ax.scatter(obs_dates, obs_sst, s=220, marker="*", color="black", edgecolor="white",
           linewidth=1.0, zorder=6, label="Real CTD checkpoints")
ax.set_ylabel("Near-surface temperature (°C)")
ax.set_xlabel("Date")
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
ax.set_title("(a) Temperature: Kondo (solid) vs. Fairall (dashed)")
bold_ticks(ax)
leg1 = ax.legend(fontsize=12.5, loc="upper left", framealpha=0.95, edgecolor=DARK,
                  handlelength=2.4, labelspacing=0.5)
bold_legend_text(leg1)

ax2 = axes[1]
for label, kondo_run, fairall_run, color in CELLS[1:]:
    dk = load_run(RUNS_KONDO, kondo_run, ["Hice"])
    df = load_run(RUNS_FAIRALL, fairall_run, ["Hice"])
    ax2.plot(dk["time"], dk["Hice"], "-", color=color, lw=2.6, label=f"{label}, Kondo")
    ax2.plot(df["time"], df["Hice"], "--", color=color, lw=2.2, label=f"{label}, Fairall")
ax2.set_ylabel("Simulated sea-ice thickness (m)")
ax2.set_xlabel("Date")
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
ax2.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
ax2.set_title("(b) Ice thickness: Kondo (solid) vs. Fairall (dashed)")
bold_ticks(ax2)

# -- OSI-SAF concentration overlay, same LB9 grid point / cadence / masking
#    rule as Figure 3 (this run covers the identical main-event window) --
osisaf_dates, osisaf_conc, osisaf_status = [], [], []
with open(str(_PKG_ROOT / "data/osisaf_lb9_main_event.csv")) as fh:
    for r in csv.DictReader(fh):
        d = datetime.date.fromisoformat(r["date"])
        osisaf_dates.append(datetime.datetime.combine(d, datetime.time(12, 0)))
        osisaf_conc.append(float(r["concentration_pct"]))
        osisaf_status.append(int(r["status_flag"]))
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
leg2 = ax2.legend(lines2 + lines2b, labels2 + labels2b, fontsize=11.5, loc="upper left")
bold_legend_text(leg2)

fig.tight_layout()
fig.savefig(f"{FIG_DIR}/F34_LB9_fairall_experiments.png", dpi=200, metadata={"Software": None})
print("Figure 6 (F34) saved")
