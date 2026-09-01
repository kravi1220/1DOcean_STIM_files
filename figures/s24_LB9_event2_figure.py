"""
LB9 event-2 generalization test figure: compares the engagement pattern
(bit-identical-to-no-ice classification) and ice-thickness behavior at
the main event (14 Feb 2019, -1.835 degC) vs. the second, milder event
(12 Feb 2020, -1.023 degC) at the same station, same 9-cell factorial.
Produces Figure 5 (F31_LB9_event2_comparison.png).
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

TAB_DIR = str(_PKG_ROOT / "tables_output")
FIG_DIR = str(_PKG_ROOT / "figures_output")
RUNS_MAIN = str(_PKG_ROOT / "model_inputs/case_LB9/runs")
RUNS_EVENT2 = str(_PKG_ROOT / "model_inputs/case_LB9/runs_event2")
RUNS_RESTART = str(_PKG_ROOT / "model_inputs/case_LB9/runs_restart2019")
RUNS_RESTART_OCT2018 = str(_PKG_ROOT / "model_inputs/case_LB9/runs_restart_oct2018")
_CLOSURES = ["keps", "komega", "kpp"]
_ICES = ["winton", "lebedev", "mylake"]


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


T19 = read_csv(f"{TAB_DIR}/T19_LB9_noice_comparison.csv")  # event 1
T24 = read_csv(f"{TAB_DIR}/T24_LB9_event2_noice_comparison.csv")  # event 2


def uniq_pairs(rows):
    seen, out = set(), []
    for r in rows:
        key = (r["closure"], r["ice"])
        if key not in seen:
            seen.add(key)
            out.append(r)
    return out


event1_pairs = uniq_pairs(T19)
event1_engaged = sum(1 for r in event1_pairs if r["identical_to_noice"] == "False")
event1_n = len(event1_pairs)

event2_pairs = uniq_pairs(T24)
event2_engaged = sum(1 for r in event2_pairs if r["identical_to_noice"] == "False")
event2_n = len(event2_pairs)

print(f"Event 1: {event1_engaged}/{event1_n} engaged")
print(f"Event 2 (continuous run): {event2_engaged}/{event2_n} engaged")

# -- restart check: same 9-cell factorial + matched no-ice controls,
#    initialized fresh from the real 25 Oct 2019 cast instead of
#    continued from 4 Aug 2018 (see s06h_make_gotm_yaml_LB9_restart2019.py) --
def load_temp(rundir, name):
    f = nc.Dataset(f"{rundir}/{name}/all.nc")
    temp = np.array(f.variables["temp"][:]).squeeze()[:, -1]
    f.close()
    return temp

restart_rows = []
restart_engaged = 0
for c in _CLOSURES:
    T_ni = load_temp(RUNS_RESTART, f"r19_noice_{c}")
    for i in _ICES:
        T = load_temp(RUNS_RESTART, f"r19_{c}_{i}")
        diff = float(np.max(np.abs(T - T_ni)))
        identical = diff <= 1e-6
        if not identical:
            restart_engaged += 1
        restart_rows.append({"run": f"r19_{c}_{i}", "closure": c, "ice": i,
                              "max_diff_vs_noice": diff, "identical_to_noice": identical})
restart_n = len(restart_rows)

with open(f"{TAB_DIR}/T31_LB9_restart2019_noice_comparison.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["run", "closure", "ice", "max_diff_vs_noice", "identical_to_noice"])
    w.writeheader()
    w.writerows(restart_rows)

print(f"Second event, RESTART from 25 Oct 2019: {restart_engaged}/{restart_n} engaged")

# -- the decisive check: same 9-cell factorial + matched no-ice controls,
#    initialized fresh from the real 13 Oct 2018 cast (already one of the
#    main run's own validation dates) instead of continued from 4 Aug
#    2018, integrated only the 124 days to the main event itself --
main_restart_rows = []
main_restart_engaged = 0
for c in _CLOSURES:
    T_ni = load_temp(RUNS_RESTART_OCT2018, f"roct18_noice_{c}")
    for i in _ICES:
        T = load_temp(RUNS_RESTART_OCT2018, f"roct18_{c}_{i}")
        diff = float(np.max(np.abs(T - T_ni)))
        identical = diff <= 1e-6
        if not identical:
            main_restart_engaged += 1
        main_restart_rows.append({"run": f"roct18_{c}_{i}", "closure": c, "ice": i,
                                   "max_diff_vs_noice": diff, "identical_to_noice": identical})
main_restart_n = len(main_restart_rows)

with open(f"{TAB_DIR}/T32_LB9_restart_oct2018_noice_comparison.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["run", "closure", "ice", "max_diff_vs_noice", "identical_to_noice"])
    w.writeheader()
    w.writerows(main_restart_rows)

print(f"Main event, RESTART from 13 Oct 2018: {main_restart_engaged}/{main_restart_n} engaged")

# ------------------------------------------------------------------
# Figure 5: 2-panel. (a) engaged-fraction bar comparison, four bars
# (main event, continuous run; main event, restarted from the real
# 13 Oct 2018 profile; second event, continuous run; second event,
# restarted from the real 25 Oct 2019 profile). (b) ice thickness time
# series for the two configurations engaged at the main event, spanning
# the full extended (continuous) run, so both events are visible on one
# timeline.
# ------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(16, 6.5))

ax = axes[0]
bars = ax.bar(["Main event,\ncontinuous run\n(14 Feb 2019)",
               "Main event,\nrestart from\n13 Oct 2018",
               "Second event,\ncontinuous run\n(12 Feb 2020)",
               "Second event,\nrestart from\n25 Oct 2019"],
              [event1_engaged, main_restart_engaged, event2_engaged, restart_engaged],
              color=["#d62728", "0.6", "#1f77b4", "0.6"],
              edgecolor=DARK, linewidth=1.5, width=0.6)
ax.set_ylim(0, 9)
ax.set_ylabel("Engaged cells (of 9)")
ax.set_title("(a) Engagement fraction by event")
bold_ticks(ax)

ax2 = axes[1]


def load_hice(rundir, name):
    f = nc.Dataset(f"{rundir}/{name}/all.nc")
    t = f.variables["time"]
    times = nc.num2date(t[:], t.units, calendar=t.calendar,
                         only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    hice = np.array(f.variables["Hice"][:]).squeeze()
    f.close()
    return np.array(times), hice

for label, name, color in [("k-ω+Winton", "e2_komega_winton", "#d62728"),
                             ("k-ω+MyLake", "e2_komega_mylake", "#1f77b4")]:
    times, hice = load_hice(RUNS_EVENT2, name)
    ax2.plot(times, hice, "-", color=color, lw=2.4, label=label)

ax2.set_ylim(0, 2.15)
ax2.axvline(np.datetime64("2019-02-14"), color="0.5", ls=":", lw=1.8)
ax2.axvline(np.datetime64("2020-02-12"), color="0.5", ls=":", lw=1.8)
ax2.text(np.datetime64("2019-02-14"), 1.5, " main\n event", fontsize=11, color="0.4", va="top")
ax2.text(np.datetime64("2020-01-20"), 1.5, "2nd\nevent", fontsize=11, color="0.4", va="top", ha="right")
ax2.set_ylabel("Simulated sea-ice thickness (m)")
ax2.set_xlabel("Date")
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
ax2.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax2.set_title("(b) Ice thickness across the continuous run")
bold_ticks(ax2)

# -- OSI-SAF concentration overlay, same masking rule as Figure 3
#    (status_flag in {0, 4}), spliced from the main-event cache (winter 1)
#    and the event-2 fetch (s05i_fetch_osisaf_LB9_event2.py, winter 2) --
osisaf_dates, osisaf_conc, osisaf_status = [], [], []
try:
    with open(str(_PKG_ROOT / "data/osisaf_lb9_main_event.csv")) as fh:
        for r in csv.DictReader(fh):
            d = datetime.date.fromisoformat(r["date"])
            osisaf_dates.append(datetime.datetime.combine(d, datetime.time(12, 0)))
            osisaf_conc.append(float(r["concentration_pct"]))
            osisaf_status.append(int(r["status_flag"]))
except FileNotFoundError:
    pass
with open(str(_PKG_ROOT / "data/osisaf_lb9_event2.csv")) as fh:
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

lines1, labels1 = ax2.get_legend_handles_labels()
lines2, labels2 = ax2b.get_legend_handles_labels()
leg = ax2.legend(lines1 + lines2, labels1 + labels2, fontsize=12.5, loc="upper left")
bold_legend_text(leg)

fig.tight_layout()
fig.savefig(f"{FIG_DIR}/F31_LB9_event2_comparison.png", dpi=200, metadata={"Software": None})
print("Figure 5 (F31) saved")

# write a small summary CSV for the manuscript to read numbers from
with open(f"{TAB_DIR}/T25_LB9_event2_summary.csv", "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["event", "date", "surface_temp_above_freezing", "n_engaged", "n_total"])
    w.writerow(["main_continuous", "2019-02-14", 0.019, event1_engaged, event1_n])
    w.writerow(["main_restart", "2019-02-14", 0.019, main_restart_engaged, main_restart_n])
    w.writerow(["second_continuous", "2020-02-12", 0.838, event2_engaged, event2_n])
    w.writerow(["second_restart", "2020-02-12", 0.838, restart_engaged, restart_n])
print("T25 written")
