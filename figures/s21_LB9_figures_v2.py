"""
Builds Figures 1, 2, 3, and 4 of the LB9 manuscript:
  Figure 1  F24_LB9_map.png                    (LB9/KG6 positions vs. OSI-SAF concentration)
  Figure 2  F23_LB9_all_experiments.png         (near-surface temperature, all representative runs)
  Figure 3  F25_LB9_ice_thickness_osisaf.png    (ice thickness vs. OSI-SAF concentration)
  Figure 4  F28_LB9_closure_vs_ice_effects.png  (closure-only vs. ice-scheme-only RMSE effects)

Shared style (large, bold text/lines) is applied throughout since this
figure set is denser (more series, more panels) than plotstyle.py's
original defaults were tuned for.
"""

from pathlib import Path
_PKG_ROOT = Path(__file__).resolve().parent.parent  # package root (this script lives one level below it)
import sys
sys.path.insert(0, str(_PKG_ROOT))

import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import datetime
import csv
import pickle

from seaice1d.plotstyle import apply_style, bold_ticks, bold_legend_text, DARK

apply_style()
# Bump substantially beyond plotstyle's defaults for this figure set,
# which is denser (more series, more panels) than plotstyle.py's
# original target figures.
plt.rcParams.update({
    "font.size": 15,
    "axes.labelsize": 17,
    "axes.labelweight": "bold",
    "axes.titlesize": 16,
    "axes.titleweight": "bold",
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 13,
    "lines.linewidth": 2.6,
    "axes.linewidth": 1.6,
    "xtick.major.width": 1.5,
    "ytick.major.width": 1.5,
    "xtick.major.size": 6,
    "ytick.major.size": 6,
})

RUNS_DIR = str(_PKG_ROOT / "model_inputs/case_LB9/runs")
FIG_DIR = str(_PKG_ROOT / "figures_output")
CTD_FILE = str(_PKG_ROOT / "data/ctd_profiles_LB9.nc")

LAT, LON = 66.14767, -27.25850
CHECKPOINT_IDX = [4, 5, 6, 7]

REPR_RUNS = {
    "k-ε": ("baseline", "0.35", "-"),
    "KPP": ("turb_kpp", "0.6", "--"),
    "k-ω + Lebedev": ("cross_komega_lebedev_lagrangian", "#2ca02c", ":"),
    "k-ω + Winton": ("turb_komega", "#d62728", "--"),
    "k-ω + MyLake": ("cross_komega_mylake_lagrangian", "#1f77b4", "-"),
}


def load_run(name, varnames):
    f = nc.Dataset(f"{RUNS_DIR}/{name}/all.nc")
    t = f.variables["time"]
    times = nc.num2date(t[:], t.units, calendar=t.calendar,
                         only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    out = {"time": np.array(times)}
    for v in varnames:
        out[v] = np.array(f.variables[v][:]).squeeze()
    f.close()
    return out


# ---------------------------------------------------------------
# load CTD data once
# ---------------------------------------------------------------
ctdf = nc.Dataset(CTD_FILE)
ctd_time_var = ctdf.variables["time"]
ctd_times = nc.num2date(ctd_time_var[:], ctd_time_var.units,
                         calendar=getattr(ctd_time_var, "calendar", "standard"),
                         only_use_cftime_datetimes=False, only_use_python_datetimes=True)
ctd_depth_all = np.array(ctdf.variables["depth"][:])
ctd_temp_all = np.array(ctdf.variables["temperature"][:])
ctd_salt_all = np.array(ctdf.variables["salinity"][:])
ndim_profile = ctd_depth_all.ndim
ctdf.close()


def ctd_cast(idx):
    if ndim_profile == 2:
        d_, t_, s_ = ctd_depth_all[idx, :], ctd_temp_all[idx, :], ctd_salt_all[idx, :]
    else:
        d_, t_, s_ = ctd_depth_all, ctd_temp_all[idx, :], ctd_salt_all[idx, :]
    return d_, t_, s_


# =================================================================
# Figure 1: map. Crop/mask logic uses cached global OSI-SAF arrays
# directly (data/osisaf_map_*.npy - see the README on provenance).
# =================================================================
lats = np.load(str(_PKG_ROOT / "data/osisaf_map_lats.npy"))
lons = np.load(str(_PKG_ROOT / "data/osisaf_map_lons.npy"))
conc_raw = np.load(str(_PKG_ROOT / "data/osisaf_map_conc_20190214.npy"))

region = (lats > 60) & (lats < 72) & (lons > -42) & (lons < -12)
iy, ix = np.where(region)
y0, y1 = iy.min(), iy.max() + 1
x0, x1 = ix.min(), ix.max() + 1
lats_c = lats[y0:y1, x0:x1]
lons_c = lons[y0:y1, x0:x1]
conc_c = np.ma.masked_invalid(conc_raw[y0:y1, x0:x1])

proj = ccrs.PlateCarree()
fig = plt.figure(figsize=(10.5, 8.5))
ax = fig.add_subplot(1, 1, 1, projection=proj)
ax.set_extent([-40, -14, 61, 71.5], crs=proj)
ax.add_feature(cfeature.OCEAN, facecolor="#eef4f8", zorder=0)
ax.add_feature(cfeature.LAND, facecolor="#ddd6c9", zorder=1)
ax.add_feature(cfeature.COASTLINE, linewidth=1.2, edgecolor=DARK, zorder=2)
ax.spines["geo"].set_edgecolor(DARK)
ax.spines["geo"].set_linewidth(1.6)

gl = ax.gridlines(draw_labels=True, linewidth=0.6, color="gray", alpha=0.6, linestyle="--")
gl.top_labels = False
gl.right_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {"size": 14, "weight": "bold", "color": DARK}
gl.ylabel_style = {"size": 14, "weight": "bold", "color": DARK}

pcm = ax.pcolormesh(lons_c, lats_c, conc_c, cmap="Blues", vmin=0, vmax=100,
                     transform=proj, zorder=3, shading="nearest")

KG6_LAT, KG6_LON = 67.58005324675324, -23.92983766233766
ax.scatter([LON], [LAT], s=420, marker="*", color="#ffcc00", edgecolor=DARK,
           linewidth=1.8, zorder=6, transform=proj, label="LB9")
ax.annotate("LB9", (LON, LAT), fontsize=16, fontweight="bold", color=DARK,
            xytext=(10, -22), textcoords="offset points", transform=proj,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=DARK, lw=0.8, alpha=0.92))
ax.scatter([KG6_LON], [KG6_LAT], s=280, marker="o", color="#2ca02c", edgecolor=DARK,
           linewidth=1.8, zorder=6, transform=proj, label="KG6")
ax.annotate("KG6", (KG6_LON, KG6_LAT), fontsize=16, fontweight="bold", color=DARK,
            xytext=(10, 10), textcoords="offset points", transform=proj,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=DARK, lw=0.8, alpha=0.92))

for name, (lo, la) in {"Greenland": (-32, 68.0), "Iceland": (-19.5, 65.0),
                        "Denmark Strait": (-24.5, 64.3)}.items():
    ax.text(lo, la, name, fontsize=14, style="italic", fontweight="bold", color="#3a3327",
            ha="center", transform=proj, zorder=6)

leg = ax.legend(fontsize=13, loc="lower left", framealpha=0.95, edgecolor=DARK)
bold_legend_text(leg)
fig.suptitle("LB9 and KG6 positions, and OSI-SAF sea-ice concentration, 14 Feb 2019", fontsize=17, fontweight="bold", y=0.97)
fig.subplots_adjust(left=0.10, right=0.90, top=0.90, bottom=0.06)

# Cartopy's set_extent forces a fixed data aspect ratio that rarely matches
# the figure's nominal aspect ratio, so matplotlib auto-shrinks the AXES
# (not the figure canvas) to fit -- a colorbar attached via ax=ax before
# this shrink is applied gets sized to the pre-shrink bounding box, which
# is why the colorbar was rendering much taller than the actual map. Fix:
# force a draw first so the real (post-shrink) axes position is known,
# then place the colorbar explicitly to match that exact height.
fig.canvas.draw()
pos = ax.get_position()
cbar_ax = fig.add_axes([pos.x1 + 0.02, pos.y0, 0.02, pos.height])
cbar = fig.colorbar(pcm, cax=cbar_ax, orientation="vertical")
cbar.set_label("OSI-SAF ice concentration (%)", fontsize=16, fontweight="bold")
cbar.ax.tick_params(labelsize=13)

fig.savefig(f"{FIG_DIR}/F24_LB9_map.png", dpi=200, metadata={"Software": None})
plt.close(fig)
print("Figure 1 (F24 map) done")

# =================================================================
# Figure 2: F23, MERGED into a single panel (was 2 side-by-side panels)
# =================================================================
fig, ax = plt.subplots(figsize=(11, 6.5))
for label, (run, color, ls) in REPR_RUNS.items():
    # NOTE: GOTM's own 'sst' output is the ICE-SURFACE SKIN temperature
    # when ice is present (correlation 0.95 with air temperature, ranging
    # to -12C), not the water temperature -- confirmed by direct
    # inspection. The physically meaningful, CTD-comparable quantity is
    # the actual water temperature in the model's shallowest layer,
    # temp[:, -1] (z is ascending, so index -1 is closest to the surface),
    # which stays correctly pinned near the freezing point under ice.
    d = load_run(run, ["temp"])
    near_surface = d["temp"][:, -1]
    ax.plot(d["time"], near_surface, ls, color=color, lw=2.8 if "MyLake" in label or "Winton" in label else 2.2,
            label=label, zorder=4 if color not in ("0.35", "0.6") else 2)

# overlay real CTD near-surface checkpoints
obs_dates, obs_sst = [], []
for idx in CHECKPOINT_IDX:
    d_, t_, s_ = ctd_cast(idx)
    valid = np.isfinite(d_) & np.isfinite(t_)
    order = np.argsort(d_[valid])
    obs_dates.append(ctd_times[idx])
    obs_sst.append(t_[valid][order][0])  # shallowest valid value
ax.scatter(obs_dates, obs_sst, s=260, marker="*", color="black", edgecolor="white",
           linewidth=1.2, zorder=6, label="Real CTD checkpoints")

ax.axhline(-1.8598, color="0.5", ls=":", lw=1.8, zorder=1)
ax.text(obs_dates[0], -1.8598, "  local freezing point", fontsize=12, color="0.4",
        va="bottom", style="italic")

ax.set_ylabel("Near-surface temperature (°C)")
ax.set_xlabel("Date")
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
bold_ticks(ax)
leg = ax.legend(fontsize=12.5, loc="upper left", ncol=1, framealpha=0.95, edgecolor=DARK,
                 handlelength=2.4, labelspacing=0.4)
bold_legend_text(leg)
fig.subplots_adjust(left=0.09, right=0.97, top=0.95, bottom=0.15)
fig.savefig(f"{FIG_DIR}/F23_LB9_all_experiments.png", dpi=200, metadata={"Software": None})
plt.close(fig)
print("Figure 2 (F23 merged single panel) done")

# =================================================================
# Figure 3: F25 ice thickness vs OSI-SAF -- style fix (legend outside)
# =================================================================
# Reads the persistent, reproducible fetch at the surveyed LB9 position
# (s05k_fetch_osisaf_LB9_main_event.py), not the earlier /tmp-cached
# series (sampled at the pre-correction position, Section 2.1).
osisaf_dates, osisaf_conc, osisaf_status = [], [], []
with open(str(_PKG_ROOT / "data/osisaf_lb9_main_event.csv")) as fh:
    for r in csv.DictReader(fh):
        d = datetime.date.fromisoformat(r["date"])
        osisaf_dates.append(datetime.datetime.combine(d, datetime.time(12, 0)))
        osisaf_conc.append(float(r["concentration_pct"]) if r["concentration_pct"] else None)
        osisaf_status.append(int(r["status_flag"]) if r["status_flag"] else None)
osisaf_conc = np.array([c if c is not None else np.nan for c in osisaf_conc])
osisaf_conc_masked = np.array([c if st in (0, 4) else np.nan for c, st in zip(osisaf_conc, osisaf_status)])

ICE_RUNS = {
    "k-ε / KPP": ("baseline", "0.5", "-"),
    "k-ω + Winton": ("turb_komega", "#d62728", "--"),
    "k-ω + MyLake": ("cross_komega_mylake_lagrangian", "#1f77b4", "-"),
    "k-ω + Lebedev": ("cross_komega_lebedev_lagrangian", "#2ca02c", "-"),
}
fig, ax1 = plt.subplots(figsize=(11, 6.2))
for label, (run, color, ls) in ICE_RUNS.items():
    d = load_run(run, ["Hice"])
    ax1.plot(d["time"], d["Hice"], ls, color=color, lw=2.8, label=label)
ax1.set_ylabel("Simulated sea-ice thickness (m)")
ax1.set_xlabel("Date")
ax1.set_ylim(bottom=0)
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
bold_ticks(ax1)

ax2 = ax1.twinx()
ax2.scatter(osisaf_dates, osisaf_conc_masked, color="k", marker="o", s=55,
            label="OSI-SAF concentration", zorder=5)
ax2.set_ylabel("OSI-SAF ice concentration (%)")
ax2.set_ylim(0, 100)
bold_ticks(ax2)

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
leg = ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=12.5, loc="upper left",
                  ncol=1, framealpha=0.95, edgecolor=DARK, handlelength=2.4, labelspacing=0.5)
bold_legend_text(leg)
fig.subplots_adjust(left=0.09, right=0.90, top=0.95, bottom=0.14)
fig.savefig(f"{FIG_DIR}/F25_LB9_ice_thickness_osisaf.png", dpi=200, metadata={"Software": None})
plt.close(fig)
print("Figure 3 (F25) done")

# =================================================================
# Figure 4: closure-vs-ice effects
# =================================================================
closure_effect = {"k-ε": 3.6707070032255453, "KPP": 3.6662765110040265}
ice_effect_komega = {"Winton\n(contaminated)": 3.608724415705572, "Lebedev": 3.436748222571578, "MyLake": 3.542879843226356}

fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))
ax = axes[0]
x = np.arange(2)
w = 0.28
ax.bar(x - w, [closure_effect[k] for k in closure_effect], width=w, label="Winton",
       color="#d62728", edgecolor=DARK, linewidth=1.3)
ax.bar(x, [closure_effect[k] for k in closure_effect], width=w, label="Lebedev",
       color="#2ca02c", edgecolor=DARK, linewidth=1.3)
ax.bar(x + w, [closure_effect[k] for k in closure_effect], width=w, label="MyLake",
       color="#1f77b4", edgecolor=DARK, linewidth=1.3)
ax.set_xticks(x)
ax.set_xticklabels(list(closure_effect.keys()))
ax.set_ylabel("RMSE$_T$ (°C)")
ax.set_title("(a) Closure-scheme effect\n(k-ε vs. KPP, both ice-inert)", fontsize=15)
ax.set_ylim(3.2, 3.8)
bold_ticks(ax)
leg = ax.legend(fontsize=12, loc="upper right")
bold_legend_text(leg)

ax2 = axes[1]
colors_ice = ["#d62728", "#2ca02c", "#1f77b4"]
ax2.bar(list(ice_effect_komega.keys()), list(ice_effect_komega.values()), color=colors_ice,
        edgecolor=DARK, linewidth=1.3)
ax2.set_ylabel("RMSE$_T$ (°C)")
ax2.set_title("(b) Sea-ice-scheme effect\n(closure = k-ω)", fontsize=15)
ax2.set_ylim(3.2, 3.8)
ax2.axhline(3.6707070032255453, color="gray", ls="--", lw=2.0, label="k-ε/KPP reference")
bold_ticks(ax2)
leg2 = ax2.legend(fontsize=12, loc="upper right")
bold_legend_text(leg2)

fig.tight_layout()
fig.savefig(f"{FIG_DIR}/F28_LB9_closure_vs_ice_effects.png", dpi=200, metadata={"Software": None})
plt.close(fig)
print("Figure 4 (F28) done")
