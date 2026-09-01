"""
Builds the new Figure 5 of the LB9 manuscript (F35_LB9_deep_profile_comparison.png):
depth profiles of real CTD temperature vs. modeled temperature at the two
validation dates bracketing the deep-water warming cited in Section 4.3
(13 Oct 2018 and 9 Aug 2019), full water column. Directly supports the
"~0.7 to 5.8 degC at 400 m" claim with a figure rather than the
unsupported "comparison with CTD observations (not shown here)" the
claim previously rested on.
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

from seaice1d.plotstyle import apply_style, bold_ticks, bold_legend_text, DARK

apply_style()
plt.rcParams.update({
    "font.size": 15, "axes.labelsize": 17, "axes.labelweight": "bold",
    "axes.titlesize": 16, "axes.titleweight": "bold",
    "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 13,
    "lines.linewidth": 2.6, "axes.linewidth": 1.6,
})

FIG_DIR = str(_PKG_ROOT / "figures_output")
CTD_FILE = str(_PKG_ROOT / "data/ctd_profiles_LB9.nc")
RUN = str(_PKG_ROOT / "model_inputs/case_LB9/runs/baseline")
VALIDATION = [(4, "13 Oct 2018"), (7, "9 Aug 2019")]

ds = nc.Dataset(CTD_FILE)
depth_ctd = ds.variables["depth"][:]
temp_ctd = ds.variables["temperature"][:]
dates_ctd = nc.num2date(ds.variables["time"][:], ds.variables["time"].units,
                         only_use_cftime_datetimes=False, only_use_python_datetimes=True)
ds.close()

dsr = nc.Dataset(f"{RUN}/all.nc")
t_time = nc.num2date(dsr.variables["time"][:], dsr.variables["time"].units,
                      only_use_cftime_datetimes=False, only_use_python_datetimes=True)
z_all = dsr.variables["z"][:]
temp_all = dsr.variables["temp"][:]
dsr.close()

fig, axes = plt.subplots(1, 2, figsize=(11, 8), sharey=True)

for ax, (idx, label) in zip(axes, VALIDATION):
    t_obs = temp_ctd[idx]
    valid = ~np.isnan(t_obs) & ~np.isnan(depth_ctd)
    dv_obs, tv_obs = depth_ctd[valid], t_obs[valid]
    order = np.argsort(dv_obs)
    dv_obs, tv_obs = dv_obs[order], tv_obs[order]

    date = dates_ctd[idx]
    tidx = int(np.argmin(np.abs([(tt - date).total_seconds() for tt in t_time])))
    z_mod = -z_all[tidx, :, 0, 0]
    t_mod = temp_all[tidx, :, 0, 0]
    order = np.argsort(z_mod)
    z_mod, t_mod = z_mod[order], t_mod[order]

    ax.plot(tv_obs, dv_obs, "-", color=DARK, lw=3.0, label="Real CTD cast")
    ax.plot(t_mod, z_mod, "--", color="#d62728", lw=2.6, label="Model (k-eps+Winton)")
    ax.axhline(400, color="0.6", ls=":", lw=1.6)
    ax.set_title(label)
    ax.set_xlabel("Temperature (°C)")
    bold_ticks(ax)

# Called ONCE, not per-axis: the two subplots share their y-axis (sharey=True),
# so calling invert_yaxis() inside the loop above inverted it twice - once per
# axis - which cancels out and silently leaves the axis non-inverted (0 m at
# the bottom, 500 m at the top - backwards from the standard oceanographic
# convention of shallow-at-top/deep-at-bottom this figure is meant to use).
axes[0].invert_yaxis()

axes[0].set_ylabel("Depth (m)")
axes[0].text(0.5, 405, "400 m", fontsize=11, color="0.4", va="top")
leg = axes[1].legend(loc="lower right", fontsize=13)
bold_legend_text(leg)
fig.suptitle("Real vs. modeled temperature profile, before and after the deep-water warming", fontsize=16, fontweight="bold", y=0.99)
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/F35_LB9_deep_profile_comparison.png", dpi=200, metadata={"Software": None})
print("Figure (F35) saved")

t400_obs_1 = np.interp(400, depth_ctd[~np.isnan(temp_ctd[4]) & ~np.isnan(depth_ctd)],
                        temp_ctd[4][~np.isnan(temp_ctd[4]) & ~np.isnan(depth_ctd)])
t400_obs_2 = np.interp(400, depth_ctd[~np.isnan(temp_ctd[7]) & ~np.isnan(depth_ctd)],
                        temp_ctd[7][~np.isnan(temp_ctd[7]) & ~np.isnan(depth_ctd)])
print(f"Real CTD at 400 m: {t400_obs_1:.3f} C (13 Oct 2018) -> {t400_obs_2:.3f} C (9 Aug 2019)")
