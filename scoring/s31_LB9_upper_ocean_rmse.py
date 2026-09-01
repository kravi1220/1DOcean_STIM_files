"""
Upper-ocean-only RMSE_T (top 20 m and top 50 m), alongside the existing
full-depth-profile RMSE_T (Section 3.4), for the main-event 9-cell
factorial and persistence baseline. Full-depth RMSE_T is dominated by a
deep-water signal (roughly 400 m warming from 0.7 to 5.8 degC between
Oct 2018 and Aug 2019, Section 4.3) that a one-dimensional column model
cannot represent regardless of turbulence closure or ice-scheme choice,
and that has nothing to do with near-freezing engagement at the surface
-- restricting the comparison to the upper ocean isolates the part of
the water column ice physics can actually act on. Same scoring
methodology as s16_LB9_analysis.py / s22_LB9_significance_check.py
(same variable, same 2 m interpolation grid, same nearest-time
matching), just truncated to the shallower depth range.
"""

from pathlib import Path
_PKG_ROOT = Path(__file__).resolve().parent.parent  # package root (this script lives one level below it)
import netCDF4 as nc
import numpy as np
import csv

CASE = str(_PKG_ROOT / "model_inputs/case_LB9")
CTD_FILE = str(_PKG_ROOT / "data/ctd_profiles_LB9.nc")
TAB_DIR = str(_PKG_ROOT / "tables_output")
VALIDATION_IDX = [4, 5, 6, 7]
DATE_LABELS = ["2018-10-13", "2019-02-14", "2019-05-15", "2019-08-09"]
DEPTH_LIMITS = {"top20": 20.0, "top50": 50.0}

ds = nc.Dataset(CTD_FILE)
depth_ctd = ds.variables["depth"][:]
temp_ctd = ds.variables["temperature"][:]
dates_ctd = nc.num2date(ds.variables["time"][:], ds.variables["time"].units,
                         only_use_cftime_datetimes=False, only_use_python_datetimes=True)
ds.close()


def ctd_profile(idx):
    t = temp_ctd[idx]
    valid = ~np.isnan(t) & ~np.isnan(depth_ctd)
    dv, tv = depth_ctd[valid], t[valid]
    order = np.argsort(dv)
    return dates_ctd[idx], dv[order], tv[order]


def score_run(name, depth_limit):
    dsr = nc.Dataset(f"{CASE}/runs/{name}/all.nc")
    t_time = nc.num2date(dsr.variables["time"][:], dsr.variables["time"].units,
                          only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    z_all = dsr.variables["z"][:]
    temp_all = dsr.variables["temp"][:]
    rmses = []
    for idx in VALIDATION_IDX:
        date, dv_obs, tv_obs = ctd_profile(idx)
        tidx = int(np.argmin(np.abs([(tt - date).total_seconds() for tt in t_time])))
        z_mod = -z_all[tidx, :, 0, 0]
        t_mod = temp_all[tidx, :, 0, 0]
        order = np.argsort(z_mod)
        z_mod, t_mod = z_mod[order], t_mod[order]
        zmax = min(dv_obs.max(), z_mod.max(), depth_limit)
        zcommon = np.arange(2.0, zmax, 2.0)
        t_obs_i = np.interp(zcommon, dv_obs, tv_obs)
        t_mod_i = np.interp(zcommon, z_mod, t_mod)
        rmses.append(float(np.sqrt(np.mean((t_mod_i - t_obs_i) ** 2))))
    dsr.close()
    return np.array(rmses)


def score_persistence(depth_limit):
    _, dv0, tv0 = ctd_profile(3)  # 2018-08-04 IC
    rmses = []
    for idx in VALIDATION_IDX:
        date, dv_obs, tv_obs = ctd_profile(idx)
        zmax = min(dv_obs.max(), dv0.max(), depth_limit)
        zcommon = np.arange(2.0, zmax, 2.0)
        t_obs_i = np.interp(zcommon, dv_obs, tv_obs)
        t_pers_i = np.interp(zcommon, dv0, tv0)
        rmses.append(float(np.sqrt(np.mean((t_pers_i - t_obs_i) ** 2))))
    return np.array(rmses)


RUNS = {"persistence": None, "k_eps_inert": "baseline", "kpp_inert": "turb_kpp",
        "komega_lebedev": "cross_komega_lebedev_lagrangian", "komega_mylake": "cross_komega_mylake_lagrangian"}

for tag, depth_limit in DEPTH_LIMITS.items():
    per_date = {}
    for key, run in RUNS.items():
        per_date[key] = score_persistence(depth_limit) if run is None else score_run(run, depth_limit)
        print(f"[{tag}, <={depth_limit:.0f}m] {key}: per-date={np.round(per_date[key],4)}  mean={per_date[key].mean():.4f}")

    with open(f"{TAB_DIR}/T33_LB9_{tag}_per_date_rmse.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["date"] + list(RUNS.keys()))
        for i, d in enumerate(DATE_LABELS):
            w.writerow([d] + [per_date[k][i] for k in RUNS])
    print(f"T33 ({tag}) written\n")

# aggregate (4-date mean) upper-ocean RMSE_T per closure, for Table 2
agg_rows = []
for closure, run in [("k-eps", "baseline"), ("k-omega", "cross_komega_lebedev_lagrangian"), ("kpp", "turb_kpp")]:
    row = {"closure": closure}
    for tag, depth_limit in DEPTH_LIMITS.items():
        row[tag] = score_run(run, depth_limit).mean()
    agg_rows.append(row)
persist_row = {"closure": "persistence"}
for tag, depth_limit in DEPTH_LIMITS.items():
    persist_row[tag] = score_persistence(depth_limit).mean()
agg_rows.append(persist_row)

with open(f"{TAB_DIR}/T34_LB9_upper_ocean_aggregate.csv", "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["closure", "top20", "top50"])
    for r in agg_rows:
        w.writerow([r["closure"], r["top20"], r["top50"]])
print("T34 written")
for r in agg_rows:
    print(r)
