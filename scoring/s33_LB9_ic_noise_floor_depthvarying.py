"""
Compute the depth-varying initial-condition ensemble's noise floor
(max-min spread of RMSE_T across the 4 members) separately for each of
the three closures it was run on (k-epsilon, k-omega, KPP + Winton),
using the identical scoring methodology as s16_LB9_analysis.py. This is
the noise floor actually used in Section 4.3's persistence/closure
comparisons - the original constant-offset ensemble (Section 3.2)
cannot perturb stratification and understates true IC sensitivity,
which is why a second, depth-varying ensemble was built
(build_lb9_ic_ensemble_depthvarying.py /
s06j_make_gotm_yaml_LB9_ic_ens_depthvarying.py).
"""

from pathlib import Path
_PKG_ROOT = Path(__file__).resolve().parent.parent  # package root (this script lives one level below it)
import netCDF4 as nc
import numpy as np
import csv

CASE = str(_PKG_ROOT / "model_inputs/case_LB9/runs")
CTD_FILE = str(_PKG_ROOT / "data/ctd_profiles_LB9.nc")
TAB_DIR = str(_PKG_ROOT / "tables_output")
VALIDATION_IDX = [4, 5, 6, 7]
CLOSURES = {"keps": "baseline", "komega": "turb_komega", "kpp": "turb_kpp"}
N_MEMBERS = 4

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


def score_run(name):
    dsr = nc.Dataset(f"{CASE}/{name}/all.nc")
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
        zmax = min(dv_obs.max(), z_mod.max())
        zcommon = np.arange(2.0, zmax, 2.0)
        t_obs_i = np.interp(zcommon, dv_obs, tv_obs)
        t_mod_i = np.interp(zcommon, z_mod, t_mod)
        rmses.append(float(np.sqrt(np.mean((t_mod_i - t_obs_i) ** 2))))
    dsr.close()
    return np.array(rmses).mean()


rows = []
for closure, baseline_name in CLOSURES.items():
    vals = [score_run(f"ic_ens_dv{m}_{closure}") for m in range(1, N_MEMBERS + 1)]
    vals = np.array(vals)
    baseline_val = score_run(baseline_name)
    noise_floor = float(vals.max() - vals.min())
    rows.append(dict(closure=closure, baseline_rmse_T=baseline_val,
                      member_rmse_T="|".join(f"{v:.6f}" for v in vals),
                      noise_floor=noise_floor))
    print(f"{closure}: baseline={baseline_val:.6f}, members={np.round(vals,6)}, "
          f"noise_floor(max-min)={noise_floor:.6f}")

with open(f"{TAB_DIR}/T35_LB9_ic_noise_floor_depthvarying.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["closure", "baseline_rmse_T", "member_rmse_T", "noise_floor"], lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print("\nT35 written")
