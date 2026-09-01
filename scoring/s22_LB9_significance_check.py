"""
Per-date RMSE breakdown and formal significance checks for the LB9
factorial: the aggregate "fractionally worse than persistence" and
"closure vs. ice-scheme effect" claims (Section 4.3) need per-date and
significance support, not just the aggregate number. Reuses the exact
scoring methodology from s16_LB9_analysis.py (same variable ['temp'],
same zcommon = np.arange(2.0, zmax, 2.0) grid, same nearest-time
matching) so results are directly comparable to T18_LB9_skill.csv /
T18_LB9_persistence.csv.
"""

from pathlib import Path
_PKG_ROOT = Path(__file__).resolve().parent.parent  # package root (this script lives one level below it)
import netCDF4 as nc
import numpy as np
import csv
from scipy import stats

CASE = str(_PKG_ROOT / "model_inputs/case_LB9")
CTD_FILE = str(_PKG_ROOT / "data/ctd_profiles_LB9.nc")
TAB_DIR = str(_PKG_ROOT / "tables_output")
VALIDATION_IDX = [4, 5, 6, 7]
DATE_LABELS = ["2018-10-13", "2019-02-14", "2019-05-15", "2019-08-09"]

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
    """Per-date RMSE_T, identical methodology to s16_LB9_analysis.py."""
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
        zmax = min(dv_obs.max(), z_mod.max())
        zcommon = np.arange(2.0, zmax, 2.0)
        t_obs_i = np.interp(zcommon, dv_obs, tv_obs)
        t_mod_i = np.interp(zcommon, z_mod, t_mod)
        rmses.append(float(np.sqrt(np.mean((t_mod_i - t_obs_i) ** 2))))
    dsr.close()
    return np.array(rmses)


def score_persistence():
    _, dv0, tv0 = ctd_profile(3)  # 2018-08-04 IC
    rmses = []
    for idx in VALIDATION_IDX:
        date, dv_obs, tv_obs = ctd_profile(idx)
        zmax = min(dv_obs.max(), dv0.max())
        zcommon = np.arange(2.0, zmax, 2.0)
        t_obs_i = np.interp(zcommon, dv_obs, tv_obs)
        t_pers_i = np.interp(zcommon, dv0, tv0)
        rmses.append(float(np.sqrt(np.mean((t_pers_i - t_obs_i) ** 2))))
    return np.array(rmses)


rmse_persist = score_persistence()
rmse_keps = score_run("baseline")
rmse_kpp = score_run("turb_kpp")
rmse_lebedev = score_run("cross_komega_lebedev_lagrangian")
rmse_mylake = score_run("cross_komega_mylake_lagrangian")

# sanity check against the existing pipeline's aggregate numbers
assert abs(rmse_persist.mean() - 3.4363) < 0.001, rmse_persist.mean()
assert abs(rmse_lebedev.mean() - 3.3279364146999013) < 1e-6, rmse_lebedev.mean()
print("Aggregate sanity check passed:", rmse_persist.mean(), rmse_lebedev.mean())

# ---------------------------------------------------------------
# T21: per-date RMSE_T table
# ---------------------------------------------------------------
with open(f"{TAB_DIR}/T21_LB9_per_date_rmse.csv", "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["date", "persistence", "k_eps_inert", "kpp_inert", "komega_lebedev", "komega_mylake"])
    for i, d in enumerate(DATE_LABELS):
        w.writerow([d, rmse_persist[i], rmse_keps[i], rmse_kpp[i], rmse_lebedev[i], rmse_mylake[i]])

# ---------------------------------------------------------------
# T22: paired significance tests
# ---------------------------------------------------------------
def paired_test(a, b, label):
    diff = a - b
    t_stat, t_p = stats.ttest_1samp(diff, 0)
    try:
        w_stat, w_p = stats.wilcoxon(diff)
    except ValueError:
        w_stat, w_p = np.nan, np.nan
    sign_consistent = bool(np.all(diff > 0) or np.all(diff < 0))
    return dict(comparison=label, mean_diff=diff.mean(), diff_by_date=list(diff),
                t_stat=t_stat, t_pvalue=t_p, wilcoxon_pvalue=w_p, sign_consistent=sign_consistent)

tests = [
    paired_test(rmse_lebedev, rmse_persist, "best_valid_vs_persistence"),
    paired_test(rmse_keps, rmse_kpp, "closure_only_keps_vs_kpp"),
    paired_test(rmse_mylake, rmse_lebedev, "ice_scheme_only_mylake_vs_lebedev_within_komega"),
]

with open(f"{TAB_DIR}/T22_LB9_significance_tests.csv", "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["comparison", "mean_diff", "diff_by_date", "t_stat", "t_pvalue", "wilcoxon_pvalue", "sign_consistent"])
    for r in tests:
        w.writerow([r["comparison"], r["mean_diff"],
                    "|".join(f"{x:+.4f}" for x in r["diff_by_date"]),
                    r["t_stat"], r["t_pvalue"], r["wilcoxon_pvalue"], r["sign_consistent"]])

print("\nT21 and T22 written.")
for r in tests:
    print(r["comparison"], "mean_diff=%.4f" % r["mean_diff"], "t_p=%.4f" % r["t_pvalue"],
          "wilcoxon_p=%.4f" % r["wilcoxon_pvalue"], "sign_consistent=", r["sign_consistent"])
    print("  per-date diffs:", ["%+.4f" % x for x in r["diff_by_date"]])
