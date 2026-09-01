"""
KG6 negative-control analysis: score the 9-cell factorial + 3 no-ice
controls at KG6 (a non-marginal site, coldest occupation 1.95 degC above
the local freezing point) against its own 4 real CTD checkpoints
(15 Oct 2018, 15 Feb 2019, 18 May 2019, 10 Aug 2019), and check whether
any cell ever differs from its own closure's no-ice control. Expected
(if the mechanism in Section 4.1 is right): 0 of 9 cells engage anywhere,
since even k-omega should not approach the freezing point at this site.
"""

from pathlib import Path
_PKG_ROOT = Path(__file__).resolve().parent.parent  # package root (this script lives one level below it)
import netCDF4 as nc
import numpy as np
import csv

RUNS = str(_PKG_ROOT / "model_inputs/case_KG6/runs")
CTD_FILE = str(_PKG_ROOT / "data/ctd_profiles_KG6.nc")
TAB_DIR = str(_PKG_ROOT / "tables_output")
VALIDATION_IDX = [4, 5, 6, 7]  # 2018-10-15, 2019-02-15, 2019-05-18, 2019-08-10

CLOSURES = ["k-eps", "k-omega", "kpp"]
ICES = ["winton", "lebedev", "mylake"]
SLUG = {"k-eps": "keps", "k-omega": "komega", "kpp": "kpp"}

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


def load_full(name):
    f = nc.Dataset(f"{RUNS}/{name}/all.nc")
    t = f.variables["time"]
    times = nc.num2date(t[:], t.units, calendar=t.calendar,
                         only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    z_all = np.array(f.variables["z"][:]).squeeze()
    temp_all = np.array(f.variables["temp"][:]).squeeze()
    hice = np.array(f.variables["Hice"][:]).squeeze() if "Hice" in f.variables else np.zeros(len(times))
    f.close()
    return np.array(times), z_all, temp_all, hice


def score(times, z_all, temp_all):
    rmses = []
    for idx in VALIDATION_IDX:
        date, dv_obs, tv_obs = ctd_profile(idx)
        tidx = int(np.argmin(np.abs([(tt - date).total_seconds() for tt in times])))
        z_mod = -z_all[tidx, :]
        t_mod = temp_all[tidx, :]
        order = np.argsort(z_mod)
        z_mod, t_mod = z_mod[order], t_mod[order]
        zmax = min(dv_obs.max(), z_mod.max())
        zcommon = np.arange(2.0, zmax, 2.0)
        t_obs_i = np.interp(zcommon, dv_obs, tv_obs)
        t_mod_i = np.interp(zcommon, z_mod, t_mod)
        rmses.append(float(np.sqrt(np.mean((t_mod_i - t_obs_i) ** 2))))
    return np.array(rmses)


rows, cache = [], {}
coldest_overall = {}
for turb in CLOSURES:
    for ice in ICES:
        name = f"kg6_{SLUG[turb]}_{ice}"
        times, z_all, temp_all, hice = load_full(name)
        cache[name] = (times, z_all, temp_all, hice)
        rmses = score(times, z_all, temp_all)
        rows.append(dict(run=name, closure=turb, ice=ice, rmse_T=rmses.mean(), max_hice=float(np.nanmax(hice))))
        print(f"{name}: RMSE_T={rmses.mean():.4f}, max_Hice={np.nanmax(hice):.4f}")

noice_cache = {}
for turb in CLOSURES:
    name = f"kg6_noice_{SLUG[turb]}"
    times, z_all, temp_all, hice = load_full(name)
    noice_cache[turb] = (times, z_all, temp_all)
    rmses = score(times, z_all, temp_all)
    rows.append(dict(run=name, closure=turb, ice="no_ice", rmse_T=rmses.mean(), max_hice=0.0))
    coldest_overall[turb] = float(np.min(temp_all[:, -1]))
    print(f"{name} (no-ice control): RMSE_T={rmses.mean():.4f}, coldest near-surface temp={coldest_overall[turb]:.3f}C")

comparison_rows = []
for turb in CLOSURES:
    times_noice, z_noice, temp_noice = noice_cache[turb]
    for ice in ICES:
        name = f"kg6_{SLUG[turb]}_{ice}"
        times, z_all, temp_all, hice = cache[name]
        assert len(times) == len(times_noice)
        max_diff = float(np.nanmax(np.abs(temp_all - temp_noice)))
        comparison_rows.append(dict(run=name, closure=turb, ice=ice,
                                     max_diff_vs_noice=max_diff, identical_to_noice=(max_diff < 1e-6)))
        print(f"{name} vs noice_{turb}: max_diff={max_diff:.6f}  identical={max_diff < 1e-6}")

with open(f"{TAB_DIR}/T29_LB9_kg6_skill.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["run", "closure", "ice", "rmse_T", "max_hice"], lineterminator="\n")
    w.writeheader(); w.writerows(rows)
with open(f"{TAB_DIR}/T30_LB9_kg6_noice_comparison.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["run", "closure", "ice", "max_diff_vs_noice", "identical_to_noice"], lineterminator="\n")
    w.writeheader(); w.writerows(comparison_rows)

n_engaged = sum(1 for r in comparison_rows if not r["identical_to_noice"])
print(f"\nKG6 (non-marginal negative control): {n_engaged} of 9 cells engaged")
print("Coldest no-ice excursion by closure:", coldest_overall)
print("T29 and T30 written.")
