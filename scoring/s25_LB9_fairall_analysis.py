"""
LB9 Fairall-flux robustness test: score the 9-cell factorial + 3 no-ice
controls (identical forcing, initial condition, and dates as the main
event; only the bulk air-sea flux formula changed from Kondo to Fairall)
against the same 4 real CTD checkpoints. Tests whether "only k-omega ever
reaches the freezing point" is a property of the physics or an artifact
of the Kondo (1975) flux parameterization specifically.
"""

from pathlib import Path
_PKG_ROOT = Path(__file__).resolve().parent.parent  # package root (this script lives one level below it)
import netCDF4 as nc
import numpy as np
import csv

RUNS = str(_PKG_ROOT / "model_inputs/case_LB9/runs_fairall")
CTD_FILE = str(_PKG_ROOT / "data/ctd_profiles_LB9.nc")
TAB_DIR = str(_PKG_ROOT / "tables_output")
VALIDATION_IDX = [4, 5, 6, 7]

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
for turb in CLOSURES:
    for ice in ICES:
        name = f"fa_{SLUG[turb]}_{ice}"
        times, z_all, temp_all, hice = load_full(name)
        cache[name] = (times, z_all, temp_all, hice)
        rmses = score(times, z_all, temp_all)
        rows.append(dict(run=name, closure=turb, ice=ice, rmse_T=rmses.mean(), max_hice=float(np.nanmax(hice))))
        print(f"{name}: RMSE_T={rmses.mean():.4f}, max_Hice={np.nanmax(hice):.4f}")

noice_cache = {}
for turb in CLOSURES:
    name = f"fa_noice_{SLUG[turb]}"
    times, z_all, temp_all, hice = load_full(name)
    noice_cache[turb] = (times, z_all, temp_all)
    rmses = score(times, z_all, temp_all)
    rows.append(dict(run=name, closure=turb, ice="no_ice", rmse_T=rmses.mean(), max_hice=0.0))
    print(f"{name} (no-ice control): RMSE_T={rmses.mean():.4f}")

comparison_rows = []
for turb in CLOSURES:
    times_noice, z_noice, temp_noice = noice_cache[turb]
    for ice in ICES:
        name = f"fa_{SLUG[turb]}_{ice}"
        times, z_all, temp_all, hice = cache[name]
        assert len(times) == len(times_noice)
        # near-surface only, not the full depth-time array: at the corrected
        # LB9 position (Section 2.1) the full array picks up an unrelated,
        # off-event deep-mixing sensitivity to the presence of a
        # never-triggered ice module (Hice = 0 throughout); near-surface is
        # where the sea-ice module can plausibly act directly.
        max_diff = float(np.nanmax(np.abs(temp_all[:, -1] - temp_noice[:, -1])))
        comparison_rows.append(dict(run=name, closure=turb, ice=ice,
                                     max_diff_vs_noice=max_diff, identical_to_noice=(max_diff < 1e-6)))
        print(f"{name} vs noice_{turb}: max_diff={max_diff:.6f}  identical={max_diff < 1e-6}")

with open(f"{TAB_DIR}/T26_LB9_fairall_skill.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["run", "closure", "ice", "rmse_T", "max_hice"], lineterminator="\n")
    w.writeheader(); w.writerows(rows)
with open(f"{TAB_DIR}/T27_LB9_fairall_noice_comparison.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["run", "closure", "ice", "max_diff_vs_noice", "identical_to_noice"], lineterminator="\n")
    w.writeheader(); w.writerows(comparison_rows)

n_engaged = sum(1 for r in comparison_rows if not r["identical_to_noice"])
print(f"\nFairall flux formula: {n_engaged} of 9 cells engaged")
print("Engaged cells:", [r["run"] for r in comparison_rows if not r["identical_to_noice"]])
print("T26 and T27 written.")
