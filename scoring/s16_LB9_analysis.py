
from pathlib import Path
_PKG_ROOT = Path(__file__).resolve().parent.parent  # package root (this script lives one level below it)
#!/usr/bin/env python
"""
STEP 16 - Full 18-cell factorial + 4-member IC ensemble + persistence for
the LB9 single-continuous-run case (4 Aug 2018 - 9 Aug 2019, 370 days).
Unlike the annual 6901911 extension, this is NOT segmented into phases -
scored as one continuous run against all 4 real CTD checkpoints that fall
within the window (13 Oct 2018, 14 Feb 2019, 15 May 2019, 9 Aug 2019 final),
with ONE aggregate RMSE per run (mean over all 4), treating this as a
single configuration scored over its full run, not as separate events.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import netCDF4 as nc
import gsw
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from seaice1d.plotstyle import apply_style, bold_ticks, bold_legend_text, DARK
apply_style()

CASE = Path(str(_PKG_ROOT / "model_inputs/case_LB9"))
CTD_FILE = Path(str(_PKG_ROOT / "data/ctd_profiles_LB9.nc"))
FIG = Path(str(_PKG_ROOT / "figures_output"))
TAB = Path(str(_PKG_ROOT / "tables_output"))

# CTD occupation indices within the run window (3=IC, 4/5/6=validation, 7=final)
ds = nc.Dataset(CTD_FILE)
time_ctd = ds.variables["time"][:]
depth_ctd = ds.variables["depth"][:]
temp_ctd = ds.variables["temperature"][:]
salt_ctd = ds.variables["salinity"][:]
units_ctd = ds.variables["time"].units
dates_ctd = nc.num2date(time_ctd, units_ctd, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
ds.close()

VALIDATION_IDX = [4, 5, 6, 7]  # 2018-10-13, 2019-02-14, 2019-05-15, 2019-08-09 (final)


def ctd_profile(idx):
    t = temp_ctd[idx]; s = salt_ctd[idx]
    valid = ~np.isnan(t) & ~np.isnan(s)
    dv, tv, sv = depth_ctd[valid], t[valid], s[valid]
    order = np.argsort(dv)
    return dates_ctd[idx], dv[order], tv[order], sv[order]


design = pd.read_csv(CASE / "design_matrix.csv")
RUNS = design.run.tolist() + ["ic_ens_pp", "ic_ens_pm", "ic_ens_mp", "ic_ens_mm"]

rows = []
sst_full_by_run = {}
for name in RUNS:
    ncfile = CASE / "runs" / name / "all.nc"
    if not ncfile.exists():
        print("MISSING:", name); continue
    ds = nc.Dataset(ncfile)
    t_time = nc.num2date(ds.variables["time"][:], ds.variables["time"].units,
                          only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    z_all = ds.variables["z"][:]; temp_all = ds.variables["temp"][:]; salt_all = ds.variables["salt"][:]
    sst_full_by_run[name] = dict(time=np.array(t_time), sst=temp_all[:, -1, 0, 0],
                                   hice=np.array(ds.variables["Hice"][:]).squeeze())

    per_cyc = []
    for idx in VALIDATION_IDX:
        date, dv_obs, tv_obs, sv_obs = ctd_profile(idx)
        tidx = int(np.argmin(np.abs([(tt - date).total_seconds() for tt in t_time])))
        z_mod = -z_all[tidx, :, 0, 0]; t_mod = temp_all[tidx, :, 0, 0]; s_mod = salt_all[tidx, :, 0, 0]
        order = np.argsort(z_mod)
        z_mod, t_mod, s_mod = z_mod[order], t_mod[order], s_mod[order]
        zmax = min(dv_obs.max(), z_mod.max())
        zcommon = np.arange(2.0, zmax, 2.0)
        t_obs_i = np.interp(zcommon, dv_obs, tv_obs); s_obs_i = np.interp(zcommon, dv_obs, sv_obs)
        t_mod_i = np.interp(zcommon, z_mod, t_mod); s_mod_i = np.interp(zcommon, z_mod, s_mod)
        rmse_t = float(np.sqrt(np.mean((t_mod_i - t_obs_i) ** 2)))
        rmse_s = float(np.sqrt(np.mean((s_mod_i - s_obs_i) ** 2)))
        per_cyc.append(dict(idx=idx, date=date, rmse_T=rmse_t, rmse_S=rmse_s))
    ds.close()

    dfc = pd.DataFrame(per_cyc)
    rows.append(dict(run=name, rmse_T=dfc.rmse_T.mean(), rmse_S=dfc.rmse_S.mean(),
                      rmse_T_by_cyc=list(dfc.rmse_T), rmse_S_by_cyc=list(dfc.rmse_S)))
    print("scored", name, "mean RMSE_T=%.4f" % dfc.rmse_T.mean())

df = pd.DataFrame(rows)
df = df.merge(design, on="run", how="left")

# ---------------------------------------------------------------------------
# instability classification (same mechanism-based criteria as the annual run)
# ---------------------------------------------------------------------------
def classify(name):
    # LB9 sits in a more weather-variable, largely ice-free-in-summer
    # regime than the EGC site; real storm-driven cooling here can
    # produce single-direction hourly steps exceeding 0.45 degC (verified
    # by direct inspection: two such steps, -0.47 and -0.58 degC, shared
    # identically across every k-eps cell regardless of ice scheme, on
    # the same two real dates - clearly a shared weather event, not
    # per-cell numerical divergence). A magnitude+recurrence threshold
    # alone therefore false-positives here. We instead require genuine
    # SIGN-REVERSING oscillation: a big jump in one direction followed
    # within 10 steps by a big jump back the other way - the actual
    # signature of the k-omega ill-conditioning documented elsewhere in
    # this project - which single-direction weather events cannot produce.
    d = sst_full_by_run[name]
    sst = d["sst"]; hice = d["hice"]
    diffs = np.diff(sst)
    idx_big = np.where(np.abs(diffs) > 0.45)[0]
    chaotic = False
    for i in range(len(idx_big) - 1):
        if idx_big[i + 1] - idx_big[i] <= 10 and np.sign(diffs[idx_big[i]]) != np.sign(diffs[idx_big[i + 1]]):
            chaotic = True
            break
    if not chaotic and len(idx_big) and idx_big[0] <= 2 and abs(diffs[idx_big[0]]) > 1.0:
        chaotic = True  # extreme (>1 degC) jump at the very first step - startup divergence
    n_contam = int(np.sum((hice > 0.05) & (sst > 2.0)))
    contaminated = (not chaotic) and (n_contam > 50)
    return chaotic, contaminated

cls = {name: classify(name) for name in RUNS if name in sst_full_by_run}
df["chaotic"] = df.run.map(lambda r: cls.get(r, (False, False))[0])
df["contaminated"] = df.run.map(lambda r: cls.get(r, (False, False))[1])
df["valid"] = ~df.chaotic & ~df.contaminated
df.drop(columns=["rmse_T_by_cyc", "rmse_S_by_cyc"]).to_csv(TAB / "T18_LB9_skill.csv", index=False)

print(f"\nChaotic ({df.chaotic.sum()}):", df[df.chaotic].run.tolist())
print(f"Contaminated ({df.contaminated.sum()}):", df[df.contaminated].run.tolist())
print(f"Valid ({df.valid.sum()}):", df[df.valid].run.tolist())

# ---------------------------------------------------------------------------
# persistence baseline
# ---------------------------------------------------------------------------
_, dv0, tv0, sv0 = ctd_profile(3)  # 2018-08-04 IC
prows = []
for idx in VALIDATION_IDX:
    date, dv_obs, tv_obs, sv_obs = ctd_profile(idx)
    zmax = min(dv_obs.max(), dv0.max())
    zcommon = np.arange(2.0, zmax, 2.0)
    t_obs_i = np.interp(zcommon, dv_obs, tv_obs); s_obs_i = np.interp(zcommon, dv_obs, sv_obs)
    t_pers_i = np.interp(zcommon, dv0, tv0); s_pers_i = np.interp(zcommon, dv0, sv0)
    prows.append(dict(idx=idx, date=date,
                       rmse_T=float(np.sqrt(np.mean((t_pers_i - t_obs_i) ** 2))),
                       rmse_S=float(np.sqrt(np.mean((s_pers_i - s_obs_i) ** 2)))))
dfp = pd.DataFrame(prows)
dfp.to_csv(TAB / "T18_LB9_persistence.csv", index=False)
print(f"\nPersistence: mean RMSE_T={dfp.rmse_T.mean():.4f}  RMSE_S={dfp.rmse_S.mean():.4f}")

# ---------------------------------------------------------------------------
# main effects (Winton-only for closure; per-level for ice; forcing expected ~0)
# ---------------------------------------------------------------------------
clean = df[df.valid].copy()
balanced = clean[clean.ice == "winton"]
gm_bal = balanced.rmse_T.mean() if len(balanced) else np.nan
main_closure = balanced.groupby("closure").rmse_T.mean() - gm_bal if len(balanced) else pd.Series(dtype=float)
gm = clean.rmse_T.mean()
main_ice = clean.groupby("ice").rmse_T.mean() - gm
main_forcing = clean.groupby("forcing").rmse_T.mean() - gm

range_closure = main_closure.max() - main_closure.min() if len(main_closure) else np.nan
range_ice = main_ice.max() - main_ice.min() if len(main_ice) else np.nan
range_forcing = main_forcing.max() - main_forcing.min() if len(main_forcing) else np.nan

ic_vals = df[df.run.str.startswith("ic_ens")].rmse_T
ic_floor = ic_vals.max() - ic_vals.min() if len(ic_vals) else np.nan

pd.DataFrame([dict(range_closure=range_closure, range_ice=range_ice, range_forcing=range_forcing,
                    ic_noise_floor=ic_floor, n_valid=len(clean), n_chaotic=int(df.chaotic.sum()),
                    n_contaminated=int(df.contaminated.sum()))]).to_csv(TAB / "T18_LB9_effects.csv", index=False)

print(f"\nMain-effect ranges: closure={range_closure:.5f}  ice={range_ice:.5f}  forcing={range_forcing:.5f}  IC_floor={ic_floor:.5f}")
print(f"(forcing main effect is expected ~0 by construction - stationary station, degenerate forcing axis)")

# ---------------------------------------------------------------------------
# figure: full-run baseline SST/SSS/ice timeseries + RMSE by checkpoint
# ---------------------------------------------------------------------------
b = sst_full_by_run["baseline"]
fig, axes = plt.subplots(3, 1, figsize=(10, 8.5), sharex=True)
axes[0].plot(b["time"], b["sst"], "-", color="#2166ac", lw=1.4, label="model (top layer)")
ctd_dates = [dates_ctd[i] for i in VALIDATION_IDX]
ctd_sst = [ctd_profile(i)[2][0] for i in VALIDATION_IDX]
axes[0].plot(ctd_dates, ctd_sst, "o", color=DARK, ms=9, zorder=5, label="CTD (real)")
axes[0].axhline(-1.8598, color="#b2182b", ls=":", lw=1.3, label="local freezing pt")
axes[0].set_ylabel("SST (°C)")
leg = axes[0].legend(fontsize=9); bold_legend_text(leg)

axes[1].plot(b["time"], b["hice"], "-", color="#4d4d4d", lw=1.6)
axes[1].set_ylabel("Ice thickness (m)")

run_rmse = df[df.run == "baseline"].iloc[0]
axes[2].plot(ctd_dates, run_rmse.rmse_T_by_cyc, "o-", color="#d6604d", lw=1.8, ms=7, label="RMSE$_T$")
axes[2].set_ylabel("RMSE$_T$ (°C)")
axes[2].set_xlabel("Date")
bold_ticks(axes[2])
for ax in axes[:2]:
    bold_ticks(ax)

fig.suptitle("LB9 single continuous run (baseline, k-ε+Winton), 4 Aug 2018 - 9 Aug 2019", fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG / "F19_LB9_timeseries.png", dpi=200, bbox_inches="tight", metadata={"Software": None})
print("\nSaved:", FIG / "F19_LB9_timeseries.png")

print("\nDONE")
