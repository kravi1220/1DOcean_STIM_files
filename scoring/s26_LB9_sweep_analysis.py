"""
LB9 synthetic forcing-offset sweep analysis: for each offset level
(-3,-2,-1,0[real main event],+1,+2,+3 degC applied to air temperature and
dewpoint), compute (a) how many of the 9 closure x ice-scheme cells
engage (differ bit-for-bit from their own closure's no-ice control), and
(b) a physically meaningful "distance from freezing" metric -- the
coldest near-surface water temperature reached by the k-epsilon no-ice
trajectory (chosen because k-epsilon never engages ice at any offset
tested, so it cleanly tracks how cold the column gets under each forcing
level without ice thermodynamics muddying the signal). Produces a
dose-response curve: engaged count vs. distance from freezing.
"""

from pathlib import Path
_PKG_ROOT = Path(__file__).resolve().parent.parent  # package root (this script lives one level below it)
import netCDF4 as nc
import numpy as np
import csv

CASE = str(_PKG_ROOT / "model_inputs/case_LB9")
TAB_DIR = str(_PKG_ROOT / "tables_output")

CLOSURES = ["k-eps", "k-omega", "kpp"]
ICES = ["winton", "lebedev", "mylake"]
SLUG = {"k-eps": "keps", "k-omega": "komega", "kpp": "kpp"}
OFFSETS = [-3, -2, -1, 0, 1, 2, 3]


def tag_of(offset):
    return f"p{offset}" if offset > 0 else (f"m{abs(offset)}" if offset < 0 else "0")


def load_temp_and_hice(rundir, name):
    f = nc.Dataset(f"{rundir}/{name}/all.nc")
    temp = np.array(f.variables["temp"][:]).squeeze()  # (time, z)
    hice = np.array(f.variables["Hice"][:]).squeeze() if "Hice" in f.variables else None
    f.close()
    return temp, hice


rows = []
for offset in OFFSETS:
    tag = tag_of(offset)
    if offset == 0:
        rundir = f"{CASE}/runs"  # the real main-event factorial already run
        name_fmt_cell = lambda turb, ice: {
            ("k-eps", "winton"): "baseline", ("k-eps", "lebedev"): "ice_lebedev", ("k-eps", "mylake"): "ice_mylake",
            ("k-omega", "winton"): "turb_komega", ("k-omega", "lebedev"): "cross_komega_lebedev_lagrangian",
            ("k-omega", "mylake"): "cross_komega_mylake_lagrangian",
            ("kpp", "winton"): "turb_kpp", ("kpp", "lebedev"): "cross_kpp_lebedev_lagrangian",
            ("kpp", "mylake"): "cross_kpp_mylake_lagrangian",
        }[(turb, ice)]
        name_fmt_noice = lambda turb: {"k-eps": "noice_keps", "k-omega": "noice_komega", "kpp": "noice_kpp"}[turb]
    else:
        rundir = f"{CASE}/runs_sweep_{tag}"
        name_fmt_cell = lambda turb, ice, tag=tag: f"sw_{tag}_{SLUG[turb]}_{ice}"
        name_fmt_noice = lambda turb, tag=tag: f"sw_{tag}_noice_{SLUG[turb]}"

    n_engaged = 0
    engaged_names = []
    noice_temps = {}
    for turb in CLOSURES:
        noice_temps[turb], _ = load_temp_and_hice(rundir, name_fmt_noice(turb))

    for turb in CLOSURES:
        for ice in ICES:
            temp, hice = load_temp_and_hice(rundir, name_fmt_cell(turb, ice))
            # near-surface only, not the full depth-time array: at the
            # corrected LB9 position (Section 2.1) the full array picks up
            # an unrelated, off-event deep-mixing sensitivity to the
            # presence of a never-triggered ice module (Hice = 0
            # throughout); near-surface is where the sea-ice module can
            # plausibly act directly.
            max_diff = float(np.nanmax(np.abs(temp[:, -1] - noice_temps[turb][:, -1])))
            if max_diff > 1e-6:
                n_engaged += 1
                engaged_names.append(f"{turb}+{ice}")

    # distance-from-freezing proxy: coldest near-surface (shallowest level) temp
    # in the k-epsilon no-ice trajectory at this offset (k-eps never engages ice
    # at any offset tested, so this cleanly isolates the forcing severity)
    coldest_keps = float(np.min(noice_temps["k-eps"][:, -1]))
    coldest_komega = float(np.min(noice_temps["k-omega"][:, -1]))

    rows.append(dict(offset=offset, n_engaged=n_engaged, engaged_cells="|".join(engaged_names),
                      coldest_keps_noice=coldest_keps, coldest_komega_noice=coldest_komega))
    print(f"offset={offset:+d}: {n_engaged}/9 engaged ({', '.join(engaged_names) or 'none'}), "
          f"coldest k-eps(no-ice)={coldest_keps:.3f}C, coldest k-omega(no-ice)={coldest_komega:.3f}C")

with open(f"{TAB_DIR}/T28_LB9_sweep_summary.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["offset", "n_engaged", "engaged_cells", "coldest_keps_noice", "coldest_komega_noice"],
                        lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print("\nT28 written.")
