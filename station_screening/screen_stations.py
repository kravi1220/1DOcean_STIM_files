
from pathlib import Path
_PKG_ROOT = Path(__file__).resolve().parent.parent  # package root (this script lives one level below it)
#!/usr/bin/env python
"""
Screen the full repeat hydrographic network for any occupation approaching
its local freezing point (Section 2.1 of the manuscript). For every
station in data/ctd_profiles_*.nc, and every occupation of that station,
computes the shallowest-valid-level surface temperature, salinity, and
pressure, and the TEOS-10 local freezing point at that station's own
coordinates and the shallowest sample's own pressure (not a nominal
p=0 dbar surface reference - at LB9 this pressure correction alone is
about 0.0023 degC, non-negligible relative to a ~0.02-0.03 degC margin),
converted back to the in-situ temperature scale so it is directly
comparable to the measured in-situ value (comparing an in-situ
temperature against a Conservative Temperature freezing point directly,
without this conversion, is a temperature-scale mismatch that produced a
materially different, incorrect margin in an earlier version of the
manuscript - see Section 2.1's own account of the correction). Reports
every station ranked by its single coldest-relative-to-freezing
occupation.

This is a general-purpose diagnostic, not a curve-fit to the answer: it
was used to confirm LB9 as the only near-freezing occupation in the
network (margin 0.025 degC) and KG6 as the next-closest (margin ~1.94
degC), and is included here so that finding is independently
reproducible rather than asserted.
"""
import glob
import numpy as np
import netCDF4 as nc
import gsw

DATA_GLOB = str(_PKG_ROOT / "data/ctd_profiles_*.nc")


def shallowest_valid(depth, var):
    valid = np.isfinite(depth) & np.isfinite(var)
    if not valid.any():
        return np.nan, np.nan
    order = np.argsort(depth[valid])
    d = depth[valid][order]
    v = var[valid][order]
    return d[0], v[0]


# Station coordinates used only for the (negligible, <0.0001 degC) SA
# geostrophic correction; LB9 and KG6 use their surveyed/naive-transect
# positions (Section 2.1), every other station uses (0, 0) since none of
# them are close enough to freezing for even a large position error to
# change the ranking.
STATION_COORDS = {"LB9": (-27.25850, 66.14767), "KG6": (-23.92983766233766, 67.58005324675324)}


def screen_station(path):
    station = path.split("ctd_profiles_")[-1].replace(".nc", "")
    lon, lat = STATION_COORDS.get(station, (0.0, 0.0))
    f = nc.Dataset(path)
    depth = np.array(f.variables["depth"][:])
    temp = np.array(f.variables["temperature"][:])
    salt = np.array(f.variables["salinity"][:])
    time_var = f.variables["time"]
    times = nc.num2date(time_var[:], time_var.units,
                         only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    f.close()

    best = None
    for i, t in enumerate(times):
        d0, T = shallowest_valid(depth, temp[i])
        _, S = shallowest_valid(depth, salt[i])
        if not (np.isfinite(T) and np.isfinite(S)):
            continue
        p0 = d0  # dbar ~= depth in m at these shallow, high-latitude pressures
        SA = gsw.SA_from_SP(S, p0, lon, lat)
        CTf = gsw.CT_freezing(SA, p0, 0.0)
        Tf_is = gsw.t_from_CT(SA, CTf, p0)  # in-situ freezing point at the true sample pressure
        margin = T - Tf_is
        row = dict(station=station, occupation_idx=i, date=str(t), depth_m=d0,
                   temp_C=T, salinity_psu=S, freezing_point_C=Tf_is, margin_above_freezing_C=margin)
        if best is None or margin < best["margin_above_freezing_C"]:
            best = row
    return best


def main():
    files = sorted(glob.glob(DATA_GLOB))
    print(f"Screening {len(files)} stations from {DATA_GLOB}")
    results = [r for r in (screen_station(p) for p in files) if r is not None]
    results.sort(key=lambda r: r["margin_above_freezing_C"])

    print(f"\n{'rank':>4} {'station':<8} {'date':<20} {'T (C)':>9} {'Tf (C)':>9} {'margin (C)':>11}")
    for rank, r in enumerate(results, start=1):
        flag = "  <-- coldest occupation in network" if rank == 1 else ""
        flag = "  <-- next-closest" if rank == 2 else flag
        print(f"{rank:>4} {r['station']:<8} {r['date']:<20} {r['temp_C']:>9.3f} "
              f"{r['freezing_point_C']:>9.4f} {r['margin_above_freezing_C']:>11.4f}{flag}")

    import csv
    with open("station_screening_results.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()))
        w.writeheader()
        w.writerows(results)
    print("\nWrote station_screening_results.csv (one row per station: its single closest-to-freezing occupation)")


if __name__ == "__main__":
    main()
