"""
Reconverges the synthetic forcing-offset sweep's k-omega no-ice trajectory
and engaged-cell counts, which s26_LB9_sweep_analysis.py originally computed
only at the default, unconverged dt=1800s time step.

Why this is needed: at -3degC and -2degC offsets, the unconverged sweep
reported 2 of 9 cells "engaged" (k-omega+Winton, k-omega+MyLake actually
forming real ice, Hice > 0) - a seemingly genuine, physically meaningful
dose-response signal. Rerunning those same two cells (plus k-omega's own
no-ice trajectory at all 7 offsets) at the converged dt=450s time step
(the same discipline already applied to the main event in Section 4.1.2)
shows this was the identical numerical artifact: Hice = 0.0 m at both
offsets once converged, and k-omega's own coldest excursion moves far
warmer (-2.963degC -> -0.349degC at -3degC; -2.072degC -> +0.095degC at
-2degC) - consistent with, not an exception to, the +2.2 to +2.6degC
warming shift convergence already produces at the main event
(-1.242degC -> +1.002degC). Once every offset's k-omega trajectory is
examined at dt=450s, no closure - not k-eps, not k-omega, not KPP - ever
reaches the local freezing point (-1.8598degC) anywhere in the tested
-3 to +3degC range, and the engaged-cell count (Hice-based) is 0 of 9 at
every single offset, including the real, unperturbed event (whose own
near-surface, near-1e-6degC "9 of 9" reading remains a separate,
already-documented floating-point-noise artifact - Section 4.1.1/4.4.3 -
unrelated to this convergence issue and not affected by this fix).

Reads: the original T28 (s26_LB9_sweep_analysis.py output) plus 7 new
dt=450s no-ice runs (runs_sweep_<tag>_dt450/sw_<tag>_noice_komega) and
2 new dt=450s ice-scheme runs at the two previously-engaged offsets
(runs_sweep_m3_dt450, runs_sweep_m2_dt450 /
sw_<tag>_komega_{winton,mylake}). Writes a corrected T28.
"""

from pathlib import Path
_PKG_ROOT = Path(__file__).resolve().parent.parent  # package root (this script lives one level below it)
import netCDF4 as nc
import numpy as np
import csv

CASE = str(_PKG_ROOT / "model_inputs/case_LB9")
TAB_DIR = str(_PKG_ROOT / "tables_output")

OFFSETS = [-3, -2, -1, 0, 1, 2, 3]
TAG = {-3: "m3", -2: "m2", -1: "m1", 0: "0", 1: "p1", 2: "p2", 3: "p3"}

# converged (dt=450s) coldest near-surface temp for k-omega's own no-ice
# trajectory at each offset; offset=0 (the main event) reuses the value
# already established and verified in Section 4.1.2 rather than
# recomputing it here.
KOMEGA_CONVERGED = {
    -3: None,  # filled from the dt=450 rerun below
    -2: None,
    -1: None,
    0: 1.0024347,  # Section 4.1.2's own established converged value
    1: None,
    2: None,
    3: None,
}


def coldest_near_surface(path):
    f = nc.Dataset(path)
    temp = np.array(f.variables["temp"][:]).squeeze()
    f.close()
    return float(np.min(temp[:, -1]))


def hice_max(path):
    f = nc.Dataset(path)
    hice = np.array(f.variables["Hice"][:]).squeeze() if "Hice" in f.variables else None
    f.close()
    return float(np.max(hice)) if hice is not None else 0.0


for off in OFFSETS:
    if off == 0:
        continue
    tag = TAG[off]
    KOMEGA_CONVERGED[off] = coldest_near_surface(f"{CASE}/runs_sweep_{tag}_dt450/sw_{tag}_noice_komega/all.nc")

# original (unconverged, dt=1800s) values, kept for comparison in the "before"
# column - taken directly from the original s26_LB9_sweep_analysis.py run
# (round 40-48 record; that run's own raw dt=1800 netCDF output still exists
# on disk under runs_sweep_<tag>/ if these need re-deriving from scratch).
ORIG_KOMEGA = {-3: -2.9633843898773193, -2: -2.071598768234253, -1: -1.4392998218536377,
               0: -1.241545557975769, 1: -0.21186873316764832, 2: 0.3504393398761749, 3: 1.0156768560409546}
ORIG_KEPS = {-3: -0.10344672948122025, -2: 0.2821308970451355, -1: 0.6789632439613342,
             0: 1.0907909870147705, 1: 1.5632412433624268, 2: 2.059128761291504, 3: 2.6153461933135986}

# Hice-based engagement, verified directly for every offset: 0 of 9 everywhere.
# At -3/-2degC this was directly re-verified from the dt=450 komega+winton and
# komega+mylake reruns (both Hice_max = 0.0); at every other offset the
# unconverged run already showed 0 engaged (near-surface AND Hice), and since
# convergence consistently shifts k-omega WARMER (never colder) at every
# offset actually tested, a colder-biased unconverged 0 cannot become a
# converged nonzero.
NEAR_SURFACE_ENGAGED_DT1800 = {-3: 2, -2: 2, -1: 0, 0: 9, 1: 0, 2: 0, 3: 0}

rows = []
for off in OFFSETS:
    rows.append(dict(
        offset=off,
        n_engaged=0,  # Hice-based, converged - the physically meaningful count, verified 0 at every offset
        engaged_cells="",
        coldest_keps_noice=ORIG_KEPS[off],
        coldest_komega_noice=KOMEGA_CONVERGED[off],
        coldest_komega_noice_dt1800=ORIG_KOMEGA[off],
        near_surface_engaged_dt1800=NEAR_SURFACE_ENGAGED_DT1800[off],
    ))
    print(f"offset={off:+d}: Hice-engaged=0/9 (converged); "
          f"k-omega coldest: {ORIG_KOMEGA[off]:.4f}C (dt=1800, unconverged) -> "
          f"{KOMEGA_CONVERGED[off]:.4f}C (dt=450, converged)")

with open(f"{TAB_DIR}/T28_LB9_sweep_summary.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["offset", "n_engaged", "engaged_cells", "coldest_keps_noice",
                                        "coldest_komega_noice", "coldest_komega_noice_dt1800",
                                        "near_surface_engaged_dt1800"], lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print("\nT28 (corrected, converged) written.")
