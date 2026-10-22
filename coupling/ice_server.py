"""Icepack column served to the Oceananigans LES over stdin/stdout (LES x Icepack coupling).

Wraps the same persistent Icepack process and Python driver used for MPAS-ADC x Icepack
(icepack_pipe.py): 32 sub-steps of 112.5 s per hour,
GOTM-identical atmosphere (ERA5 meteo.dat) and ERA5 snowfall, mushy-layer thermodynamics.

Usage (started by les_kg6_era5_icepack.jl):  python3 ice_server.py <run_dir> <start ISO>
  stdin : one line per hour "ISO sst_insitu_C sss_psu h_top_m" for the hour starting at ISO
  stdout: "READY" once, then per hour
          "aice vice vsno Tsfc Sice fhocn fswthru fresh fsalt frazil_heat aice_mean"
          (fluxes: hour means, grid-box means, positive into the ocean; MPAS-driver convention)
  "END" closes Icepack.
"""
import os, sys, datetime as dt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gotm_airsea import Meteo
from icepack_pipe import IcepackPipe, snow_series, KEYS

CASE = os.environ.get('ICE_CASE', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments', 'gotm')).rstrip('/') + '/'   # meteo.dat, meteo_snow.dat
LAT, LON = float(os.environ.get('LES_LAT', 67.58005324675324)), float(os.environ.get('LES_LON', -23.92983766233766))
run_dir, start = sys.argv[1], dt.datetime.fromisoformat(sys.argv[2])
met = Meteo(CASE + 'meteo.dat'); snow = snow_series(CASE + 'meteo_snow.dat')
pipe = IcepackPipe(run_dir, LAT, LON, 112.5, start, 32)
aice_prev = 0.0
print('READY', flush=True)
for line in sys.stdin:
    p = line.split()
    if not p or p[0] == 'END':
        break
    t = dt.datetime.fromisoformat(p[0]); sst, sss, htop = map(float, p[1:4])
    ice = pipe.hour(t, met, snow, sst, sss, htop / pipe.nsub)   # top-layer potential spread over sub-steps
    aice_mean = 0.5 * (aice_prev + ice['aice']); aice_prev = ice['aice']
    print(' '.join('%.10e' % ice[k] for k in KEYS) + ' %.10e' % aice_mean, flush=True)
pipe.close()
