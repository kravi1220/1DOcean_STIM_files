"""Salt conservation of the GOTM-Icepack exchange (Appendix A).
ocean salt change + ice salt change + dilution by freshwater entering the ocean = 0,
with dilution = rho0 * S_top * swf, swf = (1-aice)*precip + fresh_ai/1000 (as in gotm.F90).
Ice salt = rhoi * vice * Sice (layer mean = bulk, equal layers). Usage: icepack_salt_budget.py <GOTM Icepack run dir>..."""
import os, sys, numpy as np, netCDF4 as nc, datetime as dt
import paths as P
rho0, rhoi = 1027.0, 917.0
met = {}
for l in open(P.GOTM + 'meteo.dat'):
    p = l.split(); met[dt.datetime.strptime(p[0]+' '+p[1], '%Y/%m/%d %H:%M:%S')] = float(p[8])
ms = sorted(met); msec = np.array([(x-ms[0]).total_seconds() for x in ms]); mv = np.array([met[x] for x in ms])
T0 = dt.datetime(2018, 1, 1)
REF = {'C1_bl99': 4.0}          # saltflux_option 'constant' (BL99 run): ice_ref_salinity [ppt]
for run in sys.argv[1:]:
    f = nc.Dataset(run+'/all.nc')
    t = [x.replace(tzinfo=None) for x in nc.num2date(f['time'][:], f['time'].units, only_use_cftime_datetimes=False)]
    S = np.array(f['salt'][:]).squeeze(); h = np.array(f['h'][:]).squeeze(); oc = rho0*np.sum(S*h, axis=1)
    a = np.loadtxt(run+'/icepack_column.dat', comments='#')
    ia = {T0+dt.timedelta(hours=float(x)): i for i, x in enumerate(a[:, 0])}
    cum, res, ice_max = 0.0, [], 0.0
    for j in range(1, len(t)):
        i = ia.get(t[j])
        if i is None: break
        prec = np.interp((t[j]-dt.timedelta(minutes=30)-ms[0]).total_seconds(), msec, mv)
        swf = (1-a[i, 1])*prec + a[i, 9]/1000.0
        cum += rho0*S[j, -1]*swf*3600.0
        ice = rhoi*a[i, 2]*REF.get(os.path.basename(run.rstrip('/')), a[i, 6]); ice_max = max(ice_max, ice)
        res.append(oc[j]-oc[0] + ice + cum)
    res = np.array(res)
    print('%-10s hours %5d  max ice salt %7.1f g/m2  max |residual| %7.2f g/m2 (%.2f %% of max ice salt)  final %7.2f'
          % (run, len(res), ice_max, np.abs(res).max(), 100*np.abs(res).max()/max(ice_max, 1e-9), res[-1]))
