"""Like-for-like comparison from the observed 2018-10-15 cast: ERA5-driven LES vs GOTM k-eps/KPP
and MPAS-ADC (all with ERA5 Kondo fluxes from their own SST), scored at the 2019-02-15 21:35 cast
and against the LES (Section 4.3.2). Usage: python3 compare_oct15.py > ../results/tables/compare_oct15.txt"""
import os, sys, numpy as np, netCDF4 as nc, datetime as dt, gsw
import kg6_common as A, paths as P
LAT, LON = A.LAT, A.LON
LES, MC = P.LES, P.MPAS
def les(tag):
    f = LES + 'les_kg6_era5_%s_profiles.nc' % tag
    if not os.path.exists(f): return None
    d = nc.Dataset(f)
    t0 = dt.datetime(2018, 10, 15, 18); t = np.array([t0 + dt.timedelta(seconds=float(x)) for x in d['time'][:]])
    z = -np.array(d['z_aac'][:]); o = np.argsort(z)
    return dict(t=t, z=z[o], T=np.array(d['T'][:])[:, o], S=gsw.SP_from_SA(np.array(d['S'][:])[:, o], gsw.p_from_z(-z[o], LAT), LON, LAT))
def adc(dz):
    p = MC + 'oct15/KG6oct_dz%d' % dz
    if not os.path.exists(p + '/hourly/KG6oct_adc.nc'): return None
    M = A.mpas(p); return {k: (v[:-1] if isinstance(v, np.ndarray) and v.shape[:1] == M['t'].shape else v) for k, v in M.items()}
R = {}
for g in ['2m']:
    m = les('oct15_' + g)
    if m is not None: R['LES ' + g.split('_')[0]] = m
for c in ['keps', 'kpp']:
    for dz in [2, 5, 10]: R['GOTM %s %dm' % (c, dz)] = A.gotm(P.GOTM + 'oct15/runs/%s_dz%d' % (c, dz))
for dz in [2, 5, 10]:          # salt-fixed ADC build; the 2 m run stops in its first hour (Section 4.2.3)
    m = adc(dz)
    if m is not None and len(m['t']) > 24: R['ADC %dm' % dz] = m
C = A.ctd(); ct, cd, cT, cS = C
ic = [k for k, x in enumerate(ct) if x.date() == dt.date(2019, 2, 15)][0]
ok = ~np.isnan(cT[ic]) & ~np.isnan(cS[ic]); dv = cd[ok]; q = np.argsort(dv); oT, oS = cT[ic][ok][q], cS[ic][ok][q]; dv = dv[q]
tc = ct[ic]; ref = R.get('LES 2m')
print('# start 2018-10-15 18:00 (observed cast); scored at the cast %s; profile RMS over 2-150 m' % tc)
print('%-14s %10s %8s %8s %8s %9s %9s %9s  %s' % ('run', 'end', 'T0-2m', 'S0-2m', 'MLD', 'RMS_T obs', 'RMS_S obs', 'RMS_T LES2', 'min Tsurf (date)'))
print('%-14s %10s %8.3f %8.3f' % ('CTD obs', tc.strftime('%m-%d %H:%M'), oT[0], oS[0]))
zc = np.arange(2.0, 150.0, 1.0)
for n, M in R.items():
    if M['t'][-1] < tc - dt.timedelta(hours=1): print('%-14s ends %s (before the cast)' % (n, M['t'][-1])); continue
    k = int(np.argmin([abs((x - tc).total_seconds()) for x in M['t']]))
    Tm = np.interp(zc, M['z'], M['T'][k]); Sm = np.interp(zc, M['z'], M['S'][k])
    rt = np.sqrt(np.mean((Tm - np.interp(zc, dv, oT))**2)); rs = np.sqrt(np.mean((Sm - np.interp(zc, dv, oS))**2))
    rl = np.sqrt(np.mean((Tm - np.interp(zc, ref['z'], ref['T'][int(np.argmin([abs((x - tc).total_seconds()) for x in ref['t']]))]))**2)) if ref else np.nan
    mm = A.derived(dict(t=M['t'][k:k+1], z=M['z'], T=M['T'][k:k+1], S=M['S'][k:k+1]))[0][0]
    Ts = M['T'][:, 0]; j = int(np.argmin(Ts))
    print('%-14s %10s %8.3f %8.3f %8.1f %9.4f %9.4f %9.4f  %.3f (%s)' % (n, M['t'][k].strftime('%m-%d %H:%M'), M['T'][k][0], M['S'][k][0], mm, rt, rs, rl, Ts[j], M['t'][j].strftime('%m-%d')))
