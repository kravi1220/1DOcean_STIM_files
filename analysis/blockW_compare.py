"""Block W: 1D closures vs LES over the fixed-flux windows B1 (early-winter convection,
2018-11-28 +5 d) and B2 (cooling to Tf / ice growth, 2019-02-08 +5 d).
All runs: same initial CT, SA, u, v (les_kg6_<B>.nc), same prescribed hourly fluxes, 160 m column.
Reference: LES 2 m (WENO-9 implicit LES; 2 m converged, BLOCK_B_STATUS.md).
Metrics at 24, 72, 120 h: T0-10 = mean CT over 0-10 m, S0-10, MLD (sigma0 +0.03 from 10 m),
RMS_T / RMS_S vs LES over 0-150 m; heat-budget residual = dHC(0-160 m) - integral(Q + I0).
Sections 4.3.1 (Fig. 7). Usage: python3 blockW_compare.py > ../results/tables/blockW_compare.txt"""
import os, datetime as dt
import numpy as np, netCDF4 as nc, gsw
import kg6_common as A, paths as P
LES, MC = P.LES, P.MPAS
RHO0, CP = 1027.0, 3991.86795711963
zf = np.arange(0.25, 150.0, 0.5)                      # common grid for profiles (depth, m)


def les(w, dx):
    d = nc.Dataset(LES + 'les_kg6_%s_dx%.1f_profiles.nc' % (w, dx))
    t0 = dt.datetime.fromisoformat(nc.Dataset(LES + 'les_kg6_%s.nc' % w).window[:19])
    z = -np.array(d['z_aac'][:]); o = np.argsort(z)
    SA = np.array(d['S'][:])[:, o]
    return dict(t=np.array([t0 + dt.timedelta(seconds=float(s)) for s in d['time'][:]]), z=z[o],
                T=np.array(d['T'][:])[:, o], S=gsw.SP_from_SA(SA, gsw.p_from_z(-z[o], A.LAT), A.LON, A.LAT), dz=dx)


def gotm(path):
    M = A.gotm(path); M['dz'] = None; return M


def adc(path):
    M = A.mpas(path); M['dz'] = None; return M


def at(M, t):
    k = int(np.argmin([abs((x - t).total_seconds()) for x in M['t']]))
    return k, abs((M['t'][k] - t).total_seconds())


def prof(M, k, v):
    return np.interp(zf, M['z'], M[v][k])


def heat_budget(M, t0, t1, F):
    k0, _ = at(M, t0); k1, _ = at(M, t1)
    dzs = np.gradient(M['z'])
    dhc = RHO0*CP*np.sum((M['T'][k1] - M['T'][k0])*dzs)
    s0 = (M['t'][k0] - t0).total_seconds(); s1 = (M['t'][k1] - t0).total_seconds()
    ss = np.linspace(s0, s1, 2001)
    flux = np.trapezoid(np.interp(ss, F['time'], F['Q'] + F['I0']), ss)
    return (dhc - flux)/1e6, dhc/1e6                                   # MJ/m2


for w in ['B1', 'B2']:
    F = nc.Dataset(LES + 'les_kg6_%s.nc' % w); F = {k: np.array(F[k][:]) for k in ('time', 'Q', 'I0')}
    R = {'LES 2m': les(w, 2.0), 'LES 4m': les(w, 4.0)}
    if w == 'B1':
        if os.path.exists(LES + 'les_kg6_%s_dx1.0_profiles.nc' % w):   # 1 m run (24 h, B1 only)
            R['LES 1m (24h)'] = les(w, 1.0)
    for c in ['keps', 'kpp']:
        for dz in [2, 5, 10]:
            p = P.GOTM + 'blockW/gotm/%s_%s_dz%d' % (w, c, dz)
            if os.path.exists(p + '/all.nc'):
                R['GOTM %s %dm' % (c, dz)] = gotm(p)
    for dz in [2, 5, 10]:
        p = MC + 'blockW/%s_dz%d' % (w, dz)
        if os.path.exists(p + '/hourly/%s_adc.nc' % w):
            M = adc(p)
            if len(M['t']) > 2:
                R['ADC %dm' % dz] = M
    ref = R['LES 2m']; t0 = ref['t'][0]
    print('\n# %s: start %s; reference LES 2 m; T0-10 [C], MLD [m], RMS over 0-150 m vs LES 2 m' % (w, t0))
    print('%-14s %5s | %8s %8s %6s %7s %7s | %s' % ('run', 'hour', 'T0-10', 'S0-10', 'MLD', 'RMS_T', 'RMS_S', 'heat residual MJ/m2 (dHC)'))
    for n, M in R.items():
        mld, _ = A.derived(M)
        for h in (24, 72, 120):
            t = t0 + dt.timedelta(hours=h)
            if M['t'][-1] < t - dt.timedelta(minutes=30):
                continue
            k, off = at(M, t); kr, _ = at(ref, t)
            T, S = prof(M, k, 'T'), prof(M, k, 'S'); Tr, Sr = prof(ref, kr, 'T'), prof(ref, kr, 'S')
            top = zf <= 10.0
            line = '%-14s %5d | %8.3f %8.3f %6.1f %7.4f %7.4f |' % (n, h, T[top].mean(), S[top].mean(), mld[k],
                                                               np.sqrt(np.mean((T - Tr)**2)), np.sqrt(np.mean((S - Sr)**2)))
            if h == 120 or (n.startswith('LES 1m') and h == 24):
                res, dhc = heat_budget(M, t0, t, F)
                line += ' %+7.2f (%+.1f)' % (res, dhc)
            print(line)
