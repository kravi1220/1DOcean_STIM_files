"""Common diagnostics for KG6: loaders for GOTM, MPAS-ADC and the CTD casts, and the derived quantities.

Every model is reduced to the same quantities from its own T/S profiles:
- CT (GOTM 'temp' is CT; MPAS potential T -> CT with gsw; CTD in-situ T -> CT), practical S
  (GOTM 'salt' is absolute salinity -> SP; MPAS and CTD are practical);
- near-surface T = shallowest layer; date it first reaches the local freezing point
  (-1.73 C, the manuscript's KG6 TEOS-10 value) and its minimum;
- MLD: sigma0 exceeds its value at 10 m by 0.03 kg/m3 (profiles linearly interpolated);
- upper-300 m heat content anomaly relative to the initial profile
  (rho0 = 1027 kg/m3, cp = 3991.87 J/kg/K, TEOS-10 cp0);
- RMSE_T / RMSE_S at the 4 validation casts (2 m grid, s16 method).
"""
import os, glob, datetime as dt
import numpy as np, netCDF4 as nc, gsw
import paths as P

LAT, LON, TF = 67.58005324675324, -23.92983766233766, -1.73
RHO0, CP = 1027.0, 3991.86795711963
VAL = ['2018-10-15', '2019-02-15', '2019-05-18', '2019-08-10']


def gotm(run):
    f = nc.Dataset(run + '/all.nc')
    t = np.array([x.replace(tzinfo=None) for x in nc.num2date(f['time'][:], f['time'].units, only_use_cftime_datetimes=False)])
    z = -np.array(f['z'][0]).squeeze(); o = np.argsort(z)
    # GOTM (TEOS-10) 'salt' is ABSOLUTE salinity: convert to practical, as all S here
    # (fixed 2026-10-07; before, GOTM S was ~0.15 too high vs CTD and MPAS)
    SA = np.array(f['salt'][:]).squeeze()[:, o]
    return dict(t=t, z=z[o], T=np.array(f['temp'][:]).squeeze()[:, o],
                S=gsw.SP_from_SA(SA, gsw.p_from_z(-z[o], LAT), LON, LAT))


def mpas(case):
    f = nc.Dataset(glob.glob(case + '/hourly/*_adc.nc')[0])      # <site>_adc.nc written by coupling/mpas_driver.py
    t = np.array([dt.datetime(2018, 1, 1) + dt.timedelta(hours=float(h)) for h in f['time'][:]])
    z = -np.array(f['zmid'][:]); PT = np.array(f['temperature'][:]); SP = np.array(f['salinity'][:])
    p = gsw.p_from_z(-z, LAT); SA = gsw.SA_from_SP(SP, p, LON, LAT)
    return dict(t=t, z=z, T=gsw.CT_from_pt(SA, PT), S=SP)


def derived(M):
    z = M['z']; p = gsw.p_from_z(-z, LAT)
    SA = gsw.SA_from_SP(M['S'], p, LON, LAT); sig = gsw.sigma0(SA, M['T'])
    zf = np.arange(0.0, 500.0, 0.5)
    mld = np.empty(len(M['t']))
    for k in range(len(M['t'])):
        s = np.interp(zf, z, sig[k]); s10 = np.interp(10.0, z, sig[k])
        j = np.where((zf > 10.0) & (s > s10 + 0.03))[0]
        mld[k] = zf[j[0]] if j.size else zf[-1]
    w = (z <= 300.0)
    dzs = np.gradient(z)
    hc = RHO0*CP*np.sum((M['T'][:, w] - M['T'][0, w])*dzs[w], axis=1)
    return mld, hc


def ctd():
    f = P.DATA + 'ctd_profiles_KG6.nc'
    if not os.path.exists(f):
        raise SystemExit('CTD casts not found: %s (not distributed; request them from MFRI, see README)' % f)
    C = nc.Dataset(f)
    t = [x.replace(tzinfo=None) for x in nc.num2date(C['time'][:], C['time'].units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)]
    # CTD: in-situ temperature and practical salinity (they are the GOTM t/s_prof inputs, types
    # in-situ/practical); models are compared as CT, so convert (fixed 2026-10-07)
    d = np.array(C['depth'][:]); tin = np.array(C['temperature'][:]); SP = np.array(C['salinity'][:])
    p = gsw.p_from_z(-d, LAT)[None, :]
    CT = gsw.CT_from_t(gsw.SA_from_SP(SP, p, LON, LAT), tin, p)
    return t, d, CT, SP


def rmse(M, C, t_end):
    ct, cd, cT, cS = C; rt, rs, used = [], [], []
    for v in VAL:
        vd = dt.datetime.strptime(v, '%Y-%m-%d')
        i = min(range(len(ct)), key=lambda k: abs((ct[k] - vd).total_seconds()))
        if ct[i] > t_end + dt.timedelta(days=1):   # final cast is taken later on the last model day
            continue
        ok = ~np.isnan(cT[i]) & ~np.isnan(cS[i]); dv = cd[ok]; q = np.argsort(dv)
        k = int(np.argmin([abs((x - ct[i]).total_seconds()) for x in M['t']]))
        zc = np.arange(2.0, min(dv.max(), M['z'].max()), 2.0)
        rt.append(np.sqrt(np.mean((np.interp(zc, M['z'], M['T'][k]) - np.interp(zc, dv[q], cT[i][ok][q]))**2)))
        rs.append(np.sqrt(np.mean((np.interp(zc, M['z'], M['S'][k]) - np.interp(zc, dv[q], cS[i][ok][q]))**2)))
        used.append(v)
    return (np.mean(rt) if rt else np.nan), (np.mean(rs) if rs else np.nan), used


def window(M, t_end):
    m = M['t'] <= t_end
    return dict(t=M['t'][m], z=M['z'], T=M['T'][m], S=M['S'][m])
