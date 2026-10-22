"""Block R metrics (KG6; Section 4.1, convergence and initial-condition ensemble).
Usage: python3 blockR_metrics.py <GOTM run dir> [...]   (see scripts/run_analysis.sh)"""
import sys, numpy as np, netCDF4 as nc, datetime as dt, gsw
import paths as P
CTD = P.DATA + 'ctd_profiles_KG6.nc'   # not distributed (README, 'CTD observations')
VAL = ['2018-10-15', '2019-02-15', '2019-05-18', '2019-08-10']
TF = -1.73   # KG6 local freezing point (TEOS-10, manuscript value)
LAT, LON = 67.58005324675324, -23.92983766233766
C = nc.Dataset(CTD); ct = [x.replace(tzinfo=None) for x in nc.num2date(C['time'][:], C['time'].units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)]
cd, cT, cS = np.array(C['depth'][:]), np.array(C['temperature'][:]), np.array(C['salinity'][:])
# CTD holds in-situ T / practical S; GOTM (TEOS-10) 'temp' is CT and 'salt' is SA.
# Compare as CT and SP (fixed 2026-10-07, as blockA_analysis.py)
_p = gsw.p_from_z(-cd, LAT)[None, :]; cT = gsw.CT_from_t(gsw.SA_from_SP(cS, _p, LON, LAT), cT, _p)
def metrics(d):
    f = nc.Dataset(d + '/all.nc')
    t = [x.replace(tzinfo=None) for x in nc.num2date(f['time'][:], f['time'].units, only_use_cftime_datetimes=False)]
    z = -np.array(f['z'][0]).squeeze(); o = np.argsort(z); z = z[o]
    T = np.array(f['temp'][:]).squeeze()[:, o]; S = gsw.SP_from_SA(np.array(f['salt'][:]).squeeze()[:, o], gsw.p_from_z(-z, LAT), LON, LAT)
    H = np.array(f['Hice'][:]).squeeze() if 'Hice' in f.variables else np.zeros(len(t))
    Ts = T[:, 0]; ice = np.where(H > 0.01)[0]
    rt, rs = [], []
    for v in VAL:
        vd = dt.datetime.strptime(v, '%Y-%m-%d'); i = min(range(len(ct)), key=lambda k: abs((ct[k] - vd).total_seconds()))
        ok = ~np.isnan(cT[i]) & ~np.isnan(cS[i]); dv = cd[ok]; q = np.argsort(dv)
        k = int(np.argmin([abs((x - ct[i]).total_seconds()) for x in t]))
        zc = np.arange(2.0, min(dv.max(), z.max()), 2.0)
        rt.append(np.sqrt(np.mean((np.interp(zc, z, T[k]) - np.interp(zc, dv[q], cT[i][ok][q]))**2)))
        rs.append(np.sqrt(np.mean((np.interp(zc, z, S[k]) - np.interp(zc, dv[q], cS[i][ok][q]))**2)))
    frz = next((t[k] for k in range(len(t)) if Ts[k] <= TF), None)
    return dict(minTs=Ts.min(), date_minTs=t[int(Ts.argmin())], first_freeze=frz,
                Hmax=H.max(), ice_days=(H > 0.01).sum()*(t[1]-t[0]).total_seconds()/86400,
                first_ice=t[ice[0]] if ice.size else None, last_ice=t[ice[-1]] if ice.size else None,
                rmseT=np.mean(rt), rmseS=np.mean(rs), rmseT_by=rt, Tsurf=Ts, t=t)
if __name__ == '__main__':
    for d in sys.argv[1:]:
        m = metrics(d); fmt = lambda x: x.strftime('%Y-%m-%d') if x else '-'
        print('%-28s minTs %7.3f (%s) Tf-reached %s  Hmax %5.3f  ice %5.1f d  %s..%s  RMSE_T %.4f  RMSE_S %.4f' % (
            d, m['minTs'], fmt(m['date_minTs']), fmt(m['first_freeze']), m['Hmax'], m['ice_days'], fmt(m['first_ice']), fmt(m['last_ice']), m['rmseT'], m['rmseS']))
