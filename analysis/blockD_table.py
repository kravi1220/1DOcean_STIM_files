"""Block D 2x2 interaction table at uniform 5 m: {GOTM k-eps, MPAS-ADC (salt-fixed)} x {no ice, Icepack mushy}.
Interaction = (ADC_ice - ADC_noice) - (keps_ice - keps_noice) (Section 4.5, Fig. 9).
Usage: python3 blockD_table.py > ../results/tables/blockD_table.txt"""
import datetime as dt, numpy as np
import kg6_common as A, paths as P
MC, K = P.MPAS, P.GOTM
cells = {('keps', 'noice'): ('g', K + 'blockA/gotm/keps_dz5_dt112.5'), ('keps', 'ice'): ('g', K + 'blockD/keps_dz5_icepack'),
         ('ADC', 'noice'): ('m', MC + 'blockA/KG6_dz5'), ('ADC', 'ice'): ('m', MC + 'blockD/KG6_dz5_icepack')}
C = A.ctd(); T0 = dt.datetime(2018, 1, 1); res = {}
for key, (kind, path) in cells.items():
    M = A.gotm(path) if kind == 'g' else A.mpas(path)
    if kind == 'm': M = {k: (v[:-1] if isinstance(v, np.ndarray) and v.shape[:1] == M['t'].shape else v) for k, v in M.items()}
    t_end = dt.datetime(2019, 8, 10, 11); M = A.window(M, t_end); mld, hc = A.derived(M); Ts = M['T'][:, 0]
    rt, rs, _ = A.rmse(M, C, t_end)
    r = dict(minTs=Ts.min(), mld=mld.max(), hc=hc[-1] / 1e9, rt=rt, rs=rs, Ssurf_min=M['S'][:, 0].min())
    if key[1] == 'ice':
        a = np.loadtxt(path + '/icepack_column.dat', comments='#'); on = np.where(a[:, 1] > 0.01)[0]
        tt = [T0 + dt.timedelta(hours=float(h)) for h in a[:, 0]]
        r.update(ice_on=tt[on[0]], ice_off=tt[on[-1]], days=len(on) / 24, hmax=a[:, 4].max(), vmax=a[:, 2].max(), amax=a[:, 1].max())
    res[key] = r
f = lambda d: d.strftime('%m-%d')
print('# Block D, uniform 5 m, full year 2018-08-06 .. 2019-08-10 (4 casts)')
print('%-12s %8s %7s %9s %7s %7s %8s | %6s %6s %6s %6s %5s' % ('cell', 'minTs', 'MLDmax', 'dHC300', 'RMSE_T', 'RMSE_S', 'minSsfc', 'on', 'off', 'days', 'Hmax', 'amax'))
for k, r in res.items():
    s = '%-12s %8.3f %7.1f %9.3f %7.4f %7.4f %8.3f' % ('%s %s' % k, r['minTs'], r['mld'], r['hc'], r['rt'], r['rs'], r['Ssurf_min'])
    if 'ice_on' in r: s += ' | %6s %6s %6.1f %6.3f %5.2f' % (f(r['ice_on']), f(r['ice_off']), r['days'], r['hmax'], r['amax'])
    print(s)
print('\n# interaction = (ADC_ice - ADC_noice) - (keps_ice - keps_noice)')
for v in ['minTs', 'mld', 'hc', 'rt', 'rs', 'Ssurf_min']:
    ia = res['ADC', 'ice'][v] - res['ADC', 'noice'][v]; ik = res['keps', 'ice'][v] - res['keps', 'noice'][v]
    print('%-10s ice effect: keps %+8.3f  ADC %+8.3f  -> interaction %+8.3f' % (v, ik, ia, ia - ik))
print('%-10s closure effect on ice (ADC - keps): onset %+d d, end %+d d, duration %+.1f d, Hmax %+.3f m'
      % ('ice', (res['ADC', 'ice']['ice_on'] - res['keps', 'ice']['ice_on']).days, (res['ADC', 'ice']['ice_off'] - res['keps', 'ice']['ice_off']).days,
         res['ADC', 'ice']['days'] - res['keps', 'ice']['days'], res['ADC', 'ice']['hmax'] - res['keps', 'ice']['hmax']))
