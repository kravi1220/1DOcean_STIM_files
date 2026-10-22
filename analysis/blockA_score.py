"""Block A scoring (Fig. 3, Section 4.2): GOTM k-eps / KPP and MPAS-ADC (outer step 9.375 s) on uniform 2, 5, 10 m grids. ADC 2 m stops 2019-06-23
(ADC 'wt out of range'), so two windows: full year (ADC 5, 10 m) and common to 2019-06-22.
Usage: python3 blockA_score.py > ../results/tables/blockA_score.txt"""
import numpy as np, datetime as dt
import kg6_common as A, paths as P
C = A.ctd(); R = {}
for c in ['keps', 'kpp']:
    for dz in [2, 5, 10]: R['GOTM %s %dm' % (c, dz)] = A.gotm(P.GOTM + 'blockA/gotm/%s_dz%d_dt112.5' % (c, dz))
for dz in [2, 5, 10]:
    M = A.mpas(P.MPAS + 'blockA/KG6_dz%d' % dz)
    R['ADC %dm' % dz] = {k: (v[:-1] if isinstance(v, np.ndarray) and v.shape[:1] == M['t'].shape else v) for k, v in M.items()}
def table(names, t_end, title):
    print('\n# %s (window 2018-08-06 12:00 -> %s)' % (title, t_end))
    print('%-15s %8s %11s %7s %9s %8s %8s  casts' % ('run', 'minTs', 'Tf reached', 'MLDmax', 'dHC300end', 'RMSE_T', 'RMSE_S'))
    out = {}
    for n in names:
        M = A.window(R[n], t_end); mld, hc = A.derived(M); Ts = M['T'][:, 0]
        frz = next((M['t'][k] for k in range(len(Ts)) if Ts[k] <= A.TF), None); rt, rs, used = A.rmse(M, C, t_end)
        out[n] = dict(minTs=Ts.min(), frz=frz, mld=mld.max(), rt=rt, rs=rs)
        print('%-15s %8.3f %11s %7.1f %9.3f %8.4f %8.4f  %d' % (n, Ts.min(), frz.strftime('%Y-%m-%d') if frz else 'never', mld.max(), hc[-1]/1e9, rt, rs, len(used)))
    return out
end_full = R['ADC 5m']['t'][-1]; end_2m = R['ADC 2m']['t'][-1] - dt.timedelta(hours=1)
full = table([n for n in R if not n.startswith('ADC 2m')], end_full, 'FULL YEAR (ADC 5 and 10 m)')
com = table(list(R), end_2m, 'COMMON WINDOW incl. ADC 2 m')
print('\n# resolution sensitivity: spread across grids')
for s in ['GOTM keps', 'GOTM kpp', 'ADC']:
    g = [com['%s %dm' % (s, dz)] for dz in [2, 5, 10]]; f5, f10 = full['%s 5m' % s], full['%s 10m' % s]
    fr = [x['frz'] for x in g]
    print('%-10s 2/5/10 m to 06-22: minTs spread %.3f C, MLDmax spread %4.1f m, RMSE_T spread %.3f, RMSE_S spread %.3f, freeze dates %s'
          % (s, np.ptp([x['minTs'] for x in g]), np.ptp([x['mld'] for x in g]), np.ptp([x['rt'] for x in g]), np.ptp([x['rs'] for x in g]),
             [x.strftime('%m-%d') if x else 'never' for x in fr]))
    print('%-10s 5 vs 10 m full year: minTs diff %.3f C, MLDmax diff %4.1f m, freeze %s vs %s'
          % ('', abs(f5['minTs'] - f10['minTs']), abs(f5['mld'] - f10['mld']),
             f5['frz'].strftime('%m-%d') if f5['frz'] else 'never', f10['frz'].strftime('%m-%d') if f10['frz'] else 'never'))
