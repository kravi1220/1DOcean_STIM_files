"""ADC outer-time-step ladder (5 m, 2018-08-06 12:00 + 1700 h) against the 9.375 s run (Section 4.2.3).
Usage: python3 adc_dt_ladder.py > ../results/tables/adc_dt_ladder.txt  -> SST / MLD differences."""
import numpy as np
import kg6_common as A, paths as P
def load(c):
    M = A.mpas(P.MPAS + c); mld, _ = A.derived(M); return M['t'], M['T'][:, 0], mld
tr, Tr, Mr = load('blockA/KG6_dz5'); ref = {t: (a, b) for t, a, b in zip(tr, Tr, Mr)}
print('%-10s %8s %10s %10s %10s' % ('outer dt', 'hours', 'dSST end', 'max|dSST|', 'dMLD end'))
for lab, c in [('30 min', 'dt30min'), ('15 min', 'dt15min'), ('5 min', 'dt5min'), ('2.5 min', 'dt2p5min'), ('1.25 min', 'dt1p25min'), ('37.5 s', 'dt37p5s'), ('18.75 s', 'dt18p75s')]:
    t, T, M = load('blockA_dt/KG6_dz5_%s' % c)
    k = [i for i, x in enumerate(t) if x in ref]
    dT = np.array([T[i] - ref[t[i]][0] for i in k]); dM = np.array([M[i] - ref[t[i]][1] for i in k])
    print('%-10s %8d %+10.4f %10.4f %+10.1f' % (lab, len(k), dT[-1], np.abs(dT).max(), dM[-1]))
