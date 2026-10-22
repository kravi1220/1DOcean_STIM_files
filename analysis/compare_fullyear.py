"""Full-year 2 m comparison (6 Aug 2018 -> 10 Aug 2019, same cast and ERA5 forcing):
 (A) ocean only:   Oceananigans LES alone  vs  GOTM k-eps 2 m, GOTM KPP 2 m, MPAS-ADC 2 m (to 23 Jun) / 5 m
 (B) with Icepack: Oceananigans LES x Icepack  vs  GOTM k-eps 2 m x Icepack, MPAS-ADC 5 m x Icepack
All Icepack runs use the same mushy-layer configuration (Sections 4.3.3, 4.6; also imported by the figure scripts).
Usage: python3 compare_fullyear.py > ../results/tables/compare_fullyear.txt"""
import os, sys, glob, datetime as dt, numpy as np, netCDF4 as nc, gsw
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kg6_common as A, paths as P
L, MC, K = P.LES, P.MPAS, P.GOTM
T0 = dt.datetime(2018, 8, 6, 12); END = dt.datetime(2019, 8, 10, 12)


def les(files):
    """concatenate LES profile files (restarts/segments), keep the first record of duplicated times"""
    t, T, S = [], [], []
    for f in files:
        d = nc.Dataset(f)
        z = -np.array(d['z_aac'][:]); o = np.argsort(z)
        t += [T0 + dt.timedelta(seconds=float(x)) for x in d['time'][:]]
        T.append(np.array(d['T'][:])[:, o]); S.append(gsw.SP_from_SA(np.array(d['S'][:])[:, o], gsw.p_from_z(-z[o], A.LAT), A.LON, A.LAT))
    t = np.array(t); T = np.concatenate(T); S = np.concatenate(S)
    good = ~np.all(T == 0, axis=1)      # drop empty records (a record left empty when a run is interrupted while writing)
    t, T, S = t[good], T[good], S[good]
    t = np.array([x.replace(microsecond=0) + dt.timedelta(seconds=round(x.microsecond / 1e6)) for x in t])   # restart overlaps: match to the second
    _, keep = np.unique(t, return_index=True)
    return dict(t=t[keep], z=z[o], T=T[keep], S=S[keep])


def mp(c):
    M = A.mpas(MC + c); return {k: (v[:-1] if isinstance(v, np.ndarray) and v.shape[:1] == M['t'].shape else v) for k, v in M.items()}


def ice_gotm_style(path):     # icepack_column.dat (GOTM / MPAS couplings): hours since 1 Jan 2018
    a = np.loadtxt(path, comments='#')
    return np.array([dt.datetime(2018, 1, 1) + dt.timedelta(hours=float(h)) for h in a[:, 0]]), a[:, 1], a[:, 2], a[:, 3]


def ice_les(tags):
    import csv
    t, ai, vi, vs = [], [], [], []
    for tg in tags:
        f = L + 'les_kg6_icepack_%s_ice.csv' % tg
        if not os.path.exists(f):
            continue
        for r in csv.DictReader(open(f)):
            t.append(dt.datetime.fromisoformat(r['time'][:19])); ai.append(float(r['aice'])); vi.append(float(r['vice'])); vs.append(float(r['vsno']))
    t = np.array(t); _, k = np.unique(t, return_index=True)
    return t[k], np.array(ai)[k], np.array(vi)[k], np.array(vs)[k]


C = A.ctd()
casts = [x for x in C[0] if T0 <= x <= END + dt.timedelta(days=1)]
G = K + 'blockA/gotm/'
OCEAN = {'LES alone 2m': les([L + 'les_kg6_era5_aug06_2m_profiles.nc']),
         'GOTM k-eps 2m': A.gotm(G + 'keps_dz2_dt112.5'), 'GOTM KPP 2m': A.gotm(G + 'kpp_dz2_dt112.5'),
         'MPAS-ADC 2m*': mp('blockA/KG6_dz2'), 'MPAS-ADC 5m': mp('blockA/KG6_dz5')}
segs = ['cpl2m_seg1', 'cpl2m_seg1b', 'cpl2m_seg2', 'cpl2m_seg3']   # job segments of the coupled run (experiments/les/README)
ICE = {'LES x Icepack 2m': (les(sorted(f for tg in segs for f in glob.glob(L + 'les_kg6_icepack_%s_profiles.nc' % tg))), ice_les(segs)),
       'GOTM k-eps 2m x Icepack': (A.gotm(K + 'blockD/keps_dz2_icepack'), ice_gotm_style(K + 'blockD/keps_dz2_icepack/icepack_column.dat')),
       'MPAS-ADC 5m x Icepack': (mp('blockD/KG6_dz5_icepack'), ice_gotm_style(MC + 'blockD/KG6_dz5_icepack/icepack_column.dat'))}
zf = np.arange(2.0, 150, 1.0)


def at(M, day):
    k = int(np.argmin([abs((x - day).total_seconds()) for x in M['t']]))
    return k if abs((M["t"][k] - day).total_seconds()) < 12 * 3600 else None   # last cast is 9.4 h after the run end


def table(R, ref, title):
    print('\n# %s; reference %s; T, S at 2 m; RMS_T vs reference over 2-150 m' % (title, ref))
    for day in casts:
        i = C[0].index(day); ok = ~np.isnan(C[2][i])
        print('\n%s  CTD %7.3f %7.3f' % (day.strftime('%Y-%m-%d'), np.interp(2, C[1][ok], C[2][i][ok]), np.interp(2, C[1][ok], C[3][i][ok])))
        Mr = R[ref][0] if isinstance(R[ref], tuple) else R[ref]; kr = at(Mr, day)
        for n, v in R.items():
            M = v[0] if isinstance(v, tuple) else v; k = at(M, day)
            if k is None:
                print('  %-24s  (no output at this date)' % n); continue
            mld, _ = A.derived({**M, 't': M['t'][k:k + 1], 'T': M['T'][k:k + 1], 'S': M['S'][k:k + 1]})
            rms = np.sqrt(np.mean((np.interp(zf, M['z'], M['T'][k]) - np.interp(zf, Mr['z'], Mr['T'][kr]))**2)) if kr is not None else np.nan
            ro = A.rmse({**M, 't': M['t'][k:k + 1], 'T': M['T'][k:k + 1], 'S': M['S'][k:k + 1]}, C, day + dt.timedelta(days=1))
            print('  %-24s %7.3f %7.3f  MLD %6.1f  RMS_T vs ref %6.3f' % (n, np.interp(2, M['z'], M['T'][k]), np.interp(2, M['z'], M['S'][k]), mld[0], rms))
    print('\n# season summary (to %s)' % END.date())
    for n, v in R.items():
        M = v[0] if isinstance(v, tuple) else v; w = M['t'] <= END; mld, _ = A.derived({**M, 't': M['t'][w], 'T': M['T'][w], 'S': M['S'][w]})
        Ts = M['T'][w, 0]; i = np.argmin(Ts); line = '  %-24s to %s  min T_s %7.3f (%s)  max MLD %6.1f' % (n, M['t'][w][-1].strftime('%Y-%m-%d'), Ts[i], M['t'][w][i].strftime('%m-%d'), mld.max())
        if isinstance(v, tuple):
            t, ai, vi, vs = v[1]; on = np.where(ai > 0.01)[0]
            line += ('  ice %s - %s (%.1f d)  max vice %.3f m  max aice %.2f  max snow %.3f m' % (t[on[0]].strftime('%m-%d'), t[on[-1]].strftime('%m-%d'), len(on) / 24, vi.max(), ai.max(), vs.max())) if on.size else '  no ice'
        print(line)


if __name__ == '__main__':
    table(OCEAN, 'LES alone 2m', '(A) OCEAN ONLY, no ice model in any run')
    table(ICE, 'LES x Icepack 2m', '(B) WITH ICEPACK (mushy) in every run')
