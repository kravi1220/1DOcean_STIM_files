"""Paper Figures 2-6 (files K1-K5) and key numbers for the KG6 study (GOTM+STIM vs ADC, Oceananigans LES and Icepack).
Conventions of kg6_common.py (Conservative Temperature and practical salinity everywhere); salt-fixed ADC.
Needs results/tables/blockA_score.txt (blockA_score.py). Writes results/figures/K1-K5*.png and
results/tables/K_numbers.json.
AI use: written and run with Claude (Anthropic); the author must check figures and numbers.
Usage: python3 fig_main.py"""
import os, sys, json, datetime as dt
import numpy as np, netCDF4 as nc, gsw
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kg6_common as A, paths as P
K = P.GOTM.rstrip('/'); MC = P.MPAS; LES = P.LES; FIG = P.FIG; TAB = P.TAB
os.makedirs(FIG, exist_ok=True); os.makedirs(TAB, exist_ok=True)
SZ = dict(K1=(10, 10.5), K2=(10, 7.6), K3=(10, 10.5), K4=(11, 10.5), K5=(11, 10.5))
T0 = dt.datetime(2018, 1, 1)
OUT = {}
C = A.ctd()
COL = {'keps': '#1f77b4', 'kpp': '#2ca02c', 'adc': '#d62728', 'les': 'k', 'obs': '#ff7f0e'}


def mp(path):
    M = A.mpas(path)
    return {k: (v[:-1] if isinstance(v, np.ndarray) and v.shape[:1] == M['t'].shape else v) for k, v in M.items()}


def ice(path):
    a = np.loadtxt(path + '/icepack_column.dat', comments='#')
    return np.array([T0 + dt.timedelta(hours=float(h)) for h in a[:, 0]]), a


def winton_h(path):
    f = nc.Dataset(path + '/all.nc')
    t = np.array([x.replace(tzinfo=None) for x in nc.num2date(f['time'][:], f['time'].units, only_use_cftime_datetimes=False)])
    for v in ('Hice', 'hice', 'Hi'):
        if v in f.variables:
            return t, np.array(f[v][:]).squeeze()
    raise KeyError('no ice thickness in ' + path)


obs_cast = [(t, np.interp(2.0, C[1][~np.isnan(C[2][i])], C[2][i][~np.isnan(C[2][i])])) for i, t in enumerate(C[0])]
t_end = dt.datetime(2019, 8, 10, 11)

# ---------------------------------------------------------------- Block R / C numbers (corrected)
runs = {'R1 k-eps+Winton': K + '/blockR/dt112.5/R1', 'R2 KPP+Winton': K + '/blockR/dt112.5/R2',
        'R3 k-eps no ice': K + '/blockR/dt112.5/R3', 'R4 KPP no ice': K + '/blockR/dt112.5/R4',
        'C1 k-eps+Icepack BL99': K + '/blockC/C1_bl99', 'C2 k-eps+Icepack mushy': K + '/blockC/C2_mushy'}
RC = {}
for n, p in runs.items():
    M = A.gotm(p)
    te = dt.datetime(2019, 4, 4) if n.startswith('C1') else t_end
    rt, rs, used = A.rmse(A.window(M, te), C, te)
    mld, _ = A.derived(M)
    RC[n] = dict(minTs=float(M['T'][:, 0].min()), rmseT=float(rt), rmseS=float(rs), casts=used, mldmax=float(mld.max()),
                 Tsurf_aug=float(M['T'][-1, 0]))
OUT['blockRC'] = RC

# ---------------------------------------------------------------- figure style (v2, 2026-10-07)
# bold, larger text; every legend names the model that produced the line:
#   GOTM (k-ε, KPP; with STIM-Winton or Icepack coupled inside GOTM),
#   MPAS-ADC (Garanaik et al. 2024 closure, run in MPAS-Ocean, a separate model),
#   Oceananigans LES (separate 3D model)
FS = 11
plt.rcParams.update({'font.size': FS, 'font.weight': 'bold', 'axes.labelweight': 'bold', 'axes.titleweight': 'bold',
                     'axes.labelsize': FS + 0.5, 'axes.titlesize': FS + 0.5, 'xtick.labelsize': FS - 1, 'ytick.labelsize': FS - 1,
                     'legend.fontsize': FS - 2, 'axes.linewidth': 1.1, 'lines.linewidth': 1.6})
LAB = {'keps': 'GOTM k-ε', 'kpp': 'GOTM KPP', 'adc': 'MPAS-ADC', 'les': 'Oceananigans LES', 'obs': 'CTD observed'}
LEG = dict(prop={'weight': 'bold', 'size': FS - 2}, framealpha=0.9)


def datefmt(a, fmt='%b %Y'):
    import matplotlib.dates as mdates
    a.xaxis.set_major_formatter(mdates.DateFormatter(fmt))
    for l in a.get_xticklabels():
        l.set_rotation(30); l.set_ha('right')


def top(M, v, zmax=2.0):
    """mean over the top zmax m (each model's own layers)."""
    w = M['z'] <= max(zmax, M['z'].min())
    return M[v][:, w].mean(axis=1)


# ---------------------------------------------------------------- K1: Block A full year
G = {c: {dz: A.gotm(K + '/blockA/gotm/%s_dz%d_dt112.5' % (c, dz)) for dz in (2, 5, 10)} for c in ('keps', 'kpp')}
def les_aug(tag='2m'):
    f = LES + 'les_kg6_era5_aug06_%s_profiles.nc' % tag
    d = nc.Dataset(f); z = -np.array(d['z_aac'][:]); o = np.argsort(z)
    return dict(t=np.array([dt.datetime(2018, 8, 6, 12) + dt.timedelta(seconds=float(x)) for x in d['time'][:]]), z=z[o],
                T=np.array(d['T'][:])[:, o], S=gsw.SP_from_SA(np.array(d['S'][:])[:, o], gsw.p_from_z(-z[o], A.LAT), A.LON, A.LAT))
# full-year LES runs (ocean only; and coupled to Icepack), merged into Figures 1 and 5
import compare_fullyear as CF
LA = CF.OCEAN['LES alone 2m']                      # 6 Aug 2018 - 10 Aug 2019, duplicates from the restart removed
LIM, (LIt, LIa, LIv, LIs) = CF.ICE['LES x Icepack 2m']
AD = {dz: mp(MC + 'blockA/KG6_dz%d' % dz) for dz in (2, 5, 10)}
obs_t = [t for t, _ in obs_cast if dt.datetime(2018, 8, 1) <= t <= dt.datetime(2019, 8, 15)]
obs_T = [v for t, v in obs_cast if dt.datetime(2018, 8, 1) <= t <= dt.datetime(2019, 8, 15)]
obs_S = [np.interp(2.0, C[1][~np.isnan(C[3][i])], C[3][i][~np.isnan(C[3][i])]) for i, t in enumerate(C[0]) if dt.datetime(2018, 8, 1) <= t <= dt.datetime(2019, 8, 15)]
fig, ax = plt.subplots(3, 1, figsize=SZ['K1'], sharex=True, constrained_layout=True)
DZ1 = ((2, '-'),)   # same-grid (2 m) comparison with the 2 m LES
for key, D in (('keps', G['keps']), ('kpp', G['kpp']), ('adc', AD)):
    for dz, ls in DZ1:
        M = D[dz]; mld, _ = A.derived(M)
        lab = '%s %d m' % (LAB[key], dz) + (' (to 23 Jun)' if key == 'adc' and dz == 2 else '')
        ax[0].plot(M['t'], top(M, 'T', dz), ls, color=COL[key], label=lab)
        ax[1].plot(M['t'], top(M, 'S', dz), ls, color=COL[key])
        ax[2].plot(M['t'], mld, ls, color=COL[key])
# ADC 2 m stops on 23 Jun 2019: also show the full-year ADC 5 m run (dashed)
M = AD[5]; mld5, _ = A.derived(M)
ax[0].plot(M['t'], top(M, 'T', 5), '--', color=COL['adc'], lw=1.4, label=LAB['adc'] + ' 5 m (full year)')
ax[1].plot(M['t'], top(M, 'S', 5), '--', color=COL['adc'], lw=1.4); ax[2].plot(M['t'], mld5, '--', color=COL['adc'], lw=1.4)
mldL, _ = A.derived(LA)
ax[0].plot(LA['t'], top(LA, 'T', 2), color=COL['les'], lw=2.4, label=LAB['les'] + ' 2 m'); ax[1].plot(LA['t'], top(LA, 'S', 2), color=COL['les'], lw=2.4)
ax[2].plot(LA['t'], mldL, color=COL['les'], lw=2.4)
ax[0].axhline(A.TF, color='0.4', lw=1, ls='--')
ax[0].text(dt.datetime(2018, 8, 10), A.TF + 0.25, 'local freezing point', color='0.3', fontsize=FS - 1.5, fontweight='bold')
ax[0].plot(obs_t, obs_T, 'o', color=COL['obs'], ms=8, mec='k', label=LAB['obs'])
ax[1].plot(obs_t, obs_S, 'o', color=COL['obs'], ms=8, mec='k')
ax[0].set_ylabel('Surface temperature (°C)'); ax[1].set_ylabel('Surface salinity'); ax[2].set_ylabel('Mixed-layer depth (m)')
ax[2].set_ylim(120, 0)
h_, l_ = ax[0].get_legend_handles_labels(); fig.legend(h_, l_, ncol=3, loc='outside lower center', **LEG)
for a, s in zip(ax, ['(a) surface temperature', '(b) surface salinity', '(c) mixed-layer depth']):
    a.set_title(s, loc='left')
ax[2].set_xlim(dt.datetime(2018, 8, 1), dt.datetime(2019, 8, 15)); datefmt(ax[2])
fig.savefig(FIG + 'K1_KG6_closures_fullyear.png', dpi=220); plt.close(fig)

# ---------------------------------------------------------------- K2: grid sensitivity
common = dt.datetime(2019, 6, 23, 4)
fig, ax = plt.subplots(2, 2, figsize=SZ['K2'], constrained_layout=True); ax = ax.ravel()
GS = {}
# RMSE against the CTD casts on the common window (3 casts, to 23 Jun 2019, so ADC 2 m is included)
SC = {}
for ln in open(TAB + 'blockA_score.txt').read().split('# COMMON WINDOW')[1].split('\n#')[0].splitlines()[2:]:
    p = ln.split()
    if len(p) >= 8:
        key, dz = (p[1], p[2]) if p[0] == 'GOTM' else ('adc', p[1])
        SC[(key, int(dz.rstrip('m')))] = (float(p[-3]), float(p[-2]))   # RMSE_T, RMSE_S
for key, D in (('keps', G['keps']), ('kpp', G['kpp']), ('adc', AD)):
    mins, mlds = [], []
    for dz in (2, 5, 10):
        M = A.window(D[dz], common); mld, _ = A.derived(M)
        mins.append(float(M['T'][:, 0].min())); mlds.append(float(mld.max()))
    GS[LAB[key]] = dict(minTs=mins, mldmax=mlds)
    ax[0].plot([2, 5, 10], mins, 'o-', color=COL[key], ms=8, label=LAB[key]); ax[1].plot([2, 5, 10], mlds, 'o-', color=COL[key], ms=8)
    ax[2].plot([2, 5, 10], [SC[(key, d)][0] for d in (2, 5, 10)], 'o-', color=COL[key], ms=8); ax[3].plot([2, 5, 10], [SC[(key, d)][1] for d in (2, 5, 10)], 'o-', color=COL[key], ms=8)
ax[0].axhline(A.TF, color='0.4', ls='--', lw=1); ax[0].text(5.3, A.TF + 0.03, 'freezing point', fontsize=FS - 1.5, color='0.3', fontweight='bold')
for a in ax:
    a.set_xlabel('Vertical grid spacing Δz (m)'); a.set_xticks([2, 5, 10])
ax[0].set_ylabel('Min. surface T (°C)'); ax[1].set_ylabel('Max. mixed-layer depth (m)')
ax[2].set_ylabel('RMSE T vs CTD (°C)'); ax[3].set_ylabel('RMSE S vs CTD')
ax[0].legend(loc='upper left', **LEG)
for a_, t_ in zip(ax, ['(a) minimum surface temperature', '(b) maximum mixed-layer depth', '(c) temperature error vs CTD', '(d) salinity error vs CTD']):
    a_.set_title(t_, loc='left')
fig.savefig(FIG + 'K2_KG6_grid_sensitivity.png', dpi=220); plt.close(fig)
OUT['grid_sensitivity_to_2019-06-23'] = GS

# OSI SAF (OSI-450-a + operational continuation) at the grid cell nearest KG6, every 10 days, % -> fraction
import csv as _csv
OSI = [(dt.datetime.fromisoformat(r['date']), float(r['concentration_pct']) / 100) for r in _csv.DictReader(open(P.DATA + 'osisaf_kg6.csv'))]


def osi_twin(a, ms=8):
    # OSI SAF concentration (%) on a right-hand axis of the ice-thickness panel; returns its legend handle
    b = a.twinx()
    h, = b.plot([t for t, _ in OSI], [100 * v for _, v in OSI], 'o', color='k', ms=ms, clip_on=False, zorder=5, label='OSI SAF ice concentration (right axis)')
    b.set_ylim(0, 100); b.set_ylabel('OSI SAF ice concentration (%)')
    for l in b.get_yticklabels(): l.set_fontweight('bold')
    return h

# ---------------------------------------------------------------- K3: sea ice
# every ice scheme in the same 4 panels: (a) ice thickness + OSI SAF, (b) surface T, (c) surface S, (d) MLD
tw, hw = winton_h(K + '/blockR/dt112.5/R1')
ICE = [('GOTM k-ε + Icepack BL99', K + '/blockC/C1_bl99', '#9467bd', '--'),
       ('GOTM k-ε + Icepack mushy', K + '/blockC/C2_mushy', COL['keps'], '-'),
       ('MPAS-ADC + Icepack mushy', MC + 'blockD/KG6_dz5_icepack', COL['adc'], '-')]
R = K + '/blockR/dt112.5/'
K3RUNS = [  # label, ocean dict, (t, h) ice thickness or None, colour, style, width
    ('GOTM k-ε + STIM-Winton', A.gotm(R + 'R1'), (tw, hw), '0.25', '-', 1.8),
    ('GOTM k-ε + STIM-MyLake', A.gotm(R + 'R1_mylake'), winton_h(R + 'R1_mylake'), '#17becf', '-.', 2.2),
    ('GOTM k-ε + STIM-Lebedev (patched)', A.gotm(R + 'R1_lebedev_fix'), winton_h(R + 'R1_lebedev_fix'), '#bcbd22', '--', 2.2),
    ('GOTM k-ε + Icepack BL99', A.gotm(K + '/blockC/C1_bl99'), (lambda r: (r[0], r[1][:, 2]))(ice(ICE[0][1])), '#9467bd', '--', 1.8),
    ('GOTM k-ε + Icepack mushy', A.gotm(K + '/blockC/C2_mushy'), (lambda r: (r[0], r[1][:, 2]))(ice(ICE[1][1])), COL['keps'], '-', 1.8),
    ('MPAS-ADC 5 m + Icepack mushy', CF.ICE['MPAS-ADC 5m x Icepack'][0], (lambda r: (r[0], r[1][:, 2]))(ice(ICE[2][1])), COL['adc'], '-', 1.8),
    (LAB['les'] + ' 2 m + Icepack mushy', LIM, (LIt, LIv), COL['les'], '-', 2.6),
    ('GOTM k-ε, no ice', A.gotm(R + 'R3'), None, '#8c564b', ':', 2.0)]


def k3_figure():
    fig, ax = plt.subplots(4, 1, figsize=(10, 13.5), sharex=True, constrained_layout=True)
    for lab, M, ih, c, ls, lw in K3RUNS:
        lw_ = lw
        if ih is not None:
            ax[0].plot(ih[0], ih[1], ls, color=c, lw=lw_)
        mld, _ = A.derived(M)
        ax[1].plot(M['t'], top(M, 'T', 0), ls, color=c, lw=lw_, label=lab); ax[2].plot(M['t'], top(M, 'S', 0), ls, color=c, lw=lw_)   # shallowest layer, as Table 4
        ax[3].plot(M['t'], mld, ls, color=c, lw=lw_)
    ms = 8
    ax[1].plot(obs_t, obs_T, 'o', color=COL['obs'], ms=ms, mec='k', label=LAB['obs']); ax[2].plot(obs_t, obs_S, 'o', color=COL['obs'], ms=ms, mec='k')
    ax[1].axhline(A.TF, color='0.4', ls='--', lw=1)
    ho = osi_twin(ax[0], ms=8)
    ax[3].invert_yaxis()
    for a_, y_ in zip(ax, ['Ice thickness (m)', 'Surface temperature (°C)', 'Surface salinity', 'MLD (m)']):
        a_.set_ylabel(y_)
    tt = ['(a) ice thickness (left) and OSI SAF ice concentration (right)', '(b) surface ocean temperature', '(c) surface ocean salinity', '(d) mixed-layer depth']
    for a_, t_ in zip(ax, tt):
        a_.set_title(t_, loc='left')
    ax[3].set_xlim(dt.datetime(2018, 8, 1), dt.datetime(2019, 8, 15))
    datefmt(ax[3])
    h1, l1 = ax[1].get_legend_handles_labels()
    fig.legend(h1 + [ho], l1 + ['OSI SAF ice concentration (right axis)'], loc='outside lower center', ncol=2, **LEG)
    fig.savefig(FIG + 'K3_KG6_sea_ice.png', dpi=220); plt.close(fig)


k3_figure()

# ---------------------------------------------------------------- K4: LES windows (T, S, MLD)
zf = np.arange(0.25, 150, 0.5)


def les_win(w, dx):
    d = nc.Dataset(LES + 'les_kg6_%s_dx%.1f_profiles.nc' % (w, dx))
    z = -np.array(d['z_aac'][:]); o = np.argsort(z); SA = np.array(d['S'][:])[:, o]
    t0 = dt.datetime.fromisoformat(nc.Dataset(LES + 'les_kg6_%s.nc' % w).window[:19])
    return dict(t=np.array([t0 + dt.timedelta(seconds=float(s)) for s in d['time'][:]]), z=z[o], T=np.array(d['T'][:])[:, o],
                S=gsw.SP_from_SA(SA, gsw.p_from_z(-z[o], A.LAT), A.LON, A.LAT))


fig, ax = plt.subplots(3, 2, figsize=SZ['K4'], constrained_layout=True)
for j, w in enumerate(['B1', 'B2']):
    L = les_win(w, 2.0); t0 = L['t'][0]
    series = [(LAB['les'] + ' 2 m', L, COL['les'], '-'),
              (LAB['keps'] + ' 2 m', A.gotm(K + '/blockW/gotm/%s_keps_dz2' % w), COL['keps'], '-'),
              (LAB['kpp'] + ' 2 m', A.gotm(K + '/blockW/gotm/%s_kpp_dz2' % w), COL['kpp'], '-'),
              (LAB['adc'] + ' 2 m', mp(MC + 'blockW/%s_dz2' % w), COL['adc'], '-')]
    for lab, M, c, ls in series:
        h = np.array([(x - t0).total_seconds() / 3600 for x in M['t']])
        T10 = np.array([np.interp(zf[zf <= 10], M['z'], M['T'][k]).mean() for k in range(len(h))])
        S10 = np.array([np.interp(zf[zf <= 10], M['z'], M['S'][k]).mean() for k in range(len(h))])
        mld, _ = A.derived(M)
        lw = 2.4 if c == 'k' else 1.6
        ax[0, j].plot(h, T10, ls, color=c, lw=lw, label=lab); ax[1, j].plot(h, S10, ls, color=c, lw=lw); ax[2, j].plot(h, mld, ls, color=c, lw=lw)
    name = {'B1': 'B1, early-winter convection', 'B2': 'B2, cooling to freezing + ice fluxes'}[w]
    ax[0, j].set_title('(%s) %s\n%s + 5 days' % ('ab'[j], name, t0.strftime('%d %b %Y')), loc='left')
    ax[1, j].set_title('(%s)' % 'cd'[j], loc='left'); ax[2, j].set_title('(%s)' % 'ef'[j], loc='left')
    ax[2, j].set_xlabel('Hours from window start'); ax[2, j].invert_yaxis()
    for a in ax[:, j]:
        a.set_xlim(0, 120); a.set_xticks(range(0, 121, 24))
ax[0, 0].set_ylabel('Temperature 0–10 m (°C)'); ax[1, 0].set_ylabel('Salinity 0–10 m'); ax[2, 0].set_ylabel('Mixed-layer depth (m)')
ax[0, 0].legend(loc='lower left', **LEG)
fig.savefig(FIG + 'K4_KG6_LES_windows.png', dpi=220); plt.close(fig)

# ---------------------------------------------------------------- K5: Oct-15 long run (time series + profiles at the cast)
d = nc.Dataset(LES + 'les_kg6_era5_oct15_2m_profiles.nc')
z = -np.array(d['z_aac'][:]); o = np.argsort(z)
tL = [dt.datetime(2018, 10, 15, 18) + dt.timedelta(seconds=float(s)) for s in d['time'][:]]
Lo = dict(t=np.array(tL), z=z[o], T=np.array(d['T'][:])[:, o], S=gsw.SP_from_SA(np.array(d['S'][:])[:, o], gsw.p_from_z(-z[o], A.LAT), A.LON, A.LAT))
ic = [k for k, x in enumerate(C[0]) if x.date() == dt.date(2019, 2, 15)][0]; tc = C[0][ic]
i0 = [k for k, x in enumerate(C[0]) if x.date() == dt.date(2018, 10, 15)][0]
ok = ~np.isnan(C[2][ic]); zo = C[1][ok]
O = [(LAB['les'] + ' 2 m', Lo, COL['les'], '-'), (LAB['keps'] + ' 2 m', A.gotm(K + '/oct15/runs/keps_dz2'), COL['keps'], '-'),
     (LAB['kpp'] + ' 2 m', A.gotm(K + '/oct15/runs/kpp_dz2'), COL['kpp'], '-'),
     (LAB['adc'] + ' 5 m*', mp(MC + 'oct15/KG6oct_dz5'), COL['adc'], '-')]
fig = plt.figure(figsize=SZ['K5'], constrained_layout=True)
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.35])
axT, axS, axM = fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
axPT, axPS = fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[2, 1])
for lab, M, c, ls in O:
    lw = 2.4 if c == 'k' else 1.6
    mld, _ = A.derived(M)
    axT.plot(M['t'], top(M, 'T'), ls, color=c, lw=lw, label=lab); axS.plot(M['t'], top(M, 'S'), ls, color=c, lw=lw)
    axM.plot(M['t'], mld, ls, color=c, lw=lw)
    k = int(np.argmin([abs((x - tc).total_seconds()) for x in M['t']]))
    axPT.plot(M['T'][k], M['z'], ls, color=c, lw=lw); axPS.plot(M['S'][k], M['z'], ls, color=c, lw=lw)
for tt, ii in ((C[0][i0], i0), (tc, ic)):
    okk = ~np.isnan(C[2][ii])
    axT.plot(tt, np.interp(2.0, C[1][okk], C[2][ii][okk]), 'o', color=COL['obs'], ms=9, mec='k')
    axS.plot(tt, np.interp(2.0, C[1][okk], C[3][ii][okk]), 'o', color=COL['obs'], ms=9, mec='k')
axT.plot([], [], 'o', color=COL['obs'], ms=9, mec='k', label=LAB['obs'])
axPT.plot(C[2][ic][ok], zo, color=COL['obs'], lw=2.6); axPS.plot(C[3][ic][ok], zo, color=COL['obs'], lw=2.6)
axT.axhline(A.TF, color='0.4', ls='--', lw=1)
axT.set_ylabel('Surface temperature (°C)'); axS.set_ylabel('Surface salinity'); axM.set_ylabel('Mixed-layer depth (m)'); axM.invert_yaxis()
for a in (axT, axS, axM):
    datefmt(a, '%d %b')
for a in (axPT, axPS):
    a.set_ylim(200, 0); a.set_ylabel('Depth (m)')
axPT.set_xlabel('Temperature (°C)'); axPS.set_xlabel('Salinity')
axT.legend(ncol=3, loc='upper right', **LEG)
for a, s in ((axT, '(a) surface temperature from the 15 Oct 2018 cast'), (axS, '(b) surface salinity'), (axM, '(c) mixed-layer depth'),
             (axPT, '(d) temperature profile, 15 Feb 2019 cast'), (axPS, '(e) salinity profile, 15 Feb 2019 cast')):
    a.set_title(s, loc='left')
fig.savefig(FIG + 'K5_KG6_oct15_LES_comparison.png', dpi=220); plt.close(fig)


json.dump(OUT, open(TAB + 'K_numbers.json', 'w'), indent=1, default=str)
print('figures written')
