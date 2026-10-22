"""Paper Figures 7-11 (files K7-K11) for the KG6 study.
  K7_KG6_window_scores.png   5-day LES windows: 0-10 m T, MLD and RMS_T vs grid spacing (from blockW_compare.txt)
  K8_KG6_ice_season.png      ice season of every ice scheme vs OSI SAF, max thickness, Aug surface T
  K9_KG6_interaction.png     closure x ice interaction, Eq. 8 (from blockD_table.txt)
  K10_KG6_budgets.png        column heat and salt budgets: surface input, 0-50 m storage, export below 50 m
  K11_KG6_skill.png          target diagrams of all runs against the CTD casts (T and S, 2-150 m)
Needs results/tables/blockW_compare.txt and blockD_table.txt; numbers are written to results/tables/K_summary_numbers.txt.
Usage: python3 fig_summary.py
"""
import os, re, sys, datetime as dt, numpy as np, netCDF4 as nc
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.dates as mdates
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths as P
import compare_fullyear as CF          # OCEAN / ICE runs, CTD casts
A = CF.A
K = P.GOTM.rstrip('/'); MC = P.MPAS; FIG = P.FIG
os.makedirs(FIG, exist_ok=True)
FS = 11.5
plt.rcParams.update({'font.size': FS, 'font.weight': 'bold', 'axes.labelweight': 'bold', 'axes.titleweight': 'bold',
                     'axes.linewidth': 1.2, 'lines.linewidth': 2.0})
LEG = dict(prop={'weight': 'bold', 'size': FS - 2}, framealpha=0.9)
COL = {'keps': '#1f77b4', 'kpp': '#2ca02c', 'adc': '#d62728', 'les': 'k', 'obs': '#ff7f0e'}
LAB = {'keps': 'GOTM k-ε', 'kpp': 'GOTM KPP', 'adc': 'MPAS-ADC'}
OUT = open(P.TAB + 'K_summary_numbers.txt', 'w')
def log(*a):
    print(*a); print(*a, file=OUT)
RHOCP = A.RHO0 * A.CP


# =============================================================== K7: 5-day windows (from blockW_compare.txt)
def window_scores():
    """{window: {run: (T0-10 degC, MLD m, RMS_T degC vs LES 2 m)}} at 120 h, read from blockW_compare.py's output."""
    W, w = {}, None
    for ln in open(P.TAB + 'blockW_compare.txt'):
        m = re.match(r'# (B[12]):', ln)
        if m:
            w = m.group(1); W[w] = {}; continue
        p = ln.split('|')
        if w is None or len(p) < 3 or p[0].split()[-1] != '120':
            continue
        name = ' '.join(p[0].split()[:-1]); v = p[1].split()
        W[w][name] = (float(v[0]), float(v[2]), float(v[3]))
    return {b: {'les2': R['LES 2m'], 'les4': R['LES 4m'],
                'keps': [R['GOTM keps %dm' % d] for d in (2, 5, 10)], 'kpp': [R['GOTM kpp %dm' % d] for d in (2, 5, 10)],
                'adc': [R['ADC %dm' % d] for d in (2, 5, 10)]} for b, R in W.items()}
WIN = window_scores()
fig, ax = plt.subplots(2, 3, figsize=(11, 6.6), constrained_layout=True)
ttl = {'B1': 'B1 (convection)', 'B2': 'B2 (freezing)'}
DZ = np.array([2, 5, 10]); off = {'keps': -0.25, 'kpp': 0, 'adc': 0.25}
for i, b in enumerate(('B1', 'B2')):
    for j, (q, yl, nm) in enumerate([(0, '0–10 m T (°C)', '0–10 m T'), (1, 'MLD (m)', 'MLD'), (2, 'RMS T vs LES 2 m (°C)', 'RMS T vs LES')]):
        a = ax[i, j]
        if q < 2:
            a.axhline(WIN[b]['les2'][q], color='k', lw=2.6, label='LES 2 m (reference)')
            a.axhline(WIN[b]['les4'][q], color='k', lw=1.4, ls='--', label='LES 4 m')
        else:
            a.axhspan(0, WIN[b]['les4'][2], color='0.85', label='LES 4 m vs 2 m')
        for c in ('keps', 'kpp', 'adc'):
            a.plot(DZ + off[c] * (DZ / 5), [v[q] for v in WIN[b][c]], 'o-', color=COL[c], ms=8, label=LAB[c])
        a.set_xticks(DZ); a.set_xlabel('Δz (m)'); a.set_ylabel(yl)
        a.set_title('(%s) %s: %s' % ('abcdef'[3 * i + j], ttl[b], nm), loc='left', fontsize=FS - 1)
        if q == 1:
            a.invert_yaxis()
ax[0, 0].legend(loc='best', **LEG)
fig.savefig(FIG + 'K7_KG6_window_scores.png', dpi=220); plt.close(fig)


# =============================================================== K8: ice season of all schemes vs OSI SAF
def stim(run):
    f = nc.Dataset(K + '/blockR/dt112.5/' + run + '/all.nc')
    t0 = dt.datetime.fromisoformat(f['time'].units.split('since ')[1].strip())
    return np.array([t0 + dt.timedelta(seconds=float(s)) for s in f['time'][:]]), np.array(f['Hice'][:]).squeeze()

def ipk(path):
    t, ai, vi, vs = CF.ice_gotm_style(path + '/icepack_column.dat'); return t, vi

def T2(M, day):                                     # 2 m temperature nearest to a date (None if no output)
    k = CF.at(M, day)
    return None if k is None else float(np.interp(2, M['z'], M['T'][k]))

AUG = dt.datetime(2019, 8, 10, 12)
LIM, LI = CF.ICE['LES x Icepack 2m']
ROWS = [  # label, (t, h), ocean dict or None, colour
 ('k-ε + STIM-Winton', stim('R1'), A.gotm(K + '/blockR/dt112.5/R1'), '0.25'),
 ('k-ε + STIM-MyLake', stim('R1_mylake'), A.gotm(K + '/blockR/dt112.5/R1_mylake'), '#17becf'),
 ('k-ε + STIM-Lebedev*', stim('R1_lebedev_fix'), A.gotm(K + '/blockR/dt112.5/R1_lebedev_fix'), '#bcbd22'),
 ('k-ε + Icepack BL99†', ipk(K + '/blockC/C1_bl99'), None, '#9467bd'),
 ('k-ε + Icepack mushy', ipk(K + '/blockC/C2_mushy'), A.gotm(K + '/blockC/C2_mushy'), COL['keps']),
 ('k-ε 2 m + Icepack mushy', ipk(K + '/blockD/keps_dz2_icepack'), CF.ICE['GOTM k-eps 2m x Icepack'][0], COL['keps']),
 ('ADC 5 m + Icepack mushy', ipk(MC + 'blockD/KG6_dz5_icepack'), CF.ICE['MPAS-ADC 5m x Icepack'][0], COL['adc']),
 ('LES 2 m + Icepack mushy', (LI[0], LI[2]), LIM, 'k')]
z = np.load(P.DATA + 'osisaf_region_20181201_20190531.npz')
kc = np.unravel_index(np.argmin((z['lat'] - 67.447)**2 + (z['lon'] + 23.875)**2), z['lat'].shape)
od = np.array([dt.datetime.fromisoformat(s) + dt.timedelta(hours=12) for s in z['dates']]); oc = z['conc'][:, kc[0], kc[1]]
log('# K8 ice season (h > 1 cm), max grid-mean thickness, 2 m T on 10 Aug 2019 (CTD 3.34 C)')
fig = plt.figure(figsize=(11, 6.4), constrained_layout=True)
gs = fig.add_gridspec(1, 3, width_ratios=[3.2, 1, 1])
a0, a1, a2 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
cm = plt.get_cmap('Blues'); norm = matplotlib.colors.Normalize(0, 0.3)
names = [r[0] for r in ROWS] + ['OSI SAF (daily, KG6 cell)']
y = np.arange(len(names))[::-1]
for (lab, (t, h), M, c), yy in zip(ROWS, y[:-1]):
    t = np.asarray(t); h = np.asarray(h); on = h > 0.01
    for k0 in np.where(on & ~np.r_[False, on[:-1]])[0]:                       # contiguous ice segments
        k1 = k0 + np.argmax(~on[k0:]) if (~on[k0:]).any() else len(on)
        seg = slice(k0, k1)
        a0.scatter(t[seg], np.full(k1 - k0, yy), c=np.clip(h[seg], 0, 0.3), cmap=cm, norm=norm, marker='|', s=420, lw=2.5)
    a0.plot([t[0], t[-1]], [yy, yy], color='0.85', lw=0.8, zorder=0)
    i = np.where(on)[0]; Tend = T2(M, AUG) if M is not None else None
    a1.barh(yy, h.max(), color=c, height=0.6); a1.text(h.max() + 0.02, yy, '%.2f' % h.max(), va='center', fontsize=FS - 3)
    if Tend is not None:
        a2.barh(yy, Tend, color=c, height=0.6); a2.text(Tend + 0.2, yy, '%.1f' % Tend, va='center', fontsize=FS - 3)
    log('  %-24s %s - %s  (%5.1f d)  Hmax %.3f m  T2m 10 Aug %s' % (lab, t[i[0]].strftime('%d %b'), t[i[-1]].strftime('%d %b'),
        (t[i[-1]] - t[i[0]]).total_seconds() / 86400, h.max(), 'n/a' if Tend is None else '%.2f' % Tend))
a0.scatter(od[oc >= 15], np.full((oc >= 15).sum(), y[-1]), c=oc[oc >= 15], cmap='Reds', vmin=0, vmax=50, marker='s', s=40)
a0.scatter(od[(oc > 0) & (oc < 15)], np.full(((oc > 0) & (oc < 15)).sum(), y[-1]), color='#f4a582', marker='s', s=12)
log('  OSI SAF KG6 cell: %d of %d days with ice >= 15 %% (Dec-May); max %.0f %%' % ((oc >= 15).sum(), len(oc), np.nanmax(oc)))
a0.set_yticks(y); a0.set_yticklabels(names); a0.set_xlim(dt.datetime(2018, 12, 1), dt.datetime(2019, 5, 25))
a0.xaxis.set_major_formatter(mdates.DateFormatter('%b')); a0.set_title('(a) ice present (colour: thickness, m)', loc='left')
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cm); cb = fig.colorbar(sm, ax=a0, orientation='horizontal', shrink=0.6, pad=0.02)
cb.set_label('Grid-mean ice thickness (m); Winton exceeds 0.3 m')
a1.set_yticks(y); a1.set_yticklabels([]); a1.set_xlabel('Max thickness (m)'); a1.set_title('(b)', loc='left')
a2.set_yticks(y); a2.set_yticklabels([]); a2.axvline(3.34, color=COL['obs'], lw=2.5, label='CTD 10 Aug'); a2.set_xlabel('2 m T, 10 Aug (°C)')
a2.set_title('(c)', loc='left'); a2.legend(loc='lower right', **LEG)
for a in (a1, a2):
    a.set_ylim(a0.get_ylim())
fig.savefig(FIG + 'K8_KG6_ice_season.png', dpi=220); plt.close(fig)


# =============================================================== K9: closure x ice interaction (from blockD_table.txt)
BD = open(P.TAB + 'blockD_table.txt').read()
def effect(v):     # (k-eps, ADC) ice effect of quantity v
    m = re.search(r'(?m)^%s\s+ice effect: keps\s+(\S+)\s+ADC\s+(\S+)' % v, BD); return float(m.group(1)), float(m.group(2))
DAYS = [float(re.search(r'(?m)^%s ice .*\|\s+\S+\s+\S+\s+(\S+)' % c, BD).group(1)) for c in ('keps', 'ADC')]
D = {'min $T_s$ (°C)': effect('minTs'), 'max MLD (m)': effect('mld'), 'RMSE T (°C)': effect('rt')}
fig, ax = plt.subplots(1, 4, figsize=(11, 3.8), constrained_layout=True)
for a, (q, (ke, ad)) in zip(ax[:3], D.items()):
    v = [ke, ad, ad - ke]
    a.bar([0, 1, 2], v, color=[COL['keps'], COL['adc'], '0.45'], width=0.65)
    a.axhline(0, color='k', lw=1); a.set_xticks([0, 1, 2]); a.set_xticklabels(['k-ε', 'ADC', 'inter-\naction'])
    a.set_title('(%s) ice effect on\n%s' % ('abc'[list(D).index(q)], q.replace('$', '').replace('_s', 's')), loc='left', fontsize=FS - 1)
a = ax[3]
a.bar([0, 1], DAYS, color=[COL['keps'], COL['adc']], width=0.65); a.set_ylim(70, 80)
a.set_xticks([0, 1]); a.set_xticklabels(['k-ε', 'ADC']); a.set_title('(d) ice days\n(9 Feb – 26/27 Apr)', loc='left', fontsize=FS - 1)
fig.savefig(FIG + 'K9_KG6_interaction.png', dpi=220); plt.close(fig)


# =============================================================== K10: heat and salt budgets (closed columns)
def budget(M, h=50.0):
    z = M['z']; dz = np.gradient(z) if np.ptp(np.diff(z)) > 1e-6 else np.full(len(z), z[1] - z[0])
    up = z <= h
    HT = RHOCP * (M['T'] * dz).sum(1) / 1e6; HU = RHOCP * (M['T'][:, up] * dz[up]).sum(1) / 1e6   # MJ m-2
    ST = A.RHO0 * (M['S'] * dz).sum(1) / 1e3; SU = A.RHO0 * (M['S'][:, up] * dz[up]).sum(1) / 1e3  # kg salt m-2
    q = HT - HT[0]; u = HU - HU[0]; s = ST - ST[0]; su = SU - SU[0]
    return M['t'], q, u, q - u, s, su, s - su
BR = [('Oceananigans LES 2 m', CF.OCEAN['LES alone 2m'], 'k', 2.6), ('GOTM k-ε 2 m', CF.OCEAN['GOTM k-eps 2m'], COL['keps'], 1.8),
      ('GOTM KPP 2 m', CF.OCEAN['GOTM KPP 2m'], COL['kpp'], 1.8), ('MPAS-ADC 5 m', CF.OCEAN['MPAS-ADC 5m'], COL['adc'], 1.8)]
fig, ax = plt.subplots(2, 3, figsize=(11, 7.0), sharex=True, constrained_layout=True)
log('\n# K10 budgets relative to 6 Aug 2018 (heat MJ m-2, salt kg m-2); layer 0-50 m; below = surface - storage(0-50)')
for lab, M, c, lw in BR:
    t, q, u, d, s, su, sd = budget(M)
    for a, v in zip(ax.ravel(), (q, u, d, s, su, sd)):
        a.plot(t, v, color=c, lw=lw, label=lab)
    for day in (dt.datetime(2019, 2, 15, 12), dt.datetime(2019, 8, 10, 0)):
        k = int(np.argmin(np.abs([(x - day).total_seconds() for x in t])))
        log('  %-22s %s  Qsurf %8.0f  dH(0-50) %7.0f  below50 %8.0f | Fsalt %6.1f  dS(0-50) %6.1f  below50 %6.1f'
            % (lab, day.strftime('%d %b %Y'), q[k], u[k], d[k], s[k], su[k], sd[k]))
ttl = ['(a) cumulative surface heat input', '(b) heat stored in 0–50 m', '(c) heat moved below 50 m',
       '(d) cumulative surface salt input', '(e) salt stored in 0–50 m', '(f) salt moved below 50 m']
for a, t_ in zip(ax.ravel(), ttl):
    a.set_title(t_, loc='left', fontsize=FS - 1); a.axhline(0, color='0.5', lw=0.8)
    a.xaxis.set_major_formatter(mdates.DateFormatter('%b')); a.xaxis.set_major_locator(mdates.MonthLocator(bymonth=(8, 10, 12, 2, 4, 6)))
ax[0, 0].set_ylabel('Heat (MJ m$^{-2}$)'); ax[1, 0].set_ylabel('Salt (kg m$^{-2}$)'); ax[0, 0].legend(loc='lower left', **LEG)
fig.savefig(FIG + 'K10_KG6_budgets.png', dpi=220); plt.close(fig)


# =============================================================== K11: target diagrams vs CTD (2-150 m, casts after the start)
C = CF.C; ZC = np.arange(2.0, 151.0, 1.0)
def stats(M, var):
    m_all, o_all = [], []
    for day in CF.casts:
        k = CF.at(M, day)
        if k is None:
            continue
        i = C[0].index(day); ok = ~np.isnan(C[2 if var == 'T' else 3][i])
        o = np.interp(ZC, C[1][ok], C[2 if var == 'T' else 3][i][ok]); m = np.interp(ZC, M['z'], M[var][k])
        m_all.append(m); o_all.append(o)
    m, o = np.concatenate(m_all), np.concatenate(o_all)
    bias = (m - o).mean(); crmsd = np.sqrt((((m - m.mean()) - (o - o.mean()))**2).mean())
    return bias, np.sign(m.std() - o.std()) * crmsd, len(m_all)
G = K + '/blockA/gotm/'
TR = [(LAB['keps'], COL['keps'], 'o', {d: A.gotm(G + 'keps_dz%d_dt112.5' % d) for d in (2, 5, 10)}),
      (LAB['kpp'], COL['kpp'], 's', {d: A.gotm(G + 'kpp_dz%d_dt112.5' % d) for d in (2, 5, 10)}),
      (LAB['adc'], COL['adc'], '^', {2: CF.OCEAN['MPAS-ADC 2m*'], 5: CF.OCEAN['MPAS-ADC 5m'], 10: CF.mp('blockA/KG6_dz10')})]
EXTRA = [('LES 2 m', CF.OCEAN['LES alone 2m'], 'k', '*', 18), ('LES 2 m + Icepack', LIM, 'k', 'P', 12),
         ('k-ε + Winton', A.gotm(K + '/blockR/dt112.5/R1'), '0.45', 'D', 10), ('k-ε + Icepack mushy', A.gotm(K + '/blockC/C2_mushy'), COL['keps'], 'D', 10)]
fig, ax = plt.subplots(1, 2, figsize=(11, 5.4), constrained_layout=True)
log('\n# K11 skill vs CTD casts (2-150 m): bias, signed centred RMSD, total RMSD, n casts')
for j, var in enumerate(('T', 'S')):
    a = ax[j]; PTS = []
    for lab, c, mk, R in TR:
        for d, M in R.items():
            b, cr, n = stats(M, var); size = {2: 13, 5: 9, 10: 6}[d]
            a.plot(cr, b, mk, color=c, ms=size, mec='k', label=lab if d == 2 else None); PTS.append((cr, b))
            a.annotate('%d' % d, (cr, b), xytext=(6, -10 if mk == 's' else 4), textcoords='offset points', fontsize=FS - 4)
            log('  %s %-10s %2d m  bias %+.3f  cRMSD %+.3f  RMSD %.3f  n %d' % (var, lab, d, b, cr, np.hypot(b, cr), n))
    for lab, M, c, mk, ms in EXTRA:
        b, cr, n = stats(M, var)
        a.plot(cr, b, mk, color=c, ms=ms, mec='k', label=lab); PTS.append((cr, b))
        log('  %s %-20s bias %+.3f  cRMSD %+.3f  RMSD %.3f  n %d' % (var, lab, b, cr, np.hypot(b, cr), n))
    P_ = np.array(PTS); cx, cy = (P_.min(0) + P_.max(0)) / 2; half = 0.5 * max(np.ptp(P_[:, 0]), np.ptp(P_[:, 1])) + (0.1 if var == 'T' else 0.04)
    th = np.linspace(0, 2 * np.pi, 600); step = 0.1 if var == 'T' else 0.05
    for r in np.arange(step, 3.0, step):     # total RMSD = distance from a perfect model (origin, outside the zoom)
        a.plot(r * np.cos(th), r * np.sin(th), ':', color='0.7', lw=0.9, zorder=0)
    a.set_xlim(cx - half, cx + half); a.set_ylim(cy - half, cy + half); a.set_aspect('equal')
    a.set_xlabel('± centred RMSD (%s)' % ('°C' if var == 'T' else 'g/kg'))
    a.set_ylabel('Bias, model − CTD (%s)' % ('°C' if var == 'T' else 'g/kg'))
    a.set_title('(%s) %s (dotted: total RMSD, every %s)' % ('ab'[j], 'temperature' if var == 'T' else 'salinity', '0.1 °C' if var == 'T' else '0.05'), loc='left', fontsize=FS - 1)
ax[1].legend(loc='center left', bbox_to_anchor=(1.02, 0.5), **LEG)
fig.savefig(FIG + 'K11_KG6_skill.png', dpi=220); plt.close(fig)
OUT.close()
print('figures written')
