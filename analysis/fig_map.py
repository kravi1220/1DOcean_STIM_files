"""Study-area map for the KG6 paper: OSI SAF sea-ice concentration and the KG6 CTD station.
 (a) concentration on 1 Feb 2019 (the date the ice edge reached KG6, 26 % at the station)
 (b) ice occurrence: % of days with concentration >= 15 %, 1 Dec 2018 - 31 May 2019
Input: data/osisaf_region_20181201_20190531.npz (OSI-450-a daily fields, 1 Dec 2018 - 31 May 2019)
Output: results/figures/K0_KG6_osisaf_map.png (paper Figure 1). Needs cartopy.
AI use: written with Claude (Anthropic); the author must check it.
"""
import os, sys, datetime as dt, numpy as np, matplotlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths as P
matplotlib.use('Agg')
import matplotlib.pyplot as plt, cartopy.crs as ccrs, cartopy.feature as cfeature
FIG = P.FIG; os.makedirs(FIG, exist_ok=True)
FS = 12
plt.rcParams.update({'font.size': FS, 'font.weight': 'bold', 'axes.labelweight': 'bold', 'axes.titleweight': 'bold'})
KG6 = (67.58005324675324, -23.92983766233766)

z = np.load(P.DATA + 'osisaf_region_20181201_20190531.npz')
lat, lon, dates, conc = z['lat'], z['lon'], [dt.date.fromisoformat(s) for s in z['dates']], z['conc']
k = dates.index(dt.date(2019, 2, 1))
freq = 100 * np.nanmean(conc >= 15, axis=0); freq[np.all(np.isnan(conc), axis=0)] = np.nan
print('days', len(dates), 'KG6 cell freq %.1f %%' % freq[np.unravel_index(np.argmin((lat - 67.447)**2 + (lon + 23.875)**2), lat.shape)])

proj = ccrs.LambertConformal(central_longitude=-24, central_latitude=67, standard_parallels=(64, 70))
pc = ccrs.PlateCarree()
fig, ax = plt.subplots(1, 2, figsize=(12, 5.4), subplot_kw=dict(projection=proj), constrained_layout=True)
panels = [(conc[k], 'Blues', 'Ice concentration (%)', '(a) OSI SAF ice concentration, 1 Feb 2019'),
          (freq, 'PuBu', 'Days with ice ≥ 15 % (%)', '(b) Ice occurrence, 1 Dec 2018 – 31 May 2019')]
for a, (fld, cm, cl, ttl) in zip(ax, panels):
    a.set_extent([-38, -12, 63.2, 71.2], crs=pc)
    a.add_feature(cfeature.OCEAN, facecolor='#f4f8fb', zorder=0)
    pm = a.pcolormesh(lon, lat, np.ma.masked_where(~np.isfinite(fld) | (fld <= 0), fld), cmap=cm, vmin=0, vmax=100, transform=pc, shading='nearest', zorder=1)
    a.add_feature(cfeature.LAND, facecolor='#ddd6c9', zorder=2); a.add_feature(cfeature.COASTLINE, lw=0.9, zorder=3)
    if fld is conc[k]:
        a.contour(lon, lat, np.nan_to_num(fld), levels=[15], colors='k', linewidths=1.4, transform=pc, zorder=3)
    gl = a.gridlines(draw_labels=True, x_inline=False, y_inline=False, lw=0.5, color='gray', alpha=0.6, ls='--')
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = {'size': FS - 2, 'weight': 'bold', 'rotation': 0}; gl.ylabel_style = {'size': FS - 2, 'weight': 'bold', 'rotation': 0}
    a.plot(KG6[1], KG6[0], 'o', ms=13, mfc='#ff7f0e', mec='k', mew=1.8, transform=pc, zorder=5)
    a.text(KG6[1] + 0.8, KG6[0] + 0.25, 'KG6', transform=pc, zorder=5, fontsize=FS + 1,
           bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='k', lw=0.8))
    for name, (lo, la) in {'Greenland': (-31.5, 69.6), 'Iceland': (-18.5, 64.6), 'Denmark Strait': (-28.5, 65.7)}.items():
        a.text(lo, la, name, transform=pc, ha='center', style='italic', fontsize=FS - 1, color='#3a3327', zorder=5)
    a.set_title(ttl, loc='left', fontsize=FS)
    cb = fig.colorbar(pm, ax=a, orientation='horizontal', shrink=0.85, pad=0.02); cb.set_label(cl)
fig.savefig(FIG + 'K0_KG6_osisaf_map.png', dpi=220); plt.close(fig)
print('wrote', FIG + 'K0_KG6_osisaf_map.png')
