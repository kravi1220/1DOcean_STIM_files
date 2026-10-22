"""LES grid convergence for the KG6 fixed-flux windows (1, 2, 4 m; Section 3.5, Appendix A).
Usage: python3 les_convergence.py > ../results/tables/les_convergence.txt"""
import netCDF4 as nc, numpy as np, gsw
import paths as P
def load(f):
    d=nc.Dataset(f); t=np.array(d['time'][:])/3600; z=np.array(d['z_aac'][:]); T=np.array(d['T'][:]); S=np.array(d['S'][:])
    zf=np.array(d['z_aaf'][:]); w=np.array(d['w'][:]); ww=np.array(d['ww'][:])-w**2
    wT=np.array(d['wT'][:])-w*np.array([np.interp(zf,z,T[k]) for k in range(len(t))])
    dz=np.array(d['Δz_aac'][:]); return dict(t=t,z=z,T=T,S=S,zf=zf,ww=ww,wT=wT,dz=dz)
def mld(M,k):
    sig=gsw.sigma0(M['S'][k],M['T'][k]); zz=np.arange(0,160,0.1); s=np.interp(-zz,M['z'],sig); s10=np.interp(-10,M['z'],sig)
    j=np.where((zz>10)&(s>s10+0.03))[0]; return zz[j[0]] if j.size else np.nan
def row(M,h,ref=None):
    k=int(np.argmin(abs(M['t']-h))); hc=1027*3991.87*np.sum((M['T'][k]-M['T'][0])*M['dz'])/1e6
    zg=np.arange(-80,-1,0.5); Ti=np.interp(zg,M['z'],M['T'][k])
    rms=np.sqrt(np.mean((Ti-np.interp(zg,ref['z'],ref['T'][int(np.argmin(abs(ref['t']-h)))]))**2)) if ref else 0.0
    wtm=np.interp(-np.arange(5,70,0.5),M['zf'],M['wT'][k]); wwmax=M['ww'][k][M['zf']>-80].max()
    return '%5.0fh  SST %7.4f  SSS %7.4f  MLD %5.1f m  dHC %8.2f MJ/m2  max w\'w\' %.2e  T-RMS vs finest (0-80 m) %.4f C'%(h,M['T'][k][-1],M['S'][k][-1],mld(M,k),hc,wwmax,rms)
for b,res,hours in [('B1',['1.0','2.0','4.0'],[6,12,18,24]),('B2',['2.0','4.0'],[24,48,72,96,120])]:
    Ms={r:load(P.LES+'les_kg6_%s_dx%s_profiles.nc'%(b,r)) for r in res}; ref=Ms[res[0]]
    print('== %s (finest = %s m)'%(b,res[0]))
    for h in hours:
        for r in res: print('  %s m '%r+row(Ms[r],h,ref))
