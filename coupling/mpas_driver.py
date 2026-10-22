"""Hourly exchange loop: GOTM-identical bulk fluxes <-> MPAS-Ocean + ADC.

Usage: python3 mpas_driver.py <case_dir> [max_hours]
  <case_dir> holds init.nc, graph.info, namelist.ocean, streams.ocean and case.json
  (site name, start/stop time, position, forcing files relative to the case directory,
  optional run length 'hours'); icepack_in switches on the Icepack column (Block D),
  window_fluxes.nc switches to prescribed LES-window fluxes (Block W).
  Executables: $MPAS_OCEAN_MODEL (default models/mpas_adc/ocean_model) and, for Icepack,
  $ICEPACK_PIPE (default models/icepack_pipe/icepack_pipe).

Each 1-hour segment:
  1. SST = top-layer temperature of cell 1 from the latest restart (init.nc first).
  2. Fluxes = mean of gotm_airsea.bulk() at t+7.5, 22.5, 37.5, 52.5 min with
     that SST held fixed and ERA5 meteo linearly interpolated (as GOTM).
     Approximation vs GOTM: SST is updated once per hour, not every 450-1800 s.
  3. Write forcing_segment.nc, run ocean_model for 1 h (restart chain).
  4. Append cell-1 profiles (T, S, u, v, ADC moments) to hourly/<site>_adc.nc,
     log fluxes to hourly/<site>_fluxes.csv, delete the segment output and the
     non-midnight restart files the loop itself produced.
Resumable: restarts from Restart_timestamp if present.
"""
import os, re, sys, json, subprocess, datetime as dt, time
os.environ.setdefault('HDF5_USE_FILE_LOCKING', 'FALSE')   # NFS: allow readers while the driver appends
import numpy as np, netCDF4 as nc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gotm_airsea import Meteo, bulk

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.environ.get('MPAS_OCEAN_MODEL', os.path.join(HERE, '..', 'models', 'mpas_adc', 'ocean_model'))
TEMPLATE = os.path.join(HERE, 'forcing_template.nc')     # MPAS forcing-file layout, all fields zero

FMT = '%Y-%m-%d_%H:%M:%S'
RFMT = '%Y-%m-%d_%H.%M.%S'      # restart filename time format ($h.$m.$s)
SUB = [7.5, 22.5, 37.5, 52.5]          # minutes within the segment


def set_nml(path, key, value):
    s = open(path).read()
    s, n = re.subn(r'(?m)^(\s*%s\s*=\s*).*$' % key, r'\g<1>%s' % value, s)
    assert n == 1, key
    open(path, 'w').write(s)


def make_forcing_template(case):
    tmpl = TEMPLATE
    out = os.path.join(case, 'forcing_segment.nc')
    src = nc.Dataset(tmpl)
    dst = nc.Dataset(out, 'w', format=src.file_format)
    for n, d in src.dimensions.items():
        dst.createDimension(n, None if d.isunlimited() else len(d))
    for n, v in src.variables.items():
        w = dst.createVariable(n, v.dtype, v.dimensions)
        w[:] = v[:]
        if v.dtype.kind == 'f':
            w[:] = 0.0
    dst.close(); src.close()


def write_forcing(case, F, ice=None):
    """Open-water fluxes F (GOTM Kondo). With ice (Icepack hour means, grid-box means,
    positive into the ocean): open-water terms are weighted by (1 - aice) and the ice
    terms go into MPAS's sea-ice forcing fields (E3SM coupling convention)."""
    aow = 1.0 if ice is None else 1.0 - ice['aice_mean']
    d = nc.Dataset(os.path.join(case, 'forcing_segment.nc'), 'r+')
    d['windStressZonal'][:] = aow*F['taux']     # ice is motionless in 1D: no stress under ice
    d['windStressMeridional'][:] = aow*F['tauy']
    d['sensibleHeatFlux'][:] = aow*F['qh']
    d['latentHeatFlux'][:] = aow*F['qe']
    d['longWaveHeatFluxDown'][:] = aow*F['ql']  # GOTM net longwave (down - up)
    d['longWaveHeatFluxUp'][:] = 0.0
    d['shortWaveHeatFlux'][:] = aow*F['swr_net'] + (0.0 if ice is None else ice['fswthru'])
    d['evaporationFlux'][:] = 0.0               # GOTM: evap = 0 (flux_impact false)
    d['rainFlux'][:] = aow*F['precip']*1000.0   # m/s -> kg/m2/s (rain/snow on ice -> Icepack)
    if ice is not None:
        d['seaIceHeatFlux'][:] = ice['fhocn'] + ice['frazil_heat']
        d['seaIceFreshWaterFlux'][:] = ice['fresh']
        d['seaIceSalinityFlux'][:] = ice['fsalt']
    d.close()


class WindowFluxes:
    """Block W: prescribed hourly fluxes of an LES window (window_fluxes.nc, copy of
    model_dev_les/les_kg6_<B>.nc), linearly interpolated as in the LES and GOTM.
    Q non-solar and I0 net shortwave into the ocean; tau kinematic (x rho0 = 1027 -> Pa);
    upward salt flux ssf -> 'rain' P = ssf / SA_top0 (as the GOTM window runs), passed as
    rainFlux = P * rho_sw so that MPAS's virtual salt flux S_top*rainFlux/rho_sw = S_top*P."""
    RHO0, RHO_SW = 1027.0, 1026.0

    def __init__(self, path):
        d = nc.Dataset(path)
        self.t0 = dt.datetime.fromisoformat(d.window[:19])
        self.ts = np.array(d['time'][:])
        self.v = {k: np.array(d[k][:]) for k in ('Q', 'I0', 'taux', 'tauy', 'ssf')}
        self.sa0 = float(d['S'][-1])
        d.close()

    def __call__(self, t):
        acc = dict(qh=0.0, qe=0.0, ql=0.0, swr_net=0.0, heat=0.0, taux=0.0, tauy=0.0, precip=0.0, albedo=0.0)
        for m in SUB:
            s = (t + dt.timedelta(minutes=m) - self.t0).total_seconds()
            f = {k: float(np.interp(s, self.ts, v)) for k, v in self.v.items()}
            acc['qh'] += f['Q']/len(SUB); acc['heat'] += f['Q']/len(SUB); acc['swr_net'] += f['I0']/len(SUB)
            acc['taux'] += f['taux']*self.RHO0/len(SUB); acc['tauy'] += f['tauy']*self.RHO0/len(SUB)
            acc['precip'] += f['ssf']/self.sa0*self.RHO_SW/1000.0/len(SUB)   # driver writes precip*1000
        return acc


def segment_fluxes(t, sst, met, lon, lat):
    acc = None
    for m in SUB:
        tt = t + dt.timedelta(minutes=m)
        r = bulk(tt, sst, met.at(tt), lon, lat)
        acc = {k: r[k]/len(SUB) for k in r} if acc is None else {k: acc[k] + r[k]/len(SUB) for k in r}
    return acc


class Hourly:
    """Compact cell-1 output file."""
    PROF = ['temperature', 'salinity']
    ADC = ['w2', 'u2', 'v2', 'wt', 'ws', 'w3', 'uw', 'vw', 'areaFraction', 'Mc', 'eps',
           'length', 'KspsU', 'KspsD']

    def __init__(self, path, init_path):
        self.path = path
        if not os.path.exists(path):
            I = nc.Dataset(init_path)
            nz = len(I.dimensions['nVertLevels'])
            d = nc.Dataset(path, 'w')
            d.createDimension('time', None); d.createDimension('nz', nz); d.createDimension('nzp1', nz + 1)
            d.createVariable('time', 'f8', ('time',)).units = 'hours since 2018-01-01 00:00:00'
            d.createVariable('zmid', 'f8', ('nz',))[:] = I['refZMid'][:]
            zt = np.concatenate([[0.0], -np.array(I['refBottomDepth'][:])])
            d.createVariable('ztop', 'f8', ('nzp1',))[:] = zt
            for v in self.PROF + ['u', 'v']:
                d.createVariable(v, 'f8', ('time', 'nz'))
            d.createVariable('boundaryLayerDepth', 'f8', ('time',))
            d.close(); I.close()
            # least-squares edge -> (u, v) for cell 1
        I = nc.Dataset(init_path)
        e = np.array(I['edgesOnCell'][0, :I['nEdgesOnCell'][0]]) - 1
        a = np.array(I['angleEdge'][:])[e]
        self.edges = e
        self.A = np.linalg.pinv(np.vstack([np.cos(a), np.sin(a)]).T)
        I.close()

    def append(self, seg, rec, t):
        S = nc.Dataset(seg)
        d = nc.Dataset(self.path, 'r+')
        n = len(d['time'])
        d['time'][n] = (t - dt.datetime(2018, 1, 1)).total_seconds()/3600.0
        for v in self.PROF:
            d[v][n, :] = S[v][rec, 0, :]
        nv = np.array(S['normalVelocity'][rec])[self.edges, :]
        uv = self.A @ nv
        d['u'][n, :] = uv[0]; d['v'][n, :] = uv[1]
        for v in self.ADC:
            if v not in d.variables:     # vertical dim taken from MPAS (nVertLevels vs nVertLevelsP1)
                d.createVariable(v, 'f8', ('time', 'nzp1' if S[v].dimensions[-1] == 'nVertLevelsP1' else 'nz'))
            d[v][n, :] = S[v][rec, 0, :]
        d['boundaryLayerDepth'][n] = S['boundaryLayerDepth'][rec, 0]
        d.close(); S.close()


def latest_state(case, start):
    ts = os.path.join(case, 'Restart_timestamp')
    if os.path.exists(ts):
        t = dt.datetime.strptime(open(ts).read().strip(), FMT)
        return t, os.path.join(case, 'restarts', 'restart.%s.nc' % t.strftime(RFMT))
    return start, os.path.join(case, 'init.nc')


def main():
    case = os.path.abspath(sys.argv[1])
    S = json.load(open(os.path.join(case, 'case.json')))
    site = S['site']
    max_hours = int(sys.argv[2]) if len(sys.argv) > 2 else S.get('hours', 10**9)
    for d in ('restarts', 'output', 'analysis_members', 'hourly'):
        os.makedirs(os.path.join(case, d), exist_ok=True)
    start = dt.datetime.strptime(S['start'], FMT); stop = dt.datetime.strptime(S['stop'], FMT)
    met = Meteo(os.path.join(case, S['meteo']))
    wf = os.path.join(case, 'window_fluxes.nc')
    window = WindowFluxes(wf) if os.path.exists(wf) else None
    nml = os.path.join(case, 'namelist.ocean')
    if not os.path.exists(os.path.join(case, 'forcing_segment.nc')):
        make_forcing_template(case)
    H = Hourly(os.path.join(case, 'hourly', '%s_adc.nc' % site), os.path.join(case, 'init.nc'))
    flog = os.path.join(case, 'hourly', '%s_fluxes.csv' % site)
    if not os.path.exists(flog):
        open(flog, 'w').write('time,sst,qh,qe,ql,swr_net,heat_nonsolar,taux,tauy,precip,albedo,wall_s\n')
    env = dict(os.environ)
    done = 0
    pipe, ice_log, aice_prev = None, None, 0.0
    if os.path.exists(os.path.join(case, 'icepack_in')):      # Block D: Icepack column
        from icepack_pipe import IcepackPipe, snow_series
        t_now, _ = latest_state(case, start)
        if t_now != start:
            sys.exit('Icepack runs keep the ice state in memory: restart from the beginning')
        os.makedirs(os.path.join(case, 'restart'), exist_ok=True)
        pipe = IcepackPipe(case, S['lat'], S['lon'], 112.5, start, 32)
        snow = snow_series(os.path.join(case, S['snow']))
        ice_log = open(os.path.join(case, 'hourly', '%s_icepack.csv' % site), 'w')
        ice_log.write('time,aice,vice,vsno,Tsfc,Sice,fhocn_ai,fswthru_ai,fresh_ai,fsalt_ai,frazil_heat\n')
    while done < max_hours:
        t, state = latest_state(case, start)
        if t >= stop:
            print('reached stop', t); break
        with nc.Dataset(state) as D:
            sst = float(D['temperature'][0, 0, 0])
        if not np.isfinite(sst):
            sys.exit('non-finite SST at %s' % t)
        F = window(t) if window is not None else segment_fluxes(t, sst, met, S['lon'], S['lat'])
        ice = None
        if pipe is not None:
            with nc.Dataset(state) as D:
                sss = float(D['salinity'][0, 0, 0]); htop = float(D['layerThickness'][0, 0, 0])
            # ocean held for the hour: spread the top-layer freezing/melting potential over
            # the nsub Icepack sub-steps (h_top/nsub) so it is used once, not nsub times
            ice = pipe.hour(t, met, snow, sst, sss, htop/pipe.nsub)
            ice['aice_mean'] = 0.5*(aice_prev + ice['aice']); aice_prev = ice['aice']
            ice_log.write('%s,%.6f,%.6e,%.6e,%.4f,%.4f,%.6e,%.6e,%.6e,%.6e,%.6e\n' % (
                (t + dt.timedelta(hours=1)).strftime(FMT), ice['aice'], ice['vice'], ice['vsno'], ice['Tsfc'],
                ice['Sice'], ice['fhocn'], ice['fswthru'], ice['fresh'], ice['fsalt'], ice['frazil_heat']))
            ice_log.flush()
        write_forcing(case, F, ice)
        first = (t == start)
        set_nml(nml, 'config_do_restart', '.false.' if first else '.true.')
        set_nml(nml, 'config_start_time', "'%s'" % S['start'] if first else "'file'")
        err = os.path.join(case, 'log.ocean.0000.err')
        if os.path.exists(err):
            os.remove(err)              # stale error log from a previous failed attempt
        w0 = time.time()
        r = subprocess.run(['mpirun', '-n', '1', os.path.abspath(EXE), '-n', 'namelist.ocean', '-s', 'streams.ocean'],
                           cwd=case, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        wall = time.time() - w0
        err = os.path.join(case, 'log.ocean.0000.err')
        if r.returncode != 0 or os.path.exists(err):
            msg = open(err).read() if os.path.exists(err) else r.stdout.decode()[-2000:]
            os.rename(os.path.join(case, 'log.ocean.0000.out'),
                      os.path.join(case, 'log_failed_%s.out' % t.strftime(FMT)))
            sys.exit('segment %s failed (rc=%s):\n%s' % (t, r.returncode, msg[-2000:]))
        t1 = t + dt.timedelta(hours=1)
        seg = os.path.join(case, 'output', 'segment.nc')   # one file, truncated by every segment
        with nc.Dataset(seg) as Sg:
            times = [b''.join(x).decode().strip() for x in Sg['xtime'][:]]
        for attempt in range(10):          # tolerate a concurrent reader holding the file
            try:
                if first and attempt == 0:
                    H.append(seg, times.index(t.strftime(FMT)), t)
                H.append(seg, times.index(t1.strftime(FMT)), t1)
                break
            except OSError as e:
                print('append retry', attempt, e, flush=True); time.sleep(10)
        else:
            sys.exit('could not append %s' % t1)
        with open(flog, 'a') as fh:
            fh.write('%s,%.5f,%.4f,%.4f,%.4f,%.4f,%.4f,%.6f,%.6f,%.4e,%.4f,%.2f\n' % (
                t.strftime(FMT), sst, F['qh'], F['qe'], F['ql'], F['swr_net'], F['heat'],
                F['taux'], F['tauy'], F['precip'], F['albedo'], wall))
        if not first and t.hour != 0:
            old = os.path.join(case, 'restarts', 'restart.%s.nc' % t.strftime(RFMT))
            if os.path.exists(old):
                os.remove(old)
        done += 1
    if pipe is not None:
        pipe.close(); ice_log.close()
    print('done', done, 'segments; last', latest_state(case, start)[0])


if __name__ == '__main__':
    main()
