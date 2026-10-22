"""Python side of the persistent Icepack column (models/icepack_pipe).

IcepackPipe(run_dir, lat, lon, dt, start, nsub) starts the program in run_dir (which must
contain icepack_in); hour(t, met, snow, sst_insitu, sss_psu, h_top) sends nsub sub-steps
of forcing for the hour starting at t (atmosphere evaluated at the END of each sub-step,
as GOTM does) and returns a dict with the ice state and hour-mean ice->ocean fluxes.
Atmospheric quantities are computed with the GOTM-validated gotm_airsea functions.
"""
import os, subprocess, datetime as dt
from gotm_airsea import humidity_dewpoint

EXE = os.environ.get('ICEPACK_PIPE', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'models', 'icepack_pipe', 'icepack_pipe'))
KEYS = ['aice', 'vice', 'vsno', 'Tsfc', 'Sice', 'fhocn', 'fswthru', 'fresh', 'fsalt', 'frazil_heat']
RHO0, CP = 1027.0, 3991.86795711963       # GOTM rho0 and TEOS-10 cp0, as in the GOTM coupling


class IcepackPipe:
    def __init__(self, run_dir, lat, lon, dtsub, start, nsub):
        self.dt, self.nsub = dtsub, nsub
        y0 = dt.datetime(start.year, 1, 1)
        start_sec = (start - y0).total_seconds()
        self.p = subprocess.Popen([EXE], cwd=run_dir, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT, text=True, bufsize=1)
        self.p.stdin.write('%.10f %.10f %.6f %.3f %d\n' % (lat, lon, dtsub, start_sec, nsub))
        self.p.stdin.flush()
        log = []
        for line in self.p.stdout:
            log.append(line)
            if line.strip() == 'READY':
                break
            if 'STOP' in line or 'ABORT' in line:
                raise RuntimeError('icepack_pipe failed to start:\n' + ''.join(log))
        else:
            raise RuntimeError('icepack_pipe ended during start-up:\n' + ''.join(log))

    def hour(self, t, met, snow, sst_insitu, sss_psu, h_top):
        """met: gotm_airsea.Meteo; snow: callable(datetime) -> snowfall [m/s w.e.]."""
        lines = []
        for j in range(1, self.nsub + 1):
            tj = t + dt.timedelta(seconds=j * self.dt)
            m = met.at(tj)
            _, qa, _, rhoa = humidity_dewpoint(m['dew'], m['airp'], sst_insitu, m['airt'])
            lines.append('%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.6f %.6f\n' % (
                m['airt'], qa, m['u10'], m['v10'], rhoa, m['lwd'], m['swr'], m['precip'], snow(tj),
                sst_insitu, sss_psu, h_top, RHO0, CP))
        self.p.stdin.write(''.join(lines)); self.p.stdin.flush()
        info = []
        while True:                      # Icepack may print messages (e.g. restart dumps) first
            out = self.p.stdout.readline()
            if out == '':
                raise RuntimeError('icepack_pipe ended at %s:\n%s' % (t, ''.join(info[-40:])))
            try:
                v = [float(x) for x in out.split()]
                if len(v) == len(KEYS):
                    break
            except ValueError:
                pass
            info.append(out)
            if 'ABORT' in out or 'STOP' in out:
                raise RuntimeError('icepack_pipe aborted at %s:\n%s' % (t, ''.join(info[-40:]) + ''.join(self.p.stdout.readlines()[:40])))
        self.messages = info
        return dict(zip(KEYS, v))

    def close(self):
        try:
            self.p.stdin.write('END\n'); self.p.stdin.flush()
        except Exception:
            pass
        self.p.wait(timeout=60)


def snow_series(path):
    import numpy as np
    T, V = [], []
    for l in open(path):
        p = l.split()
        T.append(dt.datetime.strptime(p[0] + ' ' + p[1], '%Y/%m/%d %H:%M:%S')); V.append(float(p[2]))
    s = np.array([(x - T[0]).total_seconds() for x in T]); V = np.array(V)
    return lambda t: float(np.interp((t - T[0]).total_seconds(), s, V))
