"""GOTM v6 air-sea fluxes (Kondo 1975 bulk formulae), ported line-for-line.

Source: model/gotm/src/airsea/{kondo,humidity,airsea,albedo_water,
solar_zenith_angle}.F90 (GOTM v6.0.0-211-g84f8d8e9), configuration as in
case_LB9 / case_KG6 gotm.yaml:
  fluxes: kondo; hum type: dew_point; swr type 2 (downward); longwave
  type 2 (downward, from file); albedo: payne; flux_impact (rain_impact):
  false  -> Kondo returns evap = 0 and no rain terms; GUSTINESS only affects
  fairall (not used).
Sign convention (GOTM): fluxes positive INTO the ocean.
Returns stresses in Pa (GOTM stores tx = tau/rho0, rho0 = 1027 in gotm.yaml).
"""
import math
import datetime as _dt
import numpy as np

KELVIN = 273.15
CONST06 = 0.62198
RGAS = 287.1
CPA = 1008.0
EMISS = 0.97
BOLZ = 5.670374419e-8

_A = (6.107799961, 4.436518521e-1, 1.428945805e-2, 2.650648471e-4,
      3.031240396e-6, 2.034080948e-8, 6.136820929e-11)


def _esat_hpa(t):
    a1, a2, a3, a4, a5, a6, a7 = _A
    return a1 + t*(a2 + t*(a3 + t*(a4 + t*(a5 + t*(a6 + t*a7)))))


def humidity_dewpoint(dew, airp, tw, ta):
    """humidity.F90, hum_method=3 (dew point). Returns qs, qa, ea, rhoa."""
    es = 0.98 * _esat_hpa(tw) * 100.0
    qs = CONST06*es/(airp - 0.377*es)
    if dew >= 100.0:
        dew = dew - KELVIN
    ea = _esat_hpa(dew) * 100.0
    qa = CONST06*ea/(airp - 0.377*ea)
    rhoa = airp/(RGAS*(ta + KELVIN)*(1.0 + CONST06*qa))
    return qs, qa, ea, rhoa


_AE_D = (0., 0.771, 0.867, 1.2, 0.)
_AE_H = (0., 0.927, 1.15, 1.17, 1.652)
_AE_E = (0., 0.969, 1.18, 1.196, 1.68)
_BE_D = (1.08, 0.0858, 0.0667, 0.025, 0.073)
_BE_H = (1.185, 0.0546, 0.01, 0.0075, -0.017)
_BE_E = (1.23, 0.0521, 0.01, 0.008, -0.016)
_CE_H = (0., 0., 0., -0.00045, 0.)
_CE_E = (0., 0., 0., -0.0004, 0.)
_PE_D = (-0.15, 1., 1., 1., 1.)
_PE_H = (-0.157, 1., 1., 1., 1.)
_PE_E = (-0.16, 1., 1., 1., 1.)


def kondo(sst, airt, u10, v10, qs, qa, rhoa):
    """kondo.F90 (non-OLD_KONDO branch, rain_impact=false).
    Returns taux, tauy [Pa], qe, qh [W/m2], evap (=0)."""
    eps = 1.0e-12
    w = math.sqrt(u10*u10 + v10*v10)
    L = (2.5 - 0.00234*sst)*1.e6
    s0 = 0.25*(sst - airt)/(w + 1.0e-10)**2
    s = s0*abs(s0)/(abs(s0) + 0.01)
    if w < 2.2:
        x = math.log(w + eps)
        cdd = (_BE_D[0]*math.exp(_PE_D[0]*x))*1.0e-3
        chd = (_BE_H[0]*math.exp(_PE_H[0]*x))*1.0e-3
        ced = (_BE_E[0]*math.exp(_PE_E[0]*x))*1.0e-3
    else:
        i = 1 if w < 5.0 else 2 if w < 8.0 else 3 if w < 25.0 else 4
        x = math.exp(math.log(w + eps))
        cdd = (_AE_D[i] + _BE_D[i]*x)*1.0e-3
        chd = (_AE_H[i] + _BE_H[i]*x + _CE_H[i]*(w - 8.0)**2)*1.0e-3
        ced = (_AE_E[i] + _BE_E[i]*x + _CE_E[i]*(w - 8.0)**2)*1.0e-3
    if s < 0.0:
        x = 0.1 + 0.03*s + 0.9*math.exp(4.8*s) if s > -3.3 else 0.0
        cdd, chd, ced = x*cdd, x*chd, x*ced
    else:
        cdd = cdd*(1.0 + 0.47*math.sqrt(s))
        chd = chd*(1.0 + 0.63*math.sqrt(s))
        ced = ced*(1.0 + 0.63*math.sqrt(s))
    qh = -chd*CPA*rhoa*w*(sst - airt)
    qe = -ced*L*rhoa*w*(qs - qa)
    tmp = cdd*rhoa*w
    return tmp*u10, tmp*v10, qe, qh, 0.0


def solar_zenith_angle(yday, hh, dlon, dlat):
    """solar_zenith_angle.F90; yday = integer day of year, hh = UTC hours."""
    pi = 3.14159265358979323846
    d2r = pi/180.
    th0 = 2.*pi*yday/365.25
    sundec = (0.006918 - 0.399912*math.cos(th0) + 0.070257*math.sin(th0)
              - 0.006758*math.cos(2*th0) + 0.000907*math.sin(2*th0)
              - 0.002697*math.cos(3*th0) + 0.001480*math.sin(3*th0))
    thsun = (hh - 12.)*15.*d2r + d2r*dlon
    coszen = (math.sin(d2r*dlat)*math.sin(sundec)
              + math.cos(d2r*dlat)*math.cos(sundec)*math.cos(thsun))
    coszen = max(coszen, 0.0)
    return math.acos(coszen)/d2r


_ALB1 = (.719, .656, .603, .480, .385, .300, .250, .193, .164, .131,
         .103, .084, .071, .061, .054, .043, .039, .036, .034, .034)
_ZA = (90., 88., 86., 84., 82., 80., 78., 76., 74., 70.,
       66., 62., 58., 54., 50., 40., 30., 20., 10., 0.0)
_DZA = (2.0,)*8 + (4.0,)*6 + (10.0,)*5


def albedo_payne(zen):
    """albedo_water.F90, Payne (1972). Fortran int() truncation of jab."""
    if zen >= 74.:
        jab = int(.5*(90. - zen) + 1.)
    elif zen >= 50.:
        jab = int(.23*(74. - zen) + 9.)
    else:
        jab = int(.10*(50. - zen) + 15.)
    if jab == 20:
        return _ALB1[19]
    dzen = (_ZA[jab-1] - zen)/_DZA[jab-1]
    return _ALB1[jab-1] + dzen*(_ALB1[jab] - _ALB1[jab-1])


def yearday(t):
    """GOTM 'yearday' = day of year (1-based integer)."""
    return t.timetuple().tm_yday


def bulk(t, sst, met, dlon, dlat):
    """All surface fluxes at time t (datetime, UTC) for surface temperature
    sst [C], meteo dict met with u10,v10 [m/s], airp [Pa], airt [C],
    dew [C], cloud, precip [m/s], swr [W/m2 downward], lwd [W/m2 downward].
    Returns dict, all positive into the ocean."""
    qs, qa, ea, rhoa = humidity_dewpoint(met['dew'], met['airp'], sst, met['airt'])
    ql = met['lwd'] - BOLZ*EMISS*(sst + KELVIN)**4
    taux, tauy, qe, qh, evap = kondo(sst, met['airt'], met['u10'], met['v10'], qs, qa, rhoa)
    hh = t.hour + t.minute/60. + t.second/3600.
    zen = solar_zenith_angle(yearday(t), hh, dlon, dlat)
    alb = albedo_payne(zen)
    return dict(qh=qh, qe=qe, ql=ql, heat=ql + qe + qh, taux=taux, tauy=tauy,
                swr_down=met['swr'], albedo=alb, swr_net=met['swr']*(1.0 - alb),
                precip=met['precip'], evap=evap, rhoa=rhoa)


def read_meteo(path):
    """meteo.dat columns: date time u10 v10 airp[hPa] airt dew cloud precip swr lwd"""
    times, rows = [], []
    for line in open(path):
        p = line.split()
        if len(p) < 11:
            continue
        times.append(_dt.datetime.strptime(p[0] + ' ' + p[1], '%Y/%m/%d %H:%M:%S'))
        rows.append([float(x) for x in p[2:11]])
    return times, np.array(rows)


class Meteo:
    """Linear-in-time interpolation of meteo.dat (as GOTM's input manager)."""
    keys = ('u10', 'v10', 'airp', 'airt', 'dew', 'cloud', 'precip', 'swr', 'lwd')

    def __init__(self, path):
        self.times, self.data = read_meteo(path)
        t0 = self.times[0]
        self.sec = np.array([(t - t0).total_seconds() for t in self.times])
        self.t0 = t0

    def at(self, t):
        s = (t - self.t0).total_seconds()
        row = np.array([np.interp(s, self.sec, self.data[:, j]) for j in range(self.data.shape[1])])
        m = dict(zip(self.keys, row))
        m['airp'] *= 100.0   # hPa -> Pa (airsea.F90 scale_factor 100)
        return m
