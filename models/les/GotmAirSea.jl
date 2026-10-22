# GOTM v6 air-sea fluxes (Kondo 1975) for the LES - a line-for-line copy of
# model_dev_ADC_mushy/mpas_coupled/realforcing/gotm_airsea.py, which reproduces GOTM's own
# output fluxes exactly (validate_fluxes.py, LB9 and KG6). Configuration as the KG6 GOTM
# runs: dew-point humidity, downward LW and SW from meteo.dat, Payne albedo, evap = 0.
# Positive = into the ocean; stresses in Pa.
module GotmAirSea

using Dates
export Meteo, read_meteo, meteo_at, bulk

const KELVIN = 273.15
const CONST06 = 0.62198
const RGAS = 287.1
const CPA = 1008.0
const EMISS = 0.97
const BOLZ = 5.670374419e-8
const A = (6.107799961, 4.436518521e-1, 1.428945805e-2, 2.650648471e-4,
           3.031240396e-6, 2.034080948e-8, 6.136820929e-11)

esat_hpa(t) = A[1] + t*(A[2] + t*(A[3] + t*(A[4] + t*(A[5] + t*(A[6] + t*A[7])))))

function humidity_dewpoint(dew, airp, tw, ta)
    es = 0.98 * esat_hpa(tw) * 100.0
    qs = CONST06*es/(airp - 0.377*es)
    dew >= 100.0 && (dew -= KELVIN)
    ea = esat_hpa(dew) * 100.0
    qa = CONST06*ea/(airp - 0.377*ea)
    rhoa = airp/(RGAS*(ta + KELVIN)*(1.0 + CONST06*qa))
    return qs, qa, ea, rhoa
end

const AE_D = (0., 0.771, 0.867, 1.2, 0.);   const AE_H = (0., 0.927, 1.15, 1.17, 1.652)
const AE_E = (0., 0.969, 1.18, 1.196, 1.68); const BE_D = (1.08, 0.0858, 0.0667, 0.025, 0.073)
const BE_H = (1.185, 0.0546, 0.01, 0.0075, -0.017); const BE_E = (1.23, 0.0521, 0.01, 0.008, -0.016)
const CE_H = (0., 0., 0., -0.00045, 0.);    const CE_E = (0., 0., 0., -0.0004, 0.)
const PE_D = (-0.15, 1., 1., 1., 1.);       const PE_H = (-0.157, 1., 1., 1., 1.)
const PE_E = (-0.16, 1., 1., 1., 1.)

function kondo(sst, airt, u10, v10, qs, qa, rhoa)
    eps = 1.0e-12
    w = sqrt(u10^2 + v10^2)
    L = (2.5 - 0.00234*sst)*1.e6
    s0 = 0.25*(sst - airt)/(w + 1.0e-10)^2
    s = s0*abs(s0)/(abs(s0) + 0.01)
    if w < 2.2
        x = log(w + eps)
        cdd = BE_D[1]*exp(PE_D[1]*x)*1e-3; chd = BE_H[1]*exp(PE_H[1]*x)*1e-3; ced = BE_E[1]*exp(PE_E[1]*x)*1e-3
    else
        i = w < 5.0 ? 2 : w < 8.0 ? 3 : w < 25.0 ? 4 : 5
        x = exp(log(w + eps))
        cdd = (AE_D[i] + BE_D[i]*x)*1e-3
        chd = (AE_H[i] + BE_H[i]*x + CE_H[i]*(w - 8.0)^2)*1e-3
        ced = (AE_E[i] + BE_E[i]*x + CE_E[i]*(w - 8.0)^2)*1e-3
    end
    if s < 0.0
        x = s > -3.3 ? 0.1 + 0.03*s + 0.9*exp(4.8*s) : 0.0
        cdd *= x; chd *= x; ced *= x
    else
        cdd *= (1.0 + 0.47*sqrt(s)); chd *= (1.0 + 0.63*sqrt(s)); ced *= (1.0 + 0.63*sqrt(s))
    end
    qh = -chd*CPA*rhoa*w*(sst - airt)
    qe = -ced*L*rhoa*w*(qs - qa)
    tmp = cdd*rhoa*w
    return tmp*u10, tmp*v10, qe, qh
end

function solar_zenith_angle(yday, hh, dlon, dlat)
    pi_ = 3.14159265358979323846; d2r = pi_/180.
    th0 = 2*pi_*yday/365.25
    sundec = 0.006918 - 0.399912*cos(th0) + 0.070257*sin(th0) - 0.006758*cos(2th0) +
             0.000907*sin(2th0) - 0.002697*cos(3th0) + 0.001480*sin(3th0)
    thsun = (hh - 12.0)*15.0*d2r + d2r*dlon
    coszen = max(sin(d2r*dlat)*sin(sundec) + cos(d2r*dlat)*cos(sundec)*cos(thsun), 0.0)
    return acos(coszen)/d2r
end

const ALB1 = (.719, .656, .603, .480, .385, .300, .250, .193, .164, .131,
              .103, .084, .071, .061, .054, .043, .039, .036, .034, .034)
const ZA = (90., 88., 86., 84., 82., 80., 78., 76., 74., 70., 66., 62., 58., 54., 50., 40., 30., 20., 10., 0.0)
const DZA = (2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 10.0, 10.0, 10.0, 10.0, 10.0)

function albedo_payne(zen)
    jab = zen >= 74. ? trunc(Int, .5*(90. - zen) + 1.) :
          zen >= 50. ? trunc(Int, .23*(74. - zen) + 9.) : trunc(Int, .10*(50. - zen) + 15.)
    jab == 20 && return ALB1[20]
    dzen = (ZA[jab] - zen)/DZA[jab]
    return ALB1[jab] + dzen*(ALB1[jab+1] - ALB1[jab])
end

struct Meteo
    t0::DateTime
    sec::Vector{Float64}
    data::Matrix{Float64}     # u10 v10 airp[hPa] airt dew cloud precip swr lwd
end

function read_meteo(path)
    ts = DateTime[]; rows = Vector{Vector{Float64}}()
    for line in eachline(path)
        p = split(line)
        length(p) < 11 && continue
        push!(ts, DateTime(p[1] * " " * p[2], dateformat"yyyy/mm/dd HH:MM:SS"))
        push!(rows, parse.(Float64, p[3:11]))
    end
    sec = [Dates.value(t - ts[1]) / 1000 for t in ts]
    return Meteo(ts[1], sec, permutedims(reduce(hcat, rows)))
end

function meteo_at(m::Meteo, t::DateTime)
    s = Dates.value(t - m.t0) / 1000
    k = clamp(searchsortedlast(m.sec, s), 1, length(m.sec) - 1)
    w = clamp((s - m.sec[k]) / (m.sec[k+1] - m.sec[k]), 0.0, 1.0)
    r = (1 - w) .* m.data[k, :] .+ w .* m.data[k+1, :]
    return (u10=r[1], v10=r[2], airp=r[3]*100.0, airt=r[4], dew=r[5], cloud=r[6], precip=r[7], swr=r[8], lwd=r[9])
end

"All fluxes at time t for surface (in-situ) temperature sst [C]; positive into the ocean."
function bulk(t::DateTime, sst, met, dlon, dlat)
    qs, qa, ea, rhoa = humidity_dewpoint(met.dew, met.airp, sst, met.airt)
    ql = met.lwd - BOLZ*EMISS*(sst + KELVIN)^4
    taux, tauy, qe, qh = kondo(sst, met.airt, met.u10, met.v10, qs, qa, rhoa)
    hh = hour(t) + minute(t)/60 + second(t)/3600
    zen = solar_zenith_angle(dayofyear(t), hh, dlon, dlat)
    alb = albedo_payne(zen)
    return (qh=qh, qe=qe, ql=ql, heat=ql + qe + qh, taux=taux, tauy=tauy,
            swr_net=met.swr*(1 - alb), albedo=alb, precip=met.precip)
end

end
