import numpy as num
from pyrocko import orthodrome as od


def get_mean_x(xs):
    return num.mean(xs, axis=0)


def get_mean_x_and_gm(problem, xs, misfits):
    gms = problem.combine_misfits(misfits)
    return num.mean(xs, axis=0), num.mean(gms)


def get_best_x(problem, xs, misfits):
    gms = problem.combine_misfits(misfits)
    ibest = num.argmin(gms)
    return xs[ibest, :]


def get_best_x_and_gm(problem, xs, misfits):
    gms = problem.combine_misfits(misfits)
    ibest = num.argmin(gms)
    return xs[ibest, :], gms[ibest]


def get_mean_source(problem, xs):
    x_mean = get_mean_x(xs)
    source = problem.get_source(x_mean)
    return source


def get_best_source(problem, xs, misfits):
    x_best = get_best_x(problem, xs, misfits)
    source = problem.get_source(x_best)
    return source


def mean_latlondist(lats, lons):
    if len(lats) == 0:
        return 0., 0., 1000.
    else:
        ns, es = od.latlon_to_ne_numpy(lats[0], lons[0], lats, lons)
        n, e = num.mean(ns), num.mean(es)
        dists = num.sqrt((ns-n)**2 + (es-e)**2)
        lat, lon = od.ne_to_latlon(lats[0], lons[0], n, e)
        return float(lat), float(lon), float(num.max(dists))


def stations_mean_latlondist(stations):
    lats = num.array([s.lat for s in stations])
    lons = num.array([s.lon for s in stations])
    return mean_latlondist(lats, lons)
