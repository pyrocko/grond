import math
import numpy as num
from pyrocko import moment_tensor, orthodrome
from grond.meta import GrondError


def get_distance_mt_l2(eventi, eventj):
    '''
    L2 norm among two moment tensors, with 6 independet entries
    '''

    # ni = (eventi.mxx)**2 + (eventi.myy)**2 + (eventi.mzz)**2 + \
    #      (eventi.mxy)**2 + (eventi.mxz)**2 + (eventi.myz)**2
    # nj = (eventj.mxx)**2 + (eventj.myy)**2 + (eventj.mzz)**2 + \
    #      (eventj.mxy)**2 + (eventj.mxz)**2 + (eventj.myz)**2

    # mixx, miyy, mizz = eventi.mxx / ni, eventi.myy / ni, eventi.mzz / ni
    # mixy, mixz, miyz = eventi.mxy / ni, eventi.mxz / ni, eventi.myz / ni
    # mjxx, mjyy, mjzz = eventj.mxx / nj, eventj.myy / nj, eventj.mzz / nj
    # mjxy, mjxz, mjyz = eventj.mxy / nj, eventj.mxz / nj, eventj.myz / nj

    d = (eventi.mxx - eventj.mxx)**2 + (eventi.myy - eventj.myy)**2 + \
        (eventi.mzz - eventj.mzz)**2 + (eventi.mxy - eventj.mxy)**2 + \
        (eventi.mxz - eventj.mxz)**2 + (eventi.myz - eventj.myz)**2

    d = 0.5 * math.sqrt(d)

    return d


def get_distance_mt_l1(eventi, eventj):
    '''
    L1 norm among two moment tensors, with 6 independet entries
    '''

    # ni = abs(eventi.mxx) + abs(eventi.myy) + abs(eventi.mzz) + \
    #      abs(eventi.mxy) + abs(eventi.mxz) + abs(eventi.myz)
    # nj = abs(eventj.mxx) + abs(eventj.myy) + abs(eventj.mzz) + \
    #      abs(eventj.mxy) + abs(eventj.mxz) + abs(eventj.myz)
    #
    # mixx, miyy, mizz = eventi.mxx / ni, eventi.myy / ni, eventi.mzz / ni
    # mixy, mixz, miyz = eventi.mxy / ni, eventi.mxz / ni, eventi.myz / ni
    # mjxx, mjyy, mjzz = eventj.mxx / nj, eventj.myy / nj, eventj.mzz / nj
    # mjxy, mjxz, mjyz = eventj.mxy / nj, eventj.mxz / nj, eventj.myz / nj

    d = abs(eventi.mxx - eventj.mxx) + abs(eventi.myy - eventj.myy) + \
        abs(eventi.mzz - eventj.mzz) + abs(eventi.mxy - eventj.mxy) + \
        abs(eventi.mxz - eventj.mxz) + abs(eventi.myz - eventj.myz)

    d = 0.5 * math.sqrt(d)

    return d


def get_distance_mt_cos(eventi, eventj):
    '''
    Inner product among two moment tensors.

    According to Willemann 1993; and later to Tape & Tape, normalization in
    R^9 to ensure innerproduct between -1 and +1.
    '''

    ni = math.sqrt(
        eventi.mxx * eventi.mxx +
        eventi.myy * eventi.myy +
        eventi.mzz * eventi.mzz +
        2. * eventi.mxy * eventi.mxy +
        2. * eventi.mxz * eventi.mxz +
        2. * eventi.myz * eventi.myz)

    nj = math.sqrt(
        eventj.mxx * eventj.mxx +
        eventj.myy * eventj.myy +
        eventj.mzz * eventj.mzz +
        2. * eventj.mxy * eventj.mxy +
        2. * eventj.mxz * eventj.mxz +
        2. * eventj.myz * eventj.myz)

    nc = ni * nj
    innerproduct = (
        eventi.mxx * eventj.mxx +
        eventi.myy * eventj.myy +
        eventi.mzz * eventj.mzz +
        2. * eventi.mxy * eventj.mxy +
        2. * eventi.mxz * eventj.mxz +
        2. * eventi.myz * eventj.myz) / nc

    if innerproduct >= 1.0:
        innerproduct = 1.0
    elif innerproduct <= -1.0:
        innerproduct = -1.0

    d = 0.5 * (1 - innerproduct)

    return d


def get_distance_mt_weighted_cos(eventi, eventj, ws):
    '''
    Weighted moment tensor distance.

    According to Cesca et al. 2014 GJI
    '''
    ni = math.sqrt(
        (ws[0] * eventi.mxx)**2 +
        (ws[1] * eventi.mxy)**2 +
        (ws[2] * eventi.myy)**2 +
        (ws[3] * eventi.mxz)**2 +
        (ws[4] * eventi.myz)**2 +
        (ws[5] * eventi.mzz)**2)

    nj = math.sqrt(
        (ws[0] * eventj.mxx)**2 +
        (ws[1] * eventj.mxy)**2 +
        (ws[2] * eventj.myy)**2 +
        (ws[3] * eventj.mxz)**2 +
        (ws[4] * eventj.myz)**2 +
        (ws[5] * eventj.mzz)**2)

    nc = ni * nj
    innerproduct = (
        ws[0] * ws[0] * eventi.mxx * eventj.mxx +
        ws[1] * ws[1] * eventi.mxy * eventj.mxy +
        ws[2] * ws[2] * eventi.myy * eventj.myy +
        ws[3] * ws[3] * eventi.mxz * eventj.mxz +
        ws[4] * ws[4] * eventi.myz * eventj.myz +
        ws[5] * ws[5] * eventi.mzz * eventj.mzz) / nc

    if innerproduct >= 1.0:
        innerproduct = 1.0

    elif innerproduct <= -1.0:
        innerproduct = -1.0

    d = 0.5 * (1.0 - innerproduct)

    return d


def get_distance_dc(eventi, eventj):
    '''Normalized Kagan angle distance among DC components of moment tensors.
       Based on Kagan, Y. Y., 1991, GJI
    '''

    mti = eventi.moment_tensor
    mtj = eventj.moment_tensor

    d = moment_tensor.kagan_angle(mti, mtj) / 120.
    if d > 1.:
        d = 1.

    return d


def get_distance_hypo(eventi, eventj):
    '''
    Normalized Euclidean hypocentral distance, assuming flat earth to combine
    epicentral distance and depth difference.

    The normalization assumes largest considered distance is 1000 km.
    '''
    maxdist_km = 1000.
    a_lats, a_lons, b_lats, b_lons = \
        eventi.north, eventi.east, eventj.north, eventj.east

    a_dep, b_dep = eventi.down, eventj.down

    if (a_lats == b_lats) and (a_lons == b_lons) and (a_dep == b_dep):
        d = 0.
    else:
        distance_m = orthodrome.distance_accurate50m_numpy(
            a_lats, a_lons, b_lats, b_lons)

        distance_km = distance_m / 1000.
        ddepth = abs(eventi.down - eventj.down)
        hypo_distance_km = math.sqrt(
            distance_km * distance_km + ddepth * ddepth)

        # maxdist = float(inv_param['EUCLIDEAN_MAX'])

        d = hypo_distance_km / maxdist_km
        if d >= 1.:
            d = 1.

    return d


def get_distance_epi(eventi, eventj):
    '''Normalized Euclidean epicentral distance.
       The normalization assumes largest considered distance is 1000 km.
    '''
    maxdist_km = 1000.

    a_lats, a_lons, b_lats, b_lons = \
        eventi.north, eventi.east, eventj.north, eventj.east

    a_dep, b_dep = eventi.down, eventj.down

    if (a_lats == b_lats) and (a_lons == b_lons) and (a_dep == b_dep):
        d = 0.
    else:
        distance_m = orthodrome.distance_accurate50m_numpy(
            a_lats, a_lons, b_lats, b_lons)

        distance_km = distance_m / 1000.

        d = distance_km / maxdist_km
        if d >= 1.:
            d = 1.

    return d


def get_distance_mt_triangle_diagram(eventi, eventj):
    '''
    Scalar product among principal axes (?).
    '''

    mti = eventi.moment_tensor
    mtj = eventj.moment_tensor

    ti, pi, bi = mti.t_axis(), mti.p_axis(), mti.b_axis()
    deltabi = math.acos(abs(bi[2]))
    deltati = math.acos(abs(ti[2]))
    deltapi = math.acos(abs(pi[2]))
    tj, pj, bj = mtj.t_axis(), mtj.p_axis(), mtj.b_axis()
    deltabj = math.acos(abs(bj[2]))
    deltatj = math.acos(abs(tj[2]))
    deltapj = math.acos(abs(pj[2]))
    dotprod = deltabi * deltabj + deltati * deltatj + deltapi * deltapj

    if dotprod >= 1.:
        dotprod == 1.

    d = 1. - dotprod
    return d


metric_funcs = {
    'mt_l2norm': get_distance_mt_l2,
    'mt_l1norm': get_distance_mt_l1,
    'mt_cos': get_distance_mt_cos,
    'mt_weighted_cos': get_distance_mt_weighted_cos,
    'mt_principal_axis': get_distance_mt_triangle_diagram,
    'kagan_angle': get_distance_dc,
    'hypocentral': get_distance_hypo,
    'epicentral': get_distance_epi,
}


metrics = sorted(metric_funcs.keys())


def get_distance(eventi, eventj, metric, **kwargs):
    '''
    Compute the normalized distance among two earthquakes, calling the function
    for the chosen metric definition.
    '''

    try:
        func = metric_funcs[metric]
    except KeyError:
        raise GrondError('unknown metric: %s' % metric)

    return func(eventi, eventj)


def compute_similarity_matrix(events, metric):
    '''
    Compute and return a similarity matrix for all event pairs, according to
    the desired metric

    :param events: list of pyrocko events
    :param metric: metric type (string)

    :returns: similarity matrix as NumPy array
    '''

    nev = len(events)
    simmat = num.zeros((nev, nev), dtype=float)
    for i in range(len(events)):
        for j in range(i):
            d = get_distance(events[i], events[j], metric)
            simmat[i, j] = d
            simmat[j, i] = d

    return simmat


def load_similarity_matrix(fname):
    '''
    Load a binary similarity matrix from file
    '''

    simmat = num.load(fname)
    return simmat
