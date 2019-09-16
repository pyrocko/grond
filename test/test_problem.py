from __future__ import print_function
import os.path as op
import nose.tools as t

import numpy as num

from numpy.testing import assert_almost_equal as assert_ae
from pyrocko import gf
from grond.toy import scenario, ToyProblem, ToyTarget, ToySource


def test_combine_misfits(dump=False, reference=None):
    source, targets = scenario('wellposed', 'noisefree')

    p = ToyProblem(
        name='toy_problem',
        ranges={
            'north': gf.Range(start=-10., stop=10.),
            'east': gf.Range(start=-10., stop=10.),
            'depth': gf.Range(start=0., stop=10.)},
        base_source=source,
        targets=targets)

    ngx, ngy, ngz = 11, 11, 11
    xg = num.zeros((ngz*ngy*ngx, 3))

    xbounds = p.get_parameter_bounds()
    cx = num.linspace(xbounds[0][0], xbounds[0][1], ngx)
    cy = num.linspace(xbounds[1][0], xbounds[1][1], ngy)
    cz = num.linspace(xbounds[2][0], xbounds[2][1], ngz)

    xg[:, 0] = num.tile(cx, ngy*ngz)
    xg[:, 1] = num.tile(num.repeat(cy, ngx), ngz)
    xg[:, 2] = num.repeat(cz, ngx*ngy)

    misfitss = p.misfits_many(xg)
    # misfitss[imodel, itarget, 0], misfitss[imodel, itarget, 1]
    gms = p.combine_misfits(misfitss)
    gms_contrib = p.combine_misfits(misfitss, get_contributions=True)

    # gms[imodel]
    # gms_contrib[imodel, itarget]

    bweights = num.ones((2, p.ntargets))
    gms_2 = p.combine_misfits(misfitss, extra_weights=bweights)
    gms_2_contrib = p.combine_misfits(
        misfitss,
        extra_weights=bweights,
        get_contributions=True)

    if dump:
        num.savez(dump, gms, gms_contrib, gms_2)

    # gms_2[imodel, ibootstrap]
    # gms_2_contrib[imodel, ibootstrap, itarget]

    for ix, x in enumerate(xg):
        misfits = p.misfits(x)
        # misfits[itarget, 0], misfits[itarget, 1]
        gm = p.combine_misfits(misfits)
        # gm is scalar
        t.assert_equal(gm, gms[ix])

        gm_contrib = p.combine_misfits(
            misfits,
            get_contributions=True)

        assert_ae(gms_contrib[ix, :], gm_contrib)

        gm_2 = p.combine_misfits(misfits, extra_weights=bweights)

        assert gm_2[0] == gm
        assert gm_2[1] == gm
        assert gms_2[ix, 0] == gm
        assert gms_2[ix, 1] == gm

        gm_2_contrib = p.combine_misfits(
            misfits, extra_weights=bweights,
            get_contributions=True)

        assert_ae(gm_2_contrib[0, :], gm_contrib)
        assert_ae(gm_2_contrib[1, :], gm_contrib)
        assert_ae(gms_2_contrib[ix, 0, :], gm_contrib)
        assert_ae(gms_2_contrib[ix, 1, :], gm_contrib)

    if reference:
        ref_data = num.load(reference)

        assert_ae(ref_data['arr_0'], gms)
        assert_ae(ref_data['arr_1'], gms_contrib)
        assert_ae(ref_data['arr_2'], gms_2)


def test_combine_reference():
    ref_fn = op.join(op.dirname(__file__), 'combine_misfits-v1.2.0.npz')
    test_combine_misfits(reference=ref_fn)


def test_combine_covariance():
    nmodels = 1  # noqa
    nmisfits = 15

    target = ToyTarget(
        path='t_corr',
        north=4.,
        east=3.,
        depth=0.,
        obs_distance=0.,
        nmisfits=nmisfits)

    source = ToySource(
        north=0,
        east=0,
        depth=5.)

    p = ToyProblem(
        name='toy_problem',
        ranges={
            'north': gf.Range(start=-10., stop=10.),
            'east': gf.Range(start=-10., stop=10.),
            'depth': gf.Range(start=0., stop=10.)},
        base_source=source,
        targets=[target])

    rstate = num.random.RandomState(123)

    misfits = num.zeros((nmisfits, 2))
    misfits[:, 0] = rstate.normal(size=nmisfits)
    misfits[:, 1] = rstate.normal(size=nmisfits) + 3.

    weights = rstate.normal(size=nmisfits)

    corr = num.zeros((nmisfits, nmisfits))
    num.fill_diagonal(corr, weights)

    extra_residuals = num.array([])

    res_weights = p.combine_misfits(
        misfits,
        extra_weights=weights[num.newaxis, :],
        extra_residuals=extra_residuals)

    res_corr = p.combine_misfits(
        misfits,
        extra_correlated_weights={0: corr},
        extra_residuals=extra_residuals)

    num.testing.assert_almost_equal(res_weights, res_corr)


def dump_combine_misfits():
    test_combine_misfits(dump='combined_misfits.npz')
