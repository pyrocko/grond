from __future__ import print_function
import nose.tools as t

import numpy as num

from numpy.testing import assert_almost_equal as assert_ae
from pyrocko import gf
from grond.toy import scenario, ToyProblem


def test_combine_misfits():
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

    misfitss = p.evaluate_many(xg)
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

    # gms_2[imodel, ibootstrap]
    # gms_2_contrib[imodel, ibootstrap, itarget]

    for ix, x in enumerate(xg):
        misfits = p.evaluate(x)
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
