import math
import logging

import numpy as num

from pyrocko import gf, moment_tensor as mtm
from pyrocko.guts import Float, String, Dict, List, Int
from grond import core

guts_prefix = 'grond'

logger = logging.getLogger('grond.cmt')

km = 1000.

as_km = dict(scale_factor=km, scale_unit='km')


class RectangularProblem(core.Problem):

    parameters = [
        core.Parameter('north_shift', 'm', label='Northing', **as_km),
        core.Parameter('east_shift', 'm', label='Easting', **as_km),
        core.Parameter('depth', 'm', label='Depth', **as_km),
        core.Parameter('length', 'm', label='Length', **as_km),
        core.Parameter('width', 'm', label='Width', **as_km),
        core.Parameter('dip', 'deg', label='Dip'),
        core.Parameter('strike', 'deg', label='Strike'),
        core.Parameter('rake', 'deg', label='Rake'),
        core.Parameter('slip', 'm', label='Slip'),
        ]

    dependants = []

    targets = List.T(gf.Target.T())

    ranges = Dict.T(String.T(), gf.Range.T())
    distance_min = Float.T(default=0.0)
    nbootstrap = Int.T(default=10)

    def pack(self, source):
        return self.parameter_array(source)

    def unpack(self, x):
        source = self.base_source.clone(**self.parameter_dict(x))
        return source

    def extract(self, xs, i):
        if xs.ndim == 1:
            return self.extract(xs[num.newaxis, :], i)[0]

        return xs[:, i]

    def random_uniform(self, xbounds):
        x = num.zeros(self.nparameters)
        for i in [0, 1, 2, 3, 4, 11]:
            x[i] = num.random.uniform(xbounds[i, 0], xbounds[i, 1])

        x[5:11] = mtm.random_m6()

        return x.tolist()

    def preconstrain(self, x):
        source = self.unpack(x)
        if any(self.distance_min > source.distance_to(t)
               for t in self.targets):
            raise core.Forbidden()
        return x

    def bounds(self):
        out = []
        for p in self.parameters:
            r = self.ranges[p.name]
            out.append((r.start, r.stop))
        return out

    def dependant_bounds(self):
        print('dependent_bounds')
        out = [
            (-1., 1.),
            (-1., 1.)]
        return out

    def evaluate(self, x, result_mode='sparse', mask=None):
        source = self.unpack(x)
        engine = self.get_engine()

        for target in self.targets:
            target.set_result_mode(result_mode)

        if mask is not None:
            assert len(mask) == len(self.targets)
            targets_ok = [
                target for (target, ok) in zip(self.targets, mask) if ok]
        else:
            targets_ok = self.targets

        resp = engine.process(source, targets_ok)

        if mask is not None:
            ires_ok = 0
            results = []
            for target, ok in zip(self.targets, mask):
                if ok:
                    results.append(resp.results_list[0][ires_ok])
                    ires_ok += 1
                else:
                    results.append(
                        gf.SeismosizerError(
                            'skipped because of previous failure'))
        else:
            results = list(resp.results_list[0])

        data = []
        for target, result in zip(self.targets, results):
            if isinstance(result, gf.SeismosizerError):
                logger.debug(
                    '%s.%s.%s.%s: %s' % (target.codes + (str(result),)))

                data.append((None, None))
            else:
                data.append((result.misfit_value, result.misfit_norm))

        ms, ns = num.array(data, dtype=num.float).T
        if result_mode == 'full':
            return ms, ns, results
        else:
            return ms, ns

    def forward(self, x):
        source = self.unpack(x)
        engine = self.get_engine()
        plain_targets = [target.get_plain_target() for target in self.targets]

        resp = engine.process(source, plain_targets)
        results = []
        for target, result in zip(self.targets, resp.results_list[0]):
            if isinstance(result, gf.SeismosizerError):
                logger.debug(
                    '%s.%s.%s.%s: %s' % (target.codes + (str(result),)))

                results.append(None)
            else:
                results.append(result)

        return results

    def get_target_weights(self):
        if self._target_weights is None:
            self._target_weights = num.array(
                [target.get_combined_weight(
                    apply_balancing_weights=self.apply_balancing_weights)
                 for target in self.targets], dtype=num.float)

        return self._target_weights

    def inter_group_weights(self, ns):
        group, ngroups = self.get_group_mask()

        ws = num.zeros(self.ntargets)
        for igroup in xrange(ngroups):
            mask = group == igroup
            ws[mask] = 1. / num.sqrt(num.nansum(ns[mask]**2))

        return ws

    def inter_group_weights2(self, ns):
        group, ngroups = self.get_group_mask()
        ws = num.zeros(ns.shape)
        for igroup in xrange(ngroups):
            mask = group == igroup
            ws[:, mask] = (1. / num.sqrt(
                num.nansum(ns[:, mask]**2, axis=1)))[:, num.newaxis]

        return ws

    def bootstrap_misfit(self, ms, ns, ibootstrap=None):
        w = self.get_bootstrap_weights(ibootstrap) * \
            self.get_target_weights() * self.inter_group_weights(ns)

        if ibootstrap is None:
            return num.sqrt(
                num.nansum((w*ms[num.newaxis, :])**2, axis=1) /
                num.nansum((w*ns[num.newaxis, :])**2, axis=1))
        else:
            return num.sqrt(num.nansum((w*ms)**2) / num.nansum((w*ns)**2))

    def bootstrap_misfits(self, misfits, ibootstrap):
        w = self.get_bootstrap_weights(ibootstrap)[num.newaxis, :] * \
            self.get_target_weights()[num.newaxis, :] * \
            self.inter_group_weights2(misfits[:, :, 1])

        bms = num.sqrt(num.nansum((w*misfits[:, :, 0])**2, axis=1) /
                       num.nansum((w*misfits[:, :, 1])**2, axis=1))
        return bms

    def global_misfit(self, ms, ns):
        ws = self.get_target_weights() * self.inter_group_weights(ns)
        m = num.sqrt(num.nansum((ws*ms)**2) / num.nansum((ws*ns)**2))
        return m

    def global_misfits(self, misfits):
        ws = self.get_target_weights()[num.newaxis, :] * \
            self.inter_group_weights2(misfits[:, :, 1])
        gms = num.sqrt(num.nansum((ws*misfits[:, :, 0])**2, axis=1) /
                       num.nansum((ws*misfits[:, :, 1])**2, axis=1))
        return gms

    def global_contributions(self, misfits):
        ws = self.get_target_weights()[num.newaxis, :] * \
            self.inter_group_weights2(misfits[:, :, 1])

        gcms = (ws*misfits[:, :, 0])**2 / \
            num.nansum((ws*misfits[:, :, 1])**2, axis=1)[:, num.newaxis]

        return gcms


class RectangularProblemConfig(core.ProblemConfig):

    ranges = Dict.T(String.T(), gf.Range.T())
    name = String.T()
    distance_min = Float.T(default=0.0)
    nbootstrap = Int.T(default=0)
    apply_balancing_weights = False

    def get_problem(self, rect_source, targets):
        base_source = gf.RecangularSource.load(rect_source)

        problem = RectangularProblem(
            name=core.expand_template(self.name_template, self.name),
            apply_balancing_weights=self.apply_balancing_weights,
            base_source=base_source,
            targets=targets,
            ranges=self.ranges,
            distance_min=self.distance_min,
            nbootstrap=self.nbootstrap,
            mt_type=self.mt_type)

        return problem


__all__ = [
    'RectangularProblem',
    'RectangularProblemConfig',
]
