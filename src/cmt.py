import math
import logging

import numpy as num

from pyrocko import gf, moment_tensor as mtm, util
from pyrocko.guts import Float, String, Dict, List, Int, StringChoice
from grond import core

guts_prefix = 'grond'


logger = logging.getLogger('grond.cmt')

km = 1000.

as_km = dict(scale_factor=km, scale_unit='km')


class CMTProblem(core.Problem):

    parameters = [
        core.Parameter('time', 's', label='Time'),
        core.Parameter('north_shift', 'm', label='Northing', **as_km),
        core.Parameter('east_shift', 'm', label='Easting', **as_km),
        core.Parameter('depth', 'm', label='Depth', **as_km),
        core.Parameter('magnitude', label='Magnitude'),
        core.Parameter('rmnn', label='$m_{nn} / M_0$'),
        core.Parameter('rmee', label='$m_{ee} / M_0$'),
        core.Parameter('rmdd', label='$m_{dd} / M_0$'),
        core.Parameter('rmne', label='$m_{ne} / M_0$'),
        core.Parameter('rmnd', label='$m_{nd} / M_0$'),
        core.Parameter('rmed', label='$m_{ed} / M_0$'),
        core.Parameter('duration', 's', label='Duration')]

    dependants = [
        core.Parameter('rel_moment_iso', label='$M_{0}^{ISO}/M_{0}$'),
        core.Parameter('rel_moment_clvd', label='$M_{0}^{CLVD}/M_{0}$')]

    targets = List.T(core.MisfitTarget.T())

    ranges = Dict.T(String.T(), gf.Range.T())
    distance_min = Float.T(default=0.0)
    nbootstrap = Int.T(default=10)
    mt_type = StringChoice.T(
        default='full', choices=['full', 'deviatoric', 'dc'])

    def unpack(self, x):
        d = self.parameter_dict(x)
        rm6 = num.array([d.rmnn, d.rmee, d.rmdd, d.rmne, d.rmnd, d.rmed],
                        dtype=num.float)

        m0 = mtm.magnitude_to_moment(d.magnitude)
        m6 = rm6 * m0

        p = {}
        for k in self.base_source.keys():
            if k in d:
                p[k] = float(
                    self.ranges[k].make_relative(self.base_source[k], d[k]))

        stf = gf.HalfSinusoidSTF(duration=float(d.duration))

        source = self.base_source.clone(m6=m6, stf=stf, **p)
        return source

    def make_dependant(self, xs, pname):
        if xs.ndim == 1:
            return self.make_dependant(xs[num.newaxis, :], pname)[0]

        if pname in ('rel_moment_iso', 'rel_moment_clvd'):
            y = num.zeros(xs.shape[0])
            for i, x in enumerate(xs):
                source = self.unpack(x)
                mt = source.pyrocko_moment_tensor()
                res = mt.standard_decomposition()

                if pname == 'rel_moment_iso':
                    ratio_iso, m_iso = res[0][1:3]
                    y[i] = ratio_iso * num.sign(m_iso[0, 0])
                else:
                    ratio_clvd, m_clvd = res[2][1:3]
                    evals, evecs = mtm.eigh_check(m_clvd)
                    ii = num.argmax(num.abs(evals))
                    y[i] = ratio_clvd * num.sign(evals[ii])

            return y

        else:
            raise KeyError(pname)

    def extract(self, xs, i):
        if xs.ndim == 1:
            return self.extract(xs[num.newaxis, :], i)[0]

        if i < self.nparameters:
            return xs[:, i]
        else:
            return self.make_dependant(
                xs, self.dependants[i-self.nparameters].name)

    def pack(self, source):
        m6 = source.m6
        mt = source.pyrocko_moment_tensor()
        rm6 = m6 / mt.scalar_moment()

        x = num.array([
            source.time - self.base_source.time,
            source.north_shift,
            source.east_shift,
            source.depth,
            mt.moment_magnitude(),
            ] + rm6.tolist() + [source.stf.duration], dtype=num.float)

        return x

    def random_uniform(self, xbounds):
        x = num.zeros(self.nparameters)
        for i in [0, 1, 2, 3, 4, 11]:
            x[i] = num.random.uniform(xbounds[i, 0], xbounds[i, 1])

        x[5:11] = mtm.random_m6()

        return x.tolist()

    def preconstrain(self, x):

        d = self.parameter_dict(x)
        m6 = num.array([d.rmnn, d.rmee, d.rmdd, d.rmne, d.rmnd, d.rmed],
                       dtype=num.float)

        m9 = mtm.symmat6(*m6)
        if self.mt_type == 'deviatoric':
            trace_m = num.trace(m9)
            m_iso = num.diag([trace_m / 3., trace_m / 3., trace_m / 3.])
            m9 -= m_iso

        elif self.mt_type == 'dc':
            mt = mtm.MomentTensor(m=m9)
            m9 = mt.standard_decomposition()[1][2]

        m0_unscaled = math.sqrt(num.sum(m9.A**2)) / math.sqrt(2.)

        m9 /= m0_unscaled
        m6 = mtm.to6(m9)
        d.rmnn, d.rmee, d.rmdd, d.rmne, d.rmnd, d.rmed = m6
        x = self.parameter_array(d)

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

    def raise_invalid_norm_exponent(self):
        raise core.GrondError('invalid norm exponent' % self.norm_exponent)

    def get_sqr_sqrt(self):
        if self.norm_exponent == 2:
            def sqr(x):
                return x**2

            return sqr, num.sqrt

        elif self.norm_exponent == 1:
            def noop(x):
                return x

            return noop, num.abs

        else:
            self.raise_invalid_norm_exponent()

    def inter_group_weights(self, ns):
        sqr, sqrt = self.get_sqr_sqrt()

        group, ngroups = self.get_group_mask()

        ws = num.zeros(self.ntargets)
        for igroup in xrange(ngroups):
            mask = group == igroup
            ws[mask] = 1.0 / sqrt(num.nansum(sqr(ns[mask])))

        return ws

    def inter_group_weights2(self, ns):
        sqr, sqrt = self.get_sqr_sqrt()

        group, ngroups = self.get_group_mask()
        ws = num.zeros(ns.shape)
        for igroup in xrange(ngroups):
            mask = group == igroup
            ws[:, mask] = (1.0 / sqrt(
                num.nansum(sqr(ns[:, mask]), axis=1)))[:, num.newaxis]

        return ws

    def bootstrap_misfit(self, ms, ns, ibootstrap=None):
        sqr, sqrt = self.get_sqr_sqrt()

        w = self.get_bootstrap_weights(ibootstrap) * \
            self.get_target_weights() * self.inter_group_weights(ns)

        if ibootstrap is None:
            return sqrt(
                num.nansum(sqr(w*ms[num.newaxis, :]), axis=1) /
                num.nansum(sqr(w*ns[num.newaxis, :]), axis=1))
        else:
            return sqrt(num.nansum(sqr(w*ms)) / num.nansum(sqr(w*ns)))

    def bootstrap_misfits(self, misfits, ibootstrap):
        sqr, sqrt = self.get_sqr_sqrt()

        w = self.get_bootstrap_weights(ibootstrap)[num.newaxis, :] * \
            self.get_target_weights()[num.newaxis, :] * \
            self.inter_group_weights2(misfits[:, :, 1])

        bms = sqrt(num.nansum(sqr(w*misfits[:, :, 0]), axis=1) /
                   num.nansum(sqr(w*misfits[:, :, 1]), axis=1))
        return bms

    def global_misfit(self, ms, ns):
        sqr, sqrt = self.get_sqr_sqrt()

        ws = self.get_target_weights() * self.inter_group_weights(ns)
        m = sqrt(num.nansum(sqr(ws*ms)) / num.nansum(sqr(ws*ns)))
        return m

    def global_misfits(self, misfits):
        sqr, sqrt = self.get_sqr_sqrt()
        ws = self.get_target_weights()[num.newaxis, :] * \
            self.inter_group_weights2(misfits[:, :, 1])
        gms = sqrt(num.nansum(sqr(ws*misfits[:, :, 0]), axis=1) /
                   num.nansum(sqr(ws*misfits[:, :, 1]), axis=1))
        return gms

    def global_contributions(self, misfits):
        sqr, sqrt = self.get_sqr_sqrt()
        ws = self.get_target_weights()[num.newaxis, :] * \
            self.inter_group_weights2(misfits[:, :, 1])

        gcms = sqr(ws*misfits[:, :, 0]) / \
            num.nansum(sqr(ws*misfits[:, :, 1]), axis=1)[:, num.newaxis]

        return gcms


class CMTProblemConfig(core.ProblemConfig):

    ranges = Dict.T(String.T(), gf.Range.T())
    distance_min = Float.T(default=0.0)
    nbootstrap = Int.T(default=10)
    mt_type = StringChoice.T(choices=['full', 'deviatoric', 'dc'])

    def get_problem(self, event, targets):
        if event.depth is None:
            event.depth = 0.

        base_source = gf.MTSource.from_pyrocko_event(event)
        base_source.stf = gf.HalfSinusoidSTF(duration=event.duration or 0.0)

        subs = dict(
            event_name=event.name,
            event_time=util.time_to_str(event.time))

        problem = CMTProblem(
            name=core.expand_template(self.name_template, subs),
            apply_balancing_weights=self.apply_balancing_weights,
            norm_exponent=self.norm_exponent,
            base_source=base_source,
            targets=targets,
            ranges=self.ranges,
            distance_min=self.distance_min,
            nbootstrap=self.nbootstrap,
            mt_type=self.mt_type)

        return problem


__all__ = [
    'CMTProblem',
    'CMTProblemConfig',
]
