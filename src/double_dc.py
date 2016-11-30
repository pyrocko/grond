import logging

import numpy as num

from pyrocko import gf, util
from pyrocko.guts import Float, String, Dict, List, Int
from grond import core

guts_prefix = 'grond'


logger = logging.getLogger('grond.double_dc')

km = 1000.

as_km = dict(scale_factor=km, scale_unit='km')


class DoubleDCProblem(core.Problem):

    parameters = [
        core.Parameter('time', 's', label='Time'),
        core.Parameter('north_shift', 'm', label='Northing', **as_km),
        core.Parameter('east_shift', 'm', label='Easting', **as_km),
        core.Parameter('depth', 'm', label='Depth', **as_km),
        core.Parameter('magnitude', label='Magnitude'),
        core.Parameter('strike1', 'deg', label='Strike 1'),
        core.Parameter('dip1', 'deg', label='Dip 1'),
        core.Parameter('rake1', 'deg', label='Rake 1'),
        core.Parameter('strike2', 'deg', label='Strike 2'),
        core.Parameter('dip2', 'deg', label='Dip 2'),
        core.Parameter('rake2', 'deg', label='Rake 2'),
        core.Parameter('delta_time', 's', label='$\Delta$ Time'),
        core.Parameter('delta_depth', 'm', label='$\Delta$ Depth'),
        core.Parameter('azimuth', 'deg', label='Azimuth'),
        core.Parameter('distance', 'm', label='Distance'),
        core.Parameter('mix', label='Mix'),
        core.Parameter('duration1', 's', label='Duration 1'),
        core.Parameter('duration2', 's', label='Duration 2')]

    dependants = []

    targets = List.T(core.MisfitTarget.T())

    ranges = Dict.T(String.T(), gf.Range.T())
    distance_min = Float.T(default=0.0)
    nbootstrap = Int.T(default=100)

    def unpack(self, x):
        d = self.parameter_dict(x)
        p = {}
        for k in self.base_source.keys():
            if k in d:
                p[k] = float(
                    self.ranges[k].make_relative(self.base_source[k], d[k]))

        stf1 = gf.HalfSinusoidSTF(duration=float(d.duration1))
        stf2 = gf.HalfSinusoidSTF(duration=float(d.duration2))

        source = self.base_source.clone(stf1=stf1, stf2=stf2, **p)
        return source

    def make_dependant(self, xs, pname):
        if xs.ndim == 1:
            return self.make_dependant(xs[num.newaxis, :], pname)[0]

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
        x = num.array([
            source.time - self.base_source.time,
            source.north_shift,
            source.east_shift,
            source.depth,
            source.magnitude,
            source.strike1,
            source.dip1,
            source.rake1,
            source.strike2,
            source.dip2,
            source.rake2,
            source.delta_time,
            source.delta_depth,
            source.azimuth,
            source.distance,
            source.mix,
            source.stf1.duration,
            source.stf2.duration], dtype=num.float)

        return x

    def random_uniform(self, xbounds):
        x = num.zeros(self.nparameters)
        for i in xrange(self.nparameters):
            x[i] = num.random.uniform(xbounds[i, 0], xbounds[i, 1])

        return x.tolist()

    def preconstrain(self, x):
        source = self.unpack(x)
        if any(self.distance_min > source.distance_to(t)
               for t in self.targets):
            raise core.Forbidden()

        return num.array(x, dtype=num.float)

    def bounds(self):
        out = []
        for p in self.parameters:
            r = self.ranges[p.name]
            out.append((r.start, r.stop))

        return out

    def dependant_bounds(self):
        out = []

        return out

    def evaluate(self, x, result_mode='sparse'):
        source = self.unpack(x)
        engine = self.get_engine()

        for target in self.targets:
            target.set_result_mode(result_mode)

        resp = engine.process(source, self.targets)

        data = []
        results = []
        for target, result in zip(self.targets, resp.results_list[0]):
            if isinstance(result, gf.SeismosizerError):
                logger.debug(
                    '%s.%s.%s.%s: %s' % (target.codes + (str(result),)))

                data.append((None, None))
                results.append(result)
            else:
                data.append((result.misfit_value, result.misfit_norm))
                results.append(result)

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
            ws[mask] = 1.0 / num.sqrt(num.nansum(ns[mask]**2))

        return ws

    def inter_group_weights2(self, ns):
        group, ngroups = self.get_group_mask()
        ws = num.zeros(ns.shape)
        for igroup in xrange(ngroups):
            mask = group == igroup
            ws[:, mask] = (1.0 / num.sqrt(
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


class DoubleDCProblemConfig(core.ProblemConfig):

    ranges = Dict.T(String.T(), gf.Range.T())
    distance_min = Float.T(default=0.0)
    nbootstrap = Int.T(default=100)

    def get_problem(self, event, targets):
        if event.depth is None:
            event.depth = 0.

        base_source = gf.DoubleDCSource.from_pyrocko_event(event)
        base_source.stf = gf.HalfSinusoidSTF(duration=event.duration or 0.0)

        subs = dict(
            event_name=event.name,
            event_time=util.time_to_str(event.time))

        problem = DoubleDCProblem(
            name=core.expand_template(self.name_template, subs),
            apply_balancing_weights=self.apply_balancing_weights,
            base_source=base_source,
            targets=targets,
            ranges=self.ranges,
            distance_min=self.distance_min,
            nbootstrap=self.nbootstrap)

        return problem


__all__ = [
    'DoubleDCProblem',
    'DoubleDCProblemConfig',
]
