import numpy as num
import logging

from pyrocko import gf, util
from pyrocko.guts import String, Float, Dict, Int

from .base import Problem, ProblemConfig
from ..meta import Forbidden, expand_template, Parameter


guts_prefix = 'grond'
logger = logging.getLogger('grond.problems').getChild('double_dc')
km = 1e3
as_km = dict(scale_factor=km, scale_unit='km')


class DoubleDCProblemConfig(ProblemConfig):

    ranges = Dict.T(String.T(), gf.Range.T())
    distance_min = Float.T(default=0.0)

    def get_problem(self, event, targets):
        if event.depth is None:
            event.depth = 0.

        base_source = gf.DoubleDCSource.from_pyrocko_event(event)
        base_source.stf = gf.HalfSinusoidSTF(duration=event.duration or 0.0)

        subs = dict(
            event_name=event.name,
            event_time=util.time_to_str(event.time))

        problem = DoubleDCProblem(
            name=expand_template(self.name_template, subs),
            apply_balancing_weights=self.apply_balancing_weights,
            base_source=base_source,
            targets=targets,
            ranges=self.ranges,
            distance_min=self.distance_min,
            norm_exponent=self.norm_exponent)

        return problem


class DoubleDCProblem(Problem):

    problem_parameters = [
        Parameter('time', 's', label='Time'),
        Parameter('north_shift', 'm', label='Northing', **as_km),
        Parameter('east_shift', 'm', label='Easting', **as_km),
        Parameter('depth', 'm', label='Depth', **as_km),
        Parameter('magnitude', label='Magnitude'),
        Parameter('strike1', 'deg', label='Strike 1'),
        Parameter('dip1', 'deg', label='Dip 1'),
        Parameter('rake1', 'deg', label='Rake 1'),
        Parameter('strike2', 'deg', label='Strike 2'),
        Parameter('dip2', 'deg', label='Dip 2'),
        Parameter('rake2', 'deg', label='Rake 2'),
        Parameter('delta_time', 's', label='$\Delta$ Time'),
        Parameter('delta_depth', 'm', label='$\Delta$ Depth'),
        Parameter('azimuth', 'deg', label='Azimuth'),
        Parameter('distance', 'm', label='Distance'),
        Parameter('mix', label='Mix'),
        Parameter('duration1', 's', label='Duration 1'),
        Parameter('duration2', 's', label='Duration 2')]

    dependants = []
    distance_min = Float.T(default=0.0)

    def get_source(self, x):
        d = self.get_parameter_dict(x)
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
        for i in range(self.nparameters):
            x[i] = num.random.uniform(xbounds[i, 0], xbounds[i, 1])

        return x.tolist()

    def preconstrain(self, x):
        source = self.get_source(x)
        if any(self.distance_min > source.distance_to(t)
               for t in self.targets):
            raise Forbidden()

        return num.array(x, dtype=num.float)

    def evaluate(self, x, result_mode='sparse'):
        source = self.get_source(x)
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
        source = self.get_source(x)
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
