import numpy as num
import logging

from pyrocko import gf, util
from pyrocko.guts import String, Float, Dict, Int

from grond.meta import Forbidden, expand_template, Parameter, \
    has_get_plot_classes

from ..base import Problem, ProblemConfig


guts_prefix = 'grond'
logger = logging.getLogger('grond.problems.double_dc.problem')
km = 1e3
as_km = dict(scale_factor=km, scale_unit='km')


class DoubleDCProblemConfig(ProblemConfig):

    ranges = Dict.T(String.T(), gf.Range.T())
    distance_min = Float.T(default=0.0)
    nthreads = Int.T(default=1)

    def get_problem(self, event, target_groups, targets):
        if event.depth is None:
            event.depth = 0.

        base_source = gf.DoubleDCSource.from_pyrocko_event(event)
        base_source.stf = gf.HalfSinusoidSTF(duration=event.duration or 0.0)

        subs = dict(
            event_name=event.name,
            event_time=util.time_to_str(event.time))

        problem = DoubleDCProblem(
            name=expand_template(self.name_template, subs),
            base_source=base_source,
            target_groups=target_groups,
            targets=targets,
            ranges=self.ranges,
            distance_min=self.distance_min,
            norm_exponent=self.norm_exponent,
            nthreads=self.nthreads)

        return problem


@has_get_plot_classes
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
        Parameter('delta_time', 's', label='$\\Delta$ Time'),
        Parameter('delta_depth', 'm', label='$\\Delta$ Depth'),
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

    def pack(self, source):
        arr = self.get_parameter_array(source)
        for ip, p in enumerate(self.parameters):
            if p.name == 'time':
                arr[ip] -= self.base_source.time
            if p.name == 'duration1':
                arr[ip] = source.stf1.duration if source.stf1 else 0.0
            if p.name == 'duration2':
                arr[ip] = source.stf2.duration if source.stf2 else 0.0
        return arr

    def random_uniform(self, xbounds, rstate, fixed_magnitude=None):
        x = num.zeros(self.nparameters)
        for i in range(self.nparameters):
            x[i] = rstate.uniform(xbounds[i, 0], xbounds[i, 1])

        if fixed_magnitude is not None:
            x[4] = fixed_magnitude

        return x.tolist()

    def preconstrain(self, x):
        source = self.get_source(x)
        if any(self.distance_min > source.distance_to(t)
               for t in self.targets):
            raise Forbidden()

        return num.array(x, dtype=num.float)

    @classmethod
    def get_plot_classes(cls):
        from .. import plot
        plots = super(DoubleDCProblem, cls).get_plot_classes()
        plots.extend([plot.HudsonPlot, plot.MTDecompositionPlot,
                      plot.MTLocationPlot])
        return plots


__all__ = '''
    DoubleDCProblem
    DoubleDCProblemConfig
'''.split()
