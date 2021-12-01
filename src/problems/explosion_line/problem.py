import numpy as num
import logging

from pyrocko import gf, util
from pyrocko.guts import String, Float, Dict, Int

from grond.meta import expand_template, Parameter, has_get_plot_classes

from ..base import Problem, ProblemConfig

logger = logging.getLogger(__name__)
km = 1e3
as_km = dict(scale_factor=km, scale_unit='km')


class ExplosionLineProblemConfig(ProblemConfig):

    ranges = Dict.T(String.T(), gf.Range.T())
    distance_min = Float.T(default=0.)
    subsource_oversampling = Int.T(default=1)

    def get_problem(self, event, target_groups, targets):

        base_source = gf.ExplosionLineSource.from_pyrocko_event(
            event,
            subsource_oversampling=self.subsource_oversampling
        )

        subs = dict(
            event_name=event.name,
            event_time=util.time_to_str(event.time))

        problem = ExplosionLineProblem(
            name=expand_template(self.name_template, subs),
            base_source=base_source,
            distance_min=self.distance_min,
            target_groups=target_groups,
            targets=targets,
            ranges=self.ranges,
            norm_exponent=self.norm_exponent
        )

        return problem


@has_get_plot_classes
class ExplosionLineProblem(Problem):

    problem_parameters = [
        Parameter('time', 's', label='Time'),
        Parameter('east_shift', 'm', label='Easting', **as_km),
        Parameter('north_shift', 'm', label='Northing', **as_km),
        Parameter('depth', 'm', label='Depth', **as_km),

        Parameter('end_time', 's', label='Time'),
        Parameter('end_east_shift', 'm', label='Easting', **as_km),
        Parameter('end_north_shift', 'm', label='Northing', **as_km),
        Parameter('end_depth', 'm', label='Depth', **as_km),
    ]

    dependants = [
        Parameter('duration', 's', label='Duration'),
        Parameter('length', 'm', label='Length', **as_km),
        Parameter('velocity', 'm', label='Velocity'),
        Parameter('dip', u'\u00b0', label='Dip'),
        Parameter('strike', u'\u00b0', label='Strike')
    ]

    distance_min = Float.T(default=0.0)

    def pack(self, source):
        arr = self.get_parameter_array(source)
        for ip, p in enumerate(self.parameters):
            if p.name == 'time':
                arr[ip] -= self.base_source.time

        return arr

    def get_source(self, x):
        d = self.get_parameter_dict(x)

        p = {k: float(self.ranges[k].make_relative(self.base_source[k], d[k]))
             for k in self.base_source.keys() if k in d}

        source = self.base_source.clone(**p)
        return source

    def random_uniform(self, xbounds, rstate, fixed_magnitude=None):
        x = num.zeros(self.nparameters)
        for ip in range(self.nparameters):
            x[ip] = rstate.uniform(xbounds[ip, 0], xbounds[ip, 1])

        return x

    def preconstrain(self, x):
        return x

    @classmethod
    def get_plot_classes(cls):
        plots = super().get_plot_classes()
        return plots


__all__ = '''
    ExplosionLineProblem
    ExplosionLineProblemConfig
'''
