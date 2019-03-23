import logging

from pyrocko import gf, util
from pyrocko.guts import String, Float, Dict, Int

from grond.meta import expand_template, Parameter, \
    has_get_plot_classes

from ..base import Problem, ProblemConfig

guts_prefix = 'grond'
logger = logging.getLogger('grond.problems.volume_point')
km = 1e3
as_km = dict(scale_factor=km, scale_unit='km')
as_km3 = dict(scale_factor=km**3, scale_unit='km^3')


class VolumePointProblemConfig(ProblemConfig):
    ranges = Dict.T(String.T(), gf.Range.T())
    distance_min = Float.T(default=0.0)
    nthreads = Int.T(default=1)

    def get_problem(self, event, target_groups, targets):
        if event.depth is None:
            event.depth = 0.

        base_source = gf.ExplosionSource.from_pyrocko_event(event)
        base_source.stf = gf.HalfSinusoidSTF(duration=event.duration or 0.0)

        subs = dict(
            event_name=event.name,
            event_time=util.time_to_str(event.time))

        problem = VolumePointProblem(
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
class VolumePointProblem(Problem):

    problem_parameters = [
        Parameter('north_shift', 'm', label='Northing', **as_km),
        Parameter('east_shift', 'm', label='Easting', **as_km),
        Parameter('depth', 'm', label='Depth', **as_km),
        Parameter('volume_change', 'm^3', label='Volume Change', **as_km3)
    ]

    problem_waveform_parameters = [
            Parameter('time', 's', label='Time'),
            Parameter('duration', 's', label='Duration')
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

        stf = None
        if self.has_waveforms:
            stf = gf.HalfSinusoidSTF(duration=float(d.duration))

        source = self.base_source.clone(stf=stf, **p)
        return source

    @classmethod
    def get_plot_classes(cls):
        from . import plot
        plots = super(VolumePointProblem, cls).get_plot_classes()
        plots.extend([plot.VolumePointLocationPlot])
        return plots
