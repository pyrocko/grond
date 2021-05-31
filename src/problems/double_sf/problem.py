import numpy as num
import logging

from pyrocko import gf, util
from pyrocko.guts import String, Float, Dict, StringChoice, Int

from grond.meta import Forbidden, expand_template, Parameter, \
    has_get_plot_classes

from ..base import Problem, ProblemConfig

guts_prefix = 'grond'
logger = logging.getLogger('grond.problems.double_sf.problem')
km = 1e3
as_km = dict(scale_factor=km, scale_unit='km')


class DoubleSFProblemConfig(ProblemConfig):

    ranges = Dict.T(String.T(), gf.Range.T())
    distance_min = Float.T(default=0.0)
    nthreads = Int.T(default=1)
    force_directions = StringChoice.T(
        choices=('off', 'unidirectional', 'counterdirectional'),
        default='off')

    def get_problem(self, event, target_groups, targets):
        if event.depth is None:
            event.depth = 0.

        base_source = gf.DoubleSFSource.from_pyrocko_event(event)

        base_source.stf1 = gf.HalfSinusoidSTF(duration=event.duration or 0.0)
        base_source.stf2 = gf.HalfSinusoidSTF(duration=event.duration or 0.0)

        subs = dict(
            event_name=event.name,
            event_time=util.time_to_str(event.time))

        problem = DoubleSFProblem(
            name=expand_template(self.name_template, subs),
            base_source=base_source,
            target_groups=target_groups,
            targets=targets,
            ranges=self.ranges,
            distance_min=self.distance_min,
            norm_exponent=self.norm_exponent,
            nthreads=self.nthreads,
            force_directions=self.force_directions)

        return problem


@has_get_plot_classes
class DoubleSFProblem(Problem):

    problem_parameters = [
        Parameter('time', 's', label='Time'),
        Parameter('north_shift', 'm', label='Northing', **as_km),
        Parameter('east_shift', 'm', label='Easting', **as_km),
        Parameter('depth', 'm', label='Depth', **as_km),
        Parameter('force', 'N', label='$||F||$'),
        Parameter('rfn1', '', label='$rF_{n1}$'),
        Parameter('rfe1', '', label='$rF_{e1}$'),
        Parameter('rfd1', '', label='$rF_{d1}$'),
        Parameter('rfn2', '', label='$rF_{n2}$'),
        Parameter('rfe2', '', label='$rF_{e2}$'),
        Parameter('rfd2', '', label='$rF_{d2}$'),
        Parameter('delta_time', 's', label='$\\Delta$ Time'),
        Parameter('delta_depth', 'm', label='$\\Delta$ Depth'),
        Parameter('azimuth', 'deg', label='Azimuth'),
        Parameter('distance', 'm', label='Distance'),
        Parameter('mix', label='Mix'),
        Parameter('duration1', 's', label='Duration 1'),
        Parameter('duration2', 's', label='Duration 2')]

    dependants = [
        Parameter('fn1', 'N', label='$F_{n1}$'),
        Parameter('fe1', 'N', label='$F_{e1}$'),
        Parameter('fd1', 'N', label='$F_{d1}$'),
        Parameter('fn2', 'N', label='$F_{n2}$'),
        Parameter('fe2', 'N', label='$F_{e2}$'),
        Parameter('fd2', 'N', label='$F_{d2}$')]

    distance_min = Float.T(default=0.0)
    force_directions = String.T()

    def __init__(self, **kwargs):
        Problem.__init__(self, **kwargs)
        self.deps_cache = {}

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
        cache = self.deps_cache
        if xs.ndim == 1:
            return self.make_dependant(xs[num.newaxis, :], pname)[0]

        if pname not in self.dependant_names:
            raise KeyError(pname)

        y = num.zeros(xs.shape[0])
        for i, x in enumerate(xs):
            k = tuple(x.tolist())
            if k not in cache:
                source = self.get_source(x)
                cache[k] = source

            source = cache[k]

            y[i] = getattr(source, pname)

        return y

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

        return x.tolist()

    def preconstrain(self, x):
        source = self.get_source(x)
        if any(self.distance_min > source.distance_to(t)
               for t in self.targets):
            raise Forbidden()

        if self.force_directions == 'unidirectional':
            idx_rfn2 = self.get_parameter_index('rfn2')
            idx_rfe2 = self.get_parameter_index('rfe2')
            idx_rfd2 = self.get_parameter_index('rfd2')

            x[idx_rfn2] = source.rfn1
            x[idx_rfe2] = source.rfe1
            x[idx_rfd2] = source.rfd1

        elif self.force_directions == 'counterdirectional':
            idx_rfn2 = self.get_parameter_index('rfn2')
            idx_rfe2 = self.get_parameter_index('rfe2')
            idx_rfd2 = self.get_parameter_index('rfd2')

            x[idx_rfn2] = -source.rfn1
            x[idx_rfe2] = -source.rfe1
            x[idx_rfd2] = -source.rfd1

        return num.array(x, dtype=num.float)

    def get_dependant_bounds(self):
        range_force = self.ranges['force']
        out = [
            (-range_force.stop, range_force.stop),
            (-range_force.stop, range_force.stop),
            (-range_force.stop, range_force.stop),
            (-range_force.stop, range_force.stop),
            (-range_force.stop, range_force.stop),
            (-range_force.stop, range_force.stop)]

        return out

    @classmethod
    def get_plot_classes(cls):
        from . import plot
        from ..singleforce import plot as sfplot
        plots = super(DoubleSFProblem, cls).get_plot_classes()
        plots.extend([sfplot.SFLocationPlot, plot.SFForcePlot])
        return plots


__all__ = '''
    DoubleSFProblem
    DoubleSFProblemConfig
'''.split()
