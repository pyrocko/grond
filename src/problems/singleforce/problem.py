import numpy as num
import logging

from pyrocko import gf, util
from pyrocko.guts import String, Float, Dict, StringChoice, Int

from grond.meta import expand_template, Parameter, has_get_plot_classes

from ..base import Problem, ProblemConfig

guts_prefix = 'grond'
logger = logging.getLogger('grond.problems.singleforce.problem')
km = 1e3
as_km = dict(scale_factor=km, scale_unit='km')


class STFType(StringChoice):
    choices = ['HalfSinusoidSTF', 'ResonatorSTF']

    cls = {
        'HalfSinusoidSTF': gf.HalfSinusoidSTF,
        'ResonatorSTF': gf.ResonatorSTF}

    @classmethod
    def base_stf(cls, name):
        return cls.cls[name]()


class SFProblemConfig(ProblemConfig):

    ranges = Dict.T(String.T(), gf.Range.T())
    distance_min = Float.T(default=0.0)
    stf_type = STFType.T(default='HalfSinusoidSTF')
    nthreads = Int.T(default=1)

    def get_problem(self, event, target_groups, targets):
        if event.depth is None:
            event.depth = 0.

        base_source = gf.SFSource.from_pyrocko_event(event)

        stf = STFType.base_stf(self.stf_type)
        stf.duration = event.duration or 0.0

        base_source.stf = stf

        subs = dict(
            event_name=event.name,
            event_time=util.time_to_str(event.time))

        problem = SFProblem(
            name=expand_template(self.name_template, subs),
            base_source=base_source,
            target_groups=target_groups,
            targets=targets,
            ranges=self.ranges,
            distance_min=self.distance_min,
            stf_type=self.stf_type,
            norm_exponent=self.norm_exponent,
            nthreads=self.nthreads)

        return problem


@has_get_plot_classes
class SFProblem(Problem):

    problem_parameters = [
        Parameter('time', 's', label='Time'),
        Parameter('north_shift', 'm', label='Northing', **as_km),
        Parameter('east_shift', 'm', label='Easting', **as_km),
        Parameter('depth', 'm', label='Depth', **as_km),
        Parameter('fn', 'N', label='$F_{n}$'),
        Parameter('fe', 'N', label='$F_{e}$'),
        Parameter('fd', 'N', label='$F_{d}$')]

    problem_parameters_stf = {
        'HalfSinusoidSTF': [
            Parameter('duration', 's', label='Duration')],
        'ResonatorSTF': [
            Parameter('duration', 's', label='Duration'),
            Parameter('frequency', 'Hz', label='Frequency')]}

    dependants = []

    distance_min = Float.T(default=0.0)
    stf_type = STFType.T(default='HalfSinusoidSTF')

    def __init__(self, **kwargs):
        Problem.__init__(self, **kwargs)
        self.deps_cache = {}
        self.problem_parameters = self.problem_parameters \
            + self.problem_parameters_stf[self.stf_type]
        self._base_stf = STFType.base_stf(self.stf_type)

    def get_stf(self, d):
        d_stf = {}
        for p in self.problem_parameters_stf[self.stf_type]:
            d_stf[p.name] = float(d[p.name])

        return self._base_stf.clone(**d_stf)

    def get_source(self, x):
        d = self.get_parameter_dict(x)

        p = {}
        for k in self.base_source.keys():
            if k in d:
                p[k] = float(
                    self.ranges[k].make_relative(self.base_source[k], d[k]))

        source = self.base_source.clone(stf=self.get_stf(d), **p)
        return source

    def make_dependant(self, xs, pname):
        pass

    def pack_stf(self, stf):
        return [
            stf[p.name] for p in self.problem_parameters_stf[self.stf_type]]

    def pack(self, source):
        x = num.array([
            source.time - self.base_source.time,
            source.north_shift,
            source.east_shift,
            source.depth,
            source.fn,
            source.fe,
            source.fd,
            ] + self.pack_stf(source.stf), dtype=num.float)

        return x

    def random_uniform(self, xbounds, rstate, fixed_magnitude=None):

        x = num.zeros(self.nparameters)
        for i in range(self.nparameters):
            x[i] = rstate.uniform(xbounds[i, 0], xbounds[i, 1])

        return x.tolist()

    def preconstrain(self, x):
        return x

    def get_dependant_bounds(self):
        pass

    @classmethod
    def get_plot_classes(cls):
        from . import plot
        plots = super(SFProblem, cls).get_plot_classes()
        plots.extend([plot.SFLocationPlot])
        return plots


__all__ = '''
    SFProblem
    SFProblemConfig
'''.split()
