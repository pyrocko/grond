import numpy as num
import logging

from pyrocko import gf, util
from pyrocko.guts import String, Float, Dict, Int, Bool

from grond.meta import expand_template, Parameter, has_get_plot_classes

from ..base import Problem, ProblemConfig
from .. import CMTProblem

guts_prefix = 'grond'
logger = logging.getLogger('grond.problems.rectangular.problem')
km = 1e3
as_km = dict(scale_factor=km, scale_unit='km')


class RectangularProblemConfig(ProblemConfig):

    ranges = Dict.T(String.T(), gf.Range.T())
    decimation_factor = Int.T(default=1)
    distance_min = Float.T(default=0.)
    nthreads = Int.T(default=4)
    point_source_target_balancing = Bool.T(
        default=False,
        help='If ``True``, target balancing (if used) is performed on a '
             'moment tensor point source at the events location. It increases '
             'the speed, but might lead to not fully optimized target weights.'
        )

    def get_problem(self, event, target_groups, targets):
        if self.decimation_factor != 1:
            logger.warn(
                'Decimation factor for rectangular source set to %i. Results '
                'may be inaccurate.' % self.decimation_factor)

        base_source = gf.RectangularSource.from_pyrocko_event(
            event,
            anchor='top',
            decimation_factor=self.decimation_factor)

        subs = dict(
            event_name=event.name,
            event_time=util.time_to_str(event.time))

        cmt_problem = None
        if self.point_source_target_balancing:
            base_source_cmt = gf.MTSource.from_pyrocko_event(event)

            stf = gf.HalfSinusoidSTF()
            stf.duration = event.duration or 0.0

            base_source_cmt.stf = stf

            ranges = dict(
                time=self.ranges['time'],
                north_shift=self.ranges['north_shift'],
                east_shift=self.ranges['east_shift'],
                depth=self.ranges['depth'],
                magnitude=gf.Range(
                    start=event.magnitude - 1.,
                    stop=event.magnitude + 1.),
                duration=gf.Range(start=0., stop=stf.duration * 2.),
                rmnn=gf.Range(start=-1.41421, stop=1.41421),
                rmee=gf.Range(start=-1.41421, stop=1.41421),
                rmdd=gf.Range(start=-1.41421, stop=1.41421),
                rmne=gf.Range(start=-1., stop=1.),
                rmnd=gf.Range(start=-1., stop=1.),
                rmed=gf.Range(start=-1., stop=1.))

            cmt_problem = CMTProblem(
                name=expand_template(self.name_template, subs),
                base_source=base_source_cmt,
                distance_min=self.distance_min,
                target_groups=target_groups,
                targets=targets,
                ranges=ranges,
                mt_type='dc',
                stf_type='HalfSinusoidSTF',
                norm_exponent=self.norm_exponent,
                nthreads=self.nthreads)

        problem = RectangularProblem(
            name=expand_template(self.name_template, subs),
            base_source=base_source,
            distance_min=self.distance_min,
            target_groups=target_groups,
            targets=targets,
            ranges=self.ranges,
            norm_exponent=self.norm_exponent,
            nthreads=self.nthreads,
            cmt_problem=cmt_problem)

        return problem


@has_get_plot_classes
class RectangularProblem(Problem):

    problem_parameters = [
        Parameter('east_shift', 'm', label='Easting', **as_km),
        Parameter('north_shift', 'm', label='Northing', **as_km),
        Parameter('depth', 'm', label='Depth', **as_km),
        Parameter('length', 'm', label='Length', optional=False, **as_km),
        Parameter('width', 'm', label='Width', optional=False, **as_km),
        Parameter('slip', 'm', label='Slip', optional=False),
        Parameter('strike', 'deg', label='Strike'),
        Parameter('dip', 'deg', label='Dip'),
        Parameter('rake', 'deg', label='Rake')
    ]

    problem_waveform_parameters = [
        Parameter('nucleation_x', 'offset', label='Nucleation X'),
        Parameter('nucleation_y', 'offset', label='Nucleation Y'),
        Parameter('time', 's', label='Time'),
        Parameter('velocity', 'm/s', label='Rupture Velocity')
    ]

    dependants = []

    distance_min = Float.T(default=0.0)

    cmt_problem = Problem.T(optional=True)

    def set_engine(self, engine):
        self._engine = engine

        if self.cmt_problem is not None:
            self.cmt_problem.set_engine(engine)

    def pack(self, source):
        arr = self.get_parameter_array(source)
        for ip, p in enumerate(self.parameters):
            if p.name == 'time':
                arr[ip] -= self.base_source.time
        return arr

    def get_source(self, x):
        d = self.get_parameter_dict(x)
        p = {}
        for k in self.base_source.keys():
            if k in d:
                p[k] = float(
                    self.ranges[k].make_relative(self.base_source[k], d[k]))

        source = self.base_source.clone(**p)

        return source

    def random_uniform(self, xbounds, rstate, fixed_magnitude=None):
        x = num.zeros(self.nparameters)
        for i in range(self.nparameters):
            x[i] = rstate.uniform(xbounds[i, 0], xbounds[i, 1])

        return x

    def preconstrain(self, x, optimizer=False):
        # source = self.get_source(x)
        # if any(self.distance_min > source.distance_to(t)
        #        for t in self.targets):
        #    raise Forbidden()
        return x

    @classmethod
    def get_plot_classes(cls):
        plots = super(RectangularProblem, cls).get_plot_classes()
        return plots


__all__ = '''
    RectangularProblem
    RectangularProblemConfig
'''.split()
