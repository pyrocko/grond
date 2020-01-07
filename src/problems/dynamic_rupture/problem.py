import numpy as num
import logging

from pyrocko import gf, util
from pyrocko.gf import tractions as tr
from pyrocko.guts import String, Float, Dict, Int, Bool

from grond.meta import expand_template, Parameter, has_get_plot_classes, \
    GrondError

from ..base import Problem, ProblemConfig

guts_prefix = 'grond'

logger = logging.getLogger('grond.problems.dynamic_rupture.problem')
km = 1e3
d2r = 180./num.pi
as_km = dict(scale_factor=km, scale_unit='km')
as_gpa = dict(scale_factor=1e9, scale_unit='GPa')


class DynamicRuptureProblemConfig(ProblemConfig):

    ranges = Dict.T(String.T(), gf.Range.T())
    decimation_factor = Int.T(default=1)
    distance_min = Float.T(default=0.)
    nthreads = Int.T(default=1)
    nx = Int.T(default=10)
    ny = Int.T(default=10)

    pure_shear = Bool.T(
        default=True)

    tractions = tr.TractionField.T(
        default=tr.TractionComposition(
            components=[
                tr.HomogeneousTractions(),
                tr.RectangularTaper()
            ]),
        help='Traction field the rupture plane is exposed to.')

    def get_problem(self, event, target_groups, targets):
        if self.decimation_factor != 1:
            logger.warn(
                'Decimation factor for dynamic rupture source set to %i.'
                ' Results may be inaccurate.' % self.decimation_factor)

        # TODO: handle magnitude and slip
        base_source = gf.PseudoDynamicRupture.from_pyrocko_event(
            event,
            anchor='top',
            nx=self.nx,
            ny=self.ny,
            magnitude=None,
            decimation_factor=self.decimation_factor,
            nthreads=self.nthreads,
            tractions=self.tractions,
            pure_shear=self.pure_shear)

        subs = dict(
            event_name=event.name,
            event_time=util.time_to_str(event.time))

        problem = DynamicRuptureProblem(
            name=expand_template(self.name_template, subs),
            base_source=base_source,
            distance_min=self.distance_min,
            target_groups=target_groups,
            targets=targets,
            ranges=self.ranges,
            norm_exponent=self.norm_exponent,
            nthreads=self.nthreads,  # TODO: remove?
            )

        return problem


@has_get_plot_classes
class DynamicRuptureProblem(Problem):

    problem_parameters = [
        Parameter('east_shift', 'm', label='Easting', **as_km),
        Parameter('north_shift', 'm', label='Northing', **as_km),
        Parameter('depth', 'm', label='Depth', **as_km),
        Parameter('length', 'm', label='Length', **as_km),
        Parameter('width', 'm', label='Width', **as_km),
        Parameter('slip', 'm', label='Slip'),
        # Parameter('magnitude', label='Magnitude'),
        Parameter('strike', 'deg', label='Strike'),
        Parameter('dip', 'deg', label='Dip'),
        Parameter('rake', 'deg', label='Rake'),
        Parameter('nucleation_x', 'offset', label='Nucleation X'),
        Parameter('nucleation_y', 'offset', label='Nucleation Y'),
    ]

    problem_waveform_parameters = [
        Parameter('time', 's', label='Time'),
    ]

    dependants = []

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
             for k in self.base_source.keys()
             if k in d}

        return self.base_source.clone(**p)

    def random_uniform(self, xbounds, rstate, fixed_magnitude=None):
        if fixed_magnitude is not None:
            raise GrondError(
                'Setting fixed magnitude in random model generation not '
                'supported for this type of problem.')

        x = num.zeros(self.nparameters)
        for i in range(self.nparameters):
            x[i] = rstate.uniform(xbounds[i, 0], xbounds[i, 1])

        return x

    def preconstrain(self, x):
        # source = self.get_source(x)
        # if any(self.distance_min > source.distance_to(t)
        #        for t in self.targets):
        # raise Forbidden()
        return x

    @classmethod
    def get_plot_classes(cls):
        from . import plot
        plots = super(DynamicRuptureProblem, cls).get_plot_classes()
        plots.extend([plot.DynamicRuptureSlipMap])
        return plots


__all__ = '''
    DynamicRuptureProblem
    DynamicRuptureProblemConfig
'''.split()
