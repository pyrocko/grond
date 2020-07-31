import numpy as num
import logging

from pyrocko import gf, util
from pyrocko.guts import String, Float, Dict, Int, Bool, StringChoice

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
    decimation_factor = Int.T(
        default=1,
        help='Decimation factor of the discretized sup-faults')
    distance_min = Float.T(default=0.)
    nthreads = Int.T(default=1)
    adaptive_resolution = StringChoice.T(
        choices=('off', 'linear', 'uniqueness'),
        default='off')
    adaptive_start = Int.T(
        default=0)

    pure_shear = Bool.T(
        default=True)
    smooth_rupture = Bool.T(
        default=False)

    def get_problem(self, event, target_groups, targets):
        if self.decimation_factor != 1:
            logger.warn(
                'Decimation factor for dynamic rupture source set to %i.'
                ' Results may be inaccurate.' % self.decimation_factor)

        # TODO: handle magnitude and slip
        base_source = gf.PseudoDynamicRupture.from_pyrocko_event(
            event,
            magnitude=None,
            anchor='top',
            nx=self.ranges['nx'].start,
            ny=self.ranges['ny'].start,
            decimation_factor=self.decimation_factor,
            nthreads=self.nthreads,
            pure_shear=self.pure_shear,
            smooth_rupture=self.smooth_rupture)

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
            nthreads=self.nthreads,
            adaptive_resolution=self.adaptive_resolution,
            adaptive_start=self.adaptive_start)

        return problem


@has_get_plot_classes
class DynamicRuptureProblem(Problem):

    problem_parameters = [
        Parameter('east_shift', 'm', label='Easting', **as_km),
        Parameter('north_shift', 'm', label='Northing', **as_km),
        Parameter('depth', 'm', label='Depth', **as_km),
        Parameter('length', 'm', label='Length', **as_km),
        Parameter('width', 'm', label='Width', **as_km),
        Parameter('slip', 'm', label='Slip', optional=True),
        # Parameter('magnitude', label='Magnitude', optional=True),
        Parameter('strike', 'deg', label='Strike'),
        Parameter('dip', 'deg', label='Dip'),
        Parameter('rake', 'deg', label='Rake'),
        Parameter('gamma', 'vr/vs', label=r'$\gamma$'),
        Parameter('nx', 'patches', label='nx'),
        Parameter('ny', 'patches', label='ny')
    ]

    problem_waveform_parameters = [
        Parameter('nucleation_x', 'offset', label='Nucleation X'),
        Parameter('nucleation_y', 'offset', label='Nucleation Y'),
        Parameter('time', 's', label='Time'),
    ]

    dependants = []

    distance_min = Float.T(default=0.0)
    adaptive_resolution = String.T()
    adaptive_start = Int.T()

    def pack(self, source):
        arr = self.get_parameter_array(source)
        for ip, p in enumerate(self.parameters):
            if p.name == 'time':
                arr[ip] -= self.base_source.time
        return arr

    def get_source(self, x):
        src = self.base_source

        d = self.get_parameter_dict(x)
        p = {k: float(self.ranges[k].make_relative(src[k], d[k]))
             for k in src.keys() if k in d}

        p['nx'] = int(p['nx'])
        p['ny'] = int(p['ny'])

        if p['nx'] != src.nx or p['ny'] != src.ny:
            logger.info('refining patches to %dx%d', p['nx'], p['ny'])

        src.nx = p['nx']
        src.ny = p['ny']

        return src.clone(**p)

    def random_uniform(self, xbounds, rstate, fixed_magnitude=None):
        if fixed_magnitude is not None:
            raise GrondError(
                'Setting fixed magnitude in random model generation not '
                'supported for this type of problem.')

        x = num.zeros(self.nparameters)
        for i in range(self.nparameters):
            x[i] = rstate.uniform(xbounds[i, 0], xbounds[i, 1])

        return x

    def preconstrain(self, x, optimiser=None):
        nx = self.ranges['nx']
        ny = self.ranges['ny']

        if optimiser and self.adaptive_resolution == 'linear' and \
                optimiser.iiter >= self.adaptive_start:

            progress = (optimiser.iiter - self.adaptive_start - 1) / \
                (optimiser.niterations - self.adaptive_start)

            nx = num.floor(
                nx.start + progress * (nx.stop - nx.start + 1))
            ny = num.floor(
                ny.start + progress * (ny.stop - ny.start + 1))

        elif optimiser and self.adaptive_resolution == 'uniqueness' and \
                optimiser.iiter >= self.adaptive_start:
            raise NotImplementedError

        else:
            nx = nx.start
            ny = ny.start

        idx_nx = self.get_parameter_index('nx')
        idx_ny = self.get_parameter_index('ny')
        x[idx_nx] = round(nx)
        x[idx_ny] = round(ny)

        return x

    @classmethod
    def get_plot_classes(cls):
        from . import plot
        plots = super(DynamicRuptureProblem, cls).get_plot_classes()
        plots.extend([
            plot.DynamicRuptureSlipMap,
            plot.DynamicRuptureSTF,
            plot.DynamicRuptureMap])
        return plots


class LightDynamicRuptureProblemConfig(DynamicRuptureProblemConfig):

    nx = Int.T(
        default=2)
    ny = Int.T(
        default=2)

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
            decimation_factor=self.decimation_factor,
            nthreads=self.nthreads,
            pure_shear=self.pure_shear,
            smooth_rupture=self.smooth_rupture)

        subs = dict(
            event_name=event.name,
            event_time=util.time_to_str(event.time))

        problem = LightDynamicRuptureProblem(
            name=expand_template(self.name_template, subs),
            base_source=base_source,
            distance_min=self.distance_min,
            target_groups=target_groups,
            targets=targets,
            ranges=self.ranges,
            norm_exponent=self.norm_exponent,
            nthreads=self.nthreads,
            adaptive_resolution=self.adaptive_resolution,
            adaptive_start=self.adaptive_start)

        return problem


@has_get_plot_classes
class LightDynamicRuptureProblem(DynamicRuptureProblem):

    problem_parameters = [
        Parameter('depth', 'm', label='depth', **as_km),
        Parameter('length', 'm', label='Length', **as_km),
        Parameter('width', 'm', label='Width', **as_km),
        Parameter('gamma', 'vr/vs', label=r'$\gamma$')
    ]

    def get_source(self, x):
        src = self.base_source

        d = self.get_parameter_dict(x)
        p = {k: float(self.ranges[k].make_relative(src[k], d[k]))
             for k in src.keys() if k in d}

        return src.clone(**p)

    def preconstrain(self, x, optimiser=None):
        return x

    @classmethod
    def get_plot_classes(cls):
        from . import plot
        plots = super(LightDynamicRuptureProblem, cls).get_plot_classes()
        plots.extend([
            plot.DynamicRuptureSlipMap,
            plot.DynamicRuptureSTF,
            plot.DynamicRuptureMap])
        return plots


__all__ = '''
    DynamicRuptureProblem
    DynamicRuptureProblemConfig
    LightDynamicRuptureProblem
    LightDynamicRuptureProblemConfig
'''.split()
