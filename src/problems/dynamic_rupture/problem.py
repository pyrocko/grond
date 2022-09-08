import numpy as num
import logging

from pyrocko import gf, util
from pyrocko.guts import String, Float, Dict, Int, Bool, StringChoice
from pyrocko.guts_array import Array
from pyrocko.gf.seismosizer import map_anchor

from grond.meta import expand_template, Parameter, has_get_plot_classes, \
    GrondError

from ..base import Problem, ProblemConfig
from .. import CMTProblem

guts_prefix = 'grond'

logger = logging.getLogger('grond.problems.dynamic_rupture.problem')
km = 1e3
d2r = num.pi / 180.
as_km = dict(scale_factor=km, scale_unit='km')
as_gpa = dict(scale_factor=1e9, scale_unit='GPa')


class DynamicRuptureProblemConfig(ProblemConfig):

    ranges = Dict.T(String.T(), gf.Range.T())

    anchor = StringChoice.T(
        choices=tuple(map_anchor.keys()),
        default='top')

    decimation_factor = Int.T(
        default=1,
        help='Decimation factor of the discretized sub-faults')
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

    patch_mask = Array.T(
        optional=True,
        shape=(None,),
        dtype=num.bool,
        serialize_as='list')

    point_source_target_balancing = Bool.T(
        default=False,
        help='If ``True``, target balancing (if used) is performed on a '
             'moment tensor point source at the events location. It increases '
             'the speed, but might lead to not fully optimized target weights.'
        )

    def get_problem(self, event, target_groups, targets):
        if self.decimation_factor != 1:
            logger.warn(
                'Decimation factor for dynamic rupture source set to %i.'
                ' Results may be inaccurate.' % self.decimation_factor)

        base_source = gf.PseudoDynamicRupture.from_pyrocko_event(
            event,
            anchor=self.anchor,
            nx=int(self.ranges['nx'].start),
            ny=int(self.ranges['ny'].start),
            decimation_factor=self.decimation_factor,
            nthreads=self.nthreads,
            pure_shear=self.pure_shear,
            smooth_rupture=self.smooth_rupture,
            patch_mask=self.patch_mask,
            stf=None)

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
            adaptive_start=self.adaptive_start,
            cmt_problem=cmt_problem)

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
        x[idx_nx] = int(round(nx))
        x[idx_ny] = int(round(ny))

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


class DoublePDRProblemConfig(ProblemConfig):

    ranges = Dict.T(String.T(), gf.Range.T())

    decimation_factor = Int.T(
        default=1,
        help='Decimation factor of the discretized sub-faults')
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

    point_source_target_balancing = Bool.T(
        default=False,
        help='If ``True``, target balancing (if used) is performed on a '
             'moment tensor point source at the events location. It increases '
             'the speed, but might lead to not fully optimized target weights.'
        )

    def get_problem(self, event, target_groups, targets):
        if self.decimation_factor != 1:
            logger.warn(
                'Decimation factor for dynamic rupture source set to %i.'
                ' Results may be inaccurate.' % self.decimation_factor)

        base_source = gf.DoublePDR.from_pyrocko_event(
            event,
            nx1=int(self.ranges['nx1'].start),
            ny1=int(self.ranges['ny1'].start),
            nx2=int(self.ranges['nx2'].start),
            ny2=int(self.ranges['ny2'].start),
            decimation_factor=self.decimation_factor,
            nthreads=self.nthreads,
            pure_shear=self.pure_shear,
            smooth_rupture=self.smooth_rupture,
            stf=None)

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

        problem = DoublePDRProblem(
            name=expand_template(self.name_template, subs),
            base_source=base_source,
            distance_min=self.distance_min,
            target_groups=target_groups,
            targets=targets,
            ranges=self.ranges,
            norm_exponent=self.norm_exponent,
            nthreads=self.nthreads,
            adaptive_resolution=self.adaptive_resolution,
            adaptive_start=self.adaptive_start,
            cmt_problem=cmt_problem)

        return problem


@has_get_plot_classes
class DoublePDRProblem(Problem):

    problem_parameters = [
        Parameter('east_shift', 'm', label='Easting', **as_km),
        Parameter('north_shift', 'm', label='Northing', **as_km),
        Parameter('depth', 'm', label='Depth', **as_km),
        Parameter('length1', 'm', label='Length 1', **as_km),
        Parameter('width1', 'm', label='Width 1', **as_km),
        Parameter('slip1', 'm', label='Slip 1', optional=True),
        Parameter('strike1', 'deg', label='Strike 1'),
        Parameter('dip1', 'deg', label='Dip 1'),
        Parameter('rake1', 'deg', label='Rake 1'),
        Parameter('gamma1', 'vr/vs', label=r'$\gamma$ 1'),
        Parameter('nx1', 'patches', label='nx 1'),
        Parameter('ny1', 'patches', label='ny 1'),
        Parameter('length2', 'm', label='Length 2', **as_km),
        Parameter('width2', 'm', label='Width 2', **as_km),
        Parameter('slip2', 'm', label='Slip 2', optional=True),
        Parameter('strike2', 'deg', label='Strike 2'),
        Parameter('dip2', 'deg', label='Dip 2'),
        Parameter('rake2', 'deg', label='Rake 2'),
        Parameter('gamma2', 'vr/vs', label=r'$\gamma$ 2'),
        Parameter('nx2', 'patches', label='nx 2'),
        Parameter('ny2', 'patches', label='ny 2'),
        Parameter('delta_time', 's', label='$\\Delta$ Time'),
        Parameter('delta_depth', 'm', label='$\\Delta$ Depth'),
        Parameter('azimuth', 'deg', label='Azimuth'),
        Parameter('distance', 'm', label='Distance'),
        Parameter('mix', label='Mix')
    ]

    problem_waveform_parameters = [
        Parameter('nucleation_x1', 'offset', label='Nucleation X 1'),
        Parameter('nucleation_y1', 'offset', label='Nucleation Y 1'),
        Parameter('nucleation_x2', 'offset', label='Nucleation X 2'),
        Parameter('nucleation_y2', 'offset', label='Nucleation Y 2'),
        Parameter('time', 's', label='Time'),
    ]

    dependants = []

    distance_min = Float.T(default=0.0)
    adaptive_resolution = String.T()
    adaptive_start = Int.T()

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
        src = self.base_source

        d = self.get_parameter_dict(x)
        p = {k: float(self.ranges[k].make_relative(src[k], d[k]))
             for k in src.keys() if k in d}

        p['nx1'] = int(p['nx1'])
        p['ny1'] = int(p['ny1'])

        if p['nx1'] != src.nx1 or p['ny1'] != src.ny1:
            logger.info('refining patches to %dx%d', p['nx1'], p['ny1'])

        src.nx1 = p['nx1']
        src.ny1 = p['ny1']

        p['nx2'] = int(p['nx2'])
        p['ny2'] = int(p['ny2'])

        if p['nx2'] != src.nx2 or p['ny2'] != src.ny2:
            logger.info('refining patches to %dx%d', p['nx2'], p['ny2'])

        src.nx2 = p['nx2']
        src.ny2 = p['ny2']

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
        nx1 = self.ranges['nx1']
        ny1 = self.ranges['ny1']
        nx2 = self.ranges['nx2']
        ny2 = self.ranges['ny2']

        if optimiser and self.adaptive_resolution == 'linear' and \
                optimiser.iiter >= self.adaptive_start:

            progress = (optimiser.iiter - self.adaptive_start - 1) / \
                (optimiser.niterations - self.adaptive_start)

            nx1 = num.floor(
                nx1.start + progress * (nx1.stop - nx1.start + 1))
            ny1 = num.floor(
                ny1.start + progress * (ny1.stop - ny1.start + 1))
            nx2 = num.floor(
                nx2.start + progress * (nx2.stop - nx2.start + 1))
            ny2 = num.floor(
                ny2.start + progress * (ny2.stop - ny2.start + 1))

        elif optimiser and self.adaptive_resolution == 'uniqueness' and \
                optimiser.iiter >= self.adaptive_start:
            raise NotImplementedError

        else:
            nx1 = nx1.start
            ny1 = ny1.start
            nx2 = nx2.start
            ny2 = ny2.start

        idx_nx1 = self.get_parameter_index('nx1')
        idx_ny1 = self.get_parameter_index('ny1')
        idx_nx2 = self.get_parameter_index('nx2')
        idx_ny2 = self.get_parameter_index('ny2')
        x[idx_nx1] = int(round(nx1))
        x[idx_ny1] = int(round(ny1))
        x[idx_nx2] = int(round(nx2))
        x[idx_ny2] = int(round(ny2))

        return x

    @classmethod
    def get_plot_classes(cls):
        plots = super(DoublePDRProblem, cls).get_plot_classes()

        return plots


__all__ = '''
    DynamicRuptureProblem
    DynamicRuptureProblemConfig
    DoublePDRProblem
    DoublePDRProblemConfig
'''.split()
