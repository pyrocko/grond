import numpy as num
import logging

from scipy.interpolate import RegularGridInterpolator

from pyrocko import gf, util
from pyrocko.guts import String, Float, Dict, Int, Bool, StringChoice
from pyrocko.guts_array import Array
from pyrocko.gf.seismosizer import map_anchor
from pyrocko.gf.tractions import TractionField

from grond.meta import expand_template, Parameter, has_get_plot_classes, \
    GrondError, Path

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
            ranges = dict(
                time=self.ranges['time'],
                north_shift=self.ranges['north_shift'],
                east_shift=self.ranges['east_shift'],
                depth=self.ranges['depth'],
                magnitude=gf.Range(
                    start=event.magnitude - 1.,
                    stop=event.magnitude + 1.),
                duration=gf.Range(start=0., stop=event.duration * 2. or 0.),
                rmnn=gf.Range(start=-1.41421, stop=1.41421),
                rmee=gf.Range(start=-1.41421, stop=1.41421),
                rmdd=gf.Range(start=-1.41421, stop=1.41421),
                rmne=gf.Range(start=-1., stop=1.),
                rmnd=gf.Range(start=-1., stop=1.),
                rmed=gf.Range(start=-1., stop=1.))

            base_source_cmt = gf.MTSource.from_pyrocko_event(event)

            stf = gf.HalfSinusoidSTF()
            stf.duration = event.duration or 0.0

            base_source_cmt.stf = stf

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


# class LightDynamicRuptureProblemConfig(DynamicRuptureProblemConfig):
#     base_source_file = Path.T(
#         optional=True,
#         help='Yaml file with basesource with e.g. precalculated patches')

#     n_tractions_x = Int.T(default=5)
#     n_tractions_y = Int.T(default=2)

#     slip = Float.T(
#         optional=True,
#         help='If given, tractions are scaled to satisfy given slip')

#     def get_problem(self, event, target_groups, targets):
#         if self.decimation_factor != 1:
#             logger.warn(
#                 'Decimation factor for dynamic rupture source set to %i.'
#                 ' Results may be inaccurate.' % self.decimation_factor)

#         # TODO: handle magnitude and slip
#         if self.base_source_file is None:
#             base_source = gf.PseudoDynamicRupture.from_pyrocko_event(
#                 event,
#                 slip=self.slip,
#                 anchor=self.anchor,
#                 decimation_factor=self.decimation_factor,
#                 nthreads=self.nthreads,
#                 pure_shear=self.pure_shear,
#                 smooth_rupture=self.smooth_rupture)
#         else:
#             base_source = gf.PseudoDynamicRupture.load(
#                 filename=self.base_source_file)
#             base_source.slip = self.slip

#         subs = dict(
#             event_name=event.name,
#             event_time=util.time_to_str(event.time))

#         tractions = LightDynamicRuptureTraction(
#             length=base_source.length, width=base_source.width,
#             anchor=base_source.anchor, rake=base_source.rake,
#             n_tractions_x=self.n_tractions_x, n_tractions_y=self.n_tractions_y)

#         base_source.rake = None

#         problem = LightDynamicRuptureProblem(
#             name=expand_template(self.name_template, subs),
#             base_source=base_source,
#             tractions=tractions,
#             n_tractions_x=self.n_tractions_x,
#             n_tractions_y=self.n_tractions_y,
#             magnitude=self.magnitude,
#             distance_min=self.distance_min,
#             target_groups=target_groups,
#             targets=targets,
#             ranges=self.ranges,
#             norm_exponent=self.norm_exponent,
#             nthreads=self.nthreads,
#             adaptive_resolution=self.adaptive_resolution,
#             adaptive_start=self.adaptive_start)

#         return problem


# class LightDynamicRuptureTraction(TractionField):
#     n_tractions_x = Int.T(default=5)
#     n_tractions_y = Int.T(default=2)

#     length = Float.T(default=10000.)
#     width = Float.T(default=5000.)
#     anchor = StringChoice.T(
#         choices=tuple(map_anchor.keys()),
#         default='top')

#     tractions = Array.T(
#         optional=True,
#         dtype=num.float,
#         shape=(None, None))

#     rake = Float.T(default=0.)

#     def get_tractions(self, nx, ny, patches=None):
#         if patches is None:
#             raise ValueError('Patches need to be given')

#         anch_x, anch_y = map_anchor[self.anchor]

#         patch_coords = num.array([
#             (float(p.ix), float(p.iy))
#             for p in patches]).reshape(-1, 2)

#         patch_coords[:, 1] -= self.width / 2.

#         pln = self.length / self.n_tractions_x
#         plw = self.width / self.n_tractions_y

#         coords_x = num.zeros(self.n_tractions_x + 2)
#         coords_x[1:-1] = [
#             pln * i + pln/2. for i in range(self.n_tractions_x)]
#         coords_x[0] = 0
#         coords_x[-1] = self.length
#         coords_x -= self.length / 2.

#         coords_y = num.zeros(self.n_tractions_y + 2)
#         coords_y[1:-1] = [
#             plw * i + plw/2. for i in range(self.n_tractions_y)]
#         coords_y[0] = 0
#         coords_y[-1] = self.width
#         coords_y -= self.width / 2.

#         tractions = num.zeros((coords_x.shape[0], coords_y.shape[0]))
#         tractions[1:-1, 1:-1] = self.tractions
#         tractions[0, 1:-1] = self.tractions[0, :]
#         tractions[-1, 1:-1] = self.tractions[-1, :]
#         tractions[1:-1, 0] = self.tractions[:, 0]
#         tractions[1:-1, -1] = self.tractions[:, -1]

#         tractions[0, 0] = self.tractions[0, 0]
#         tractions[0, -1] = self.tractions[0, -1]
#         tractions[-1, 0] = self.tractions[-1, 0]
#         tractions[-1, -1] = self.tractions[-1, -1]

#         interpolator = RegularGridInterpolator(
#             (coords_x, coords_y),
#             tractions)

#         tractions_norm = interpolator(patch_coords)

#         tractions = num.zeros((tractions_norm.shape[0], 3))
#         tractions[:, 0] = num.cos(self.rake*d2r) * tractions_norm
#         tractions[:, 1] = num.sin(self.rake*d2r) * tractions_norm
#         tractions[:, 2] = 0.

#         return tractions


# @has_get_plot_classes
# class LightDynamicRuptureProblem(DynamicRuptureProblem):

#     n_tractions_x = Int.T(default=5)
#     n_tractions_y = Int.T(default=2)

#     tractions = LightDynamicRuptureTraction.T(optional=True)

#     problem_parameters = [
#         Parameter('t00', 'Pa', label='$t_{00}$'),
#         Parameter('t01', 'Pa', label='$t_{01}$'),
#         Parameter('t02', 'Pa', label='$t_{02}$'),
#         Parameter('t03', 'Pa', label='$t_{03}$'),
#         Parameter('t04', 'Pa', label='$t_{04}$'),
#         Parameter('t10', 'Pa', label='$t_{10}$'),
#         Parameter('t11', 'Pa', label='$t_{11}$'),
#         Parameter('t12', 'Pa', label='$t_{12}$'),
#         Parameter('t13', 'Pa', label='$t_{13}$'),
#         Parameter('t14', 'Pa', label='$t_{14}$')]

#     def get_source(self, x):
#         src = self.base_source

#         d = self.get_parameter_dict(x)
#         p = {k: float(self.ranges[k].make_relative(src[k], d[k]))
#              for k in src.keys() if k in d}
#         p['rake'] = None
#         p['magnitude'] = None
#         p['slip'] = None

#         tractions = num.zeros((self.n_tractions_y, self.n_tractions_x))
#         tractions = num.array(
#             [[float(
#                 self.ranges['t%i%i' % (i, j)].make_relative(
#                     0., d['t%i%i' % (i, j)]))
#                 for i in range(self.n_tractions_y)]
#                 for j in range(self.n_tractions_x)])

#         self.tractions.tractions = tractions

#         return src.clone(tractions=self.tractions, **p)

#     def random_uniform(self, xbounds, rstate, fixed_magnitude=None):
#         if fixed_magnitude is not None:
#             raise GrondError(
#                 'Setting fixed magnitude in random model generation not '
#                 'supported for this type of problem.')

#         x = num.zeros(self.nparameters)
#         for i in range(self.nparameters):
#             x[i] = rstate.uniform(xbounds[i, 0], xbounds[i, 1])

#         return x

#     def preconstrain(self, x, optimiser=None):
#         return x

#     @classmethod
#     def get_plot_classes(cls):
#         from . import plot
#         plots = super(LightDynamicRuptureProblem, cls).get_plot_classes()
#         plots.extend([
#             plot.DynamicRuptureSlipMap,
#             plot.DynamicRuptureSTF,
#             plot.DynamicRuptureMap])
#         return plots


__all__ = '''
    DynamicRuptureProblem
    DynamicRuptureProblemConfig
'''.split()
