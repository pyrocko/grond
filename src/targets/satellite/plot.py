import logging
import numpy as num
from matplotlib import cm, gridspec

from grond.plot.config import PlotConfig
from grond.plot.collection import PlotItem

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, FuncFormatter
from matplotlib import patches
from pyrocko.guts import Tuple, Float, String, Int, Bool, StringChoice

logger = logging.getLogger('grond.targets.satellite.plot')

km = 1e3
d2r = num.pi/180.
guts_prefix = 'grond'


def drape_displacements(
        displacement, shad_data, mappable,
        shad_lim=(.4, .98), contrast=1., mask=None):
    '''Map color data (displacement) on shaded relief.'''

    from scipy.ndimage import convolve as im_conv
    # Light source from somewhere above - psychologically the best choice
    # from upper left
    ramp = num.array([[1, 0], [0, -1.]]) * contrast

    # convolution of two 2D arrays
    shad = im_conv(shad_data*km, ramp.T)
    shad *= -1.

    # if there are strong artifical edges in the data, shades get
    # dominated by them. Cutting off the largest and smallest 2% of
    # shades helps
    percentile2 = num.quantile(shad, 0.02)
    percentile98 = num.quantile(shad, 0.98)
    shad[shad > percentile98] = percentile98
    shad[shad < percentile2] = percentile2

    # normalize shading
    shad -= num.nanmin(shad)
    shad /= num.nanmax(shad)

    if mask is not None:
        shad[mask] = num.nan

    # reduce range to balance gray color
    shad *= shad_lim[1] - shad_lim[0]
    shad += shad_lim[0]

    rgb_map = mappable.to_rgba(displacement)
    rgb_map[num.isnan(displacement)] = 1.
    rgb_map[:, :, :3] *= shad[:, :, num.newaxis]

    return rgb_map


def displ2rad(displ, wavelength):
    return (displ % wavelength) / wavelength * num.pi


def scale_axes(axis, scale, offset=0., suffix=''):
    from matplotlib.ticker import ScalarFormatter

    class FormatScaled(ScalarFormatter):

        @staticmethod
        def __call__(value, pos):
            return '{:,.1f}{:}'.format((offset + value) * scale, suffix)\
                .replace(',', ' ')

    axis.set_major_formatter(FormatScaled())


class SatelliteTargetDisplacement(PlotConfig):
    ''' Maps showing surface displacements from satellite and modelled data '''

    name = 'satellite'
    dpi = Int.T(
        default=250)
    size_cm = Tuple.T(
        2, Float.T(),
        default=(22., 12.))
    colormap = String.T(
        default='RdBu',
        help='Colormap for the surface displacements')
    relative_coordinates = Bool.T(
        default=False,
        help='Show relative coordinates, initial location centered at 0N, 0E')
    fit = StringChoice.T(
        default='best', choices=['best', 'mean'],
        help='Show the \'best\' or \'mean\' fits and source model from the'
             ' ensamble.')

    show_topo = Bool.T(
        default=True,
        help='Drape displacements over the topography.')
    displacement_unit = StringChoice.T(
        default='m',
        choices=['m', 'mm', 'cm', 'rad'],
        help="Show results in 'm', 'cm', 'mm' or 'rad' for radians.")
    show_leaf_centres = Bool.T(
        default=True,
        help='show the center points of Quadtree leaves')
    source_outline_color = String.T(
        default='grey',
        help='Choose color of source outline from named matplotlib Colors')
    common_color_scale = Bool.T(
        default=True,
        help='Results shown with common color scale for all satellite '
             'data sets (based on the data)')
    map_limits = Tuple.T(
        4, Float.T(),
        optional=True,
        help='Overwrite map limits in native coordinates. '
             'Use (xmin, xmax, ymin, ymax)')
    nticks_x = Int.T(
        optional=True,
        help='Number of ticks on the x-axis.')

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history(subset='harvest')
        optimiser = environ.get_optimiser()
        ds = environ.get_dataset()

        environ.setup_modelling()

        cm.create_group_mpl(
            self,
            self.draw_static_fits(ds, history, optimiser),
            title=u'InSAR Displacements',
            section='fits',
            feather_icon='navigation',
            description=u'''
Maps showing subsampled surface displacements as observed, modelled and the
residual (observed minus modelled).

The displacement values predicted by the orbit-ambiguity ramps are added to the
modelled displacements (middle panels). The color shows the LOS displacement
values associated with, and the extent of, every quadtree box. The light grey
dots show the focal point of pixels combined in the quadtree box. This point
corresponds to the position of the modelled data point.

The large dark grey dot shows the reference source position. The grey filled
box shows the surface projection of the modelled source, with the thick-lined
edge marking the upper fault edge. Complete data extent is shown.
''')

    def draw_static_fits(self, ds, history, optimiser, closeup=False):
        from pyrocko.orthodrome import latlon_to_ne_numpy
        problem = history.problem

        sat_targets = problem.satellite_targets
        for target in sat_targets:
            target.set_dataset(ds)

        if self.fit == 'best':
            source = history.get_best_source()
            model = history.get_best_model()
        elif self.fit == 'mean':
            source = history.get_mean_source()
            model = history.get_mean_model()

        results = problem.evaluate(model, targets=sat_targets)

        def init_axes(ax, scene, title, last_axes=False):
            ax.set_title(title, fontsize=self.font_size)
            ax.tick_params(length=2)

            if scene.frame.isMeter():
                import utm
                ax.set_xlabel('Easting [km]', fontsize=self.font_size)
                scale_x = dict(scale=1./km)
                scale_y = dict(scale=1./km)
                utm_E, utm_N, utm_zone, utm_zone_letter =\
                    utm.from_latlon(source.effective_lat,
                                    source.effective_lon)
                scale_x['offset'] = utm_E
                scale_y['offset'] = utm_N

                if last_axes:
                    ax.text(0.975, 0.025,
                            'UTM Zone %d%s' % (utm_zone, utm_zone_letter),
                            va='bottom', ha='right',
                            fontsize=8, alpha=.7,
                            transform=ax.transAxes)
                ax.set_aspect('equal')

            elif scene.frame.isDegree():
                scale_x = dict(scale=1., suffix='°')
                scale_y = dict(scale=1., suffix='°')
                scale_x['offset'] = source.effective_lon
                scale_y['offset'] = source.effective_lat

                ax.set_aspect(1./num.cos(source.effective_lat*d2r))

            if self.relative_coordinates:
                scale_x['offset'] = 0.
                scale_y['offset'] = 0.

            nticks_x = 4 if abs(scene.frame.llLon) >= 100 else 5

            ax.xaxis.set_major_locator(MaxNLocator(self.nticks_x or nticks_x))
            ax.yaxis.set_major_locator(MaxNLocator(5))

            ax.scale_x = scale_x
            ax.scale_y = scale_y

            scale_axes(ax.get_xaxis(), **scale_x)
            scale_axes(ax.get_yaxis(), **scale_y)

        def draw_source(ax, scene):
            if scene.frame.isMeter():
                fn, fe = source.outline(cs='xy').T
                fn -= fn.mean()
                fe -= fe.mean()
            elif scene.frame.isDegree():
                fn, fe = source.outline(cs='latlon').T
                fn -= source.effective_lat
                fe -= source.effective_lon

            # source is centered
            ax.scatter(0., 0., color='black', s=3, alpha=.5, marker='o')
            ax.fill(fe, fn,
                    edgecolor=(0., 0., 0.),
                    facecolor=self.source_outline_color,
                    alpha=0.7)
            ax.plot(fe[0:2], fn[0:2], 'k', linewidth=1.3)

        def get_displacement_rgba(displacements, scene, mappable):
            arr = num.full_like(scene.displacement, fill_value=num.nan)
            qt = scene.quadtree

            for syn_v, leaf in zip(displacements, qt.leaves):
                arr[leaf._slice_rows, leaf._slice_cols] = syn_v

            arr[scene.displacement_mask] = num.nan

            if not self.common_color_scale \
                    and not self.displacement_unit == 'rad':
                abs_displ = num.abs(displacements).max()
                mappable.set_clim(-abs_displ, abs_displ)

            if self.show_topo:
                try:
                    elevation = scene.get_elevation()
                    return drape_displacements(arr, elevation, mappable)
                except Exception as e:
                    logger.warning('could not plot hillshaded topo')
                    logger.exception(e)
            rgb_arr = mappable.to_rgba(arr)
            rgb_arr[num.isnan(arr)] = 1.
            rgb_arr[scene.displacement_mask] = 1.

            return rgb_arr

        def draw_leaves(ax, scene, offset_e=0., offset_n=0.):
            rects = scene.quadtree.getMPLRectangles()
            for r in rects:
                r.set_edgecolor((.4, .4, .4))
                r.set_linewidth(.5)
                r.set_facecolor('none')
                r.set_x(r.get_x() - offset_e)
                r.set_y(r.get_y() - offset_n)
            map(ax.add_artist, rects)

            if self.show_leaf_centres:
                ax.scatter(scene.quadtree.leaf_coordinates[:, 0] - offset_e,
                           scene.quadtree.leaf_coordinates[:, 1] - offset_n,
                           s=.25, c='black', alpha=.1)

        def add_arrow(ax, scene):
            phi = num.nanmean(scene.phi)
            los_dx = num.cos(phi + num.pi) * .0625
            los_dy = num.sin(phi + num.pi) * .0625

            az_dx = num.cos(phi - num.pi/2) * .125
            az_dy = num.sin(phi - num.pi/2) * .125

            anchor_x = .9 if los_dx < 0 else .1
            anchor_y = .85 if los_dx < 0 else .975

            az_arrow = patches.FancyArrow(
                x=anchor_x-az_dx, y=anchor_y-az_dy,
                dx=az_dx, dy=az_dy,
                head_width=.025,
                alpha=.5, fc='k',
                head_starts_at_zero=False,
                length_includes_head=True,
                transform=ax.transAxes)

            los_arrow = patches.FancyArrow(
                x=anchor_x-az_dx/2, y=anchor_y-az_dy/2,
                dx=los_dx, dy=los_dy,
                head_width=.02,
                alpha=.5, fc='k',
                head_starts_at_zero=False,
                length_includes_head=True,
                transform=ax.transAxes)

            ax.add_artist(az_arrow)
            ax.add_artist(los_arrow)

        urE, urN, llE, llN = (0., 0., 0., 0.)
        for target in sat_targets:

            if target.scene.frame.isMeter():
                off_n, off_e = map(float, latlon_to_ne_numpy(
                    target.scene.frame.llLat, target.scene.frame.llLon,
                    source.effective_lat, source.effective_lon))
            if target.scene.frame.isDegree():
                off_n = source.effective_lat - target.scene.frame.llLat
                off_e = source.effective_lon - target.scene.frame.llLon

            turE, turN, tllE, tllN = zip(
                *[(leaf.gridE.max()-off_e,
                   leaf.gridN.max()-off_n,
                   leaf.gridE.min()-off_e,
                   leaf.gridN.min()-off_n)
                  for leaf in target.scene.quadtree.leaves])

            turE, turN = map(max, (turE, turN))
            tllE, tllN = map(min, (tllE, tllN))
            urE, urN = map(max, ((turE, urE), (urN, turN)))
            llE, llN = map(min, ((tllE, llE), (llN, tllN)))

        def generate_plot(sat_target, result, ifig):

            scene = sat_target.scene

            fig = plt.figure()
            fig.set_size_inches(*self.size_inch)
            gs = gridspec.GridSpec(
                2, 3,
                wspace=.15, hspace=.2,
                left=.1, right=.975, top=.95,
                height_ratios=[12, 1])

            item = PlotItem(
                name='fig_%i' % ifig,
                attributes={'targets': [sat_target.path]},
                title=u'Satellite Surface Displacements - %s'
                      % scene.meta.scene_title,
                description=u'''
Surface displacements derived from satellite data.
(Left) the input data, (center) the modelled
data and (right) the model residual.
''')

            stat_obs = result.statics_obs['displacement.los']
            stat_syn = result.statics_syn['displacement.los']
            res = stat_obs - stat_syn

            if scene.frame.isMeter():
                offset_n, offset_e = map(float, latlon_to_ne_numpy(
                    scene.frame.llLat, scene.frame.llLon,
                    source.effective_lat, source.effective_lon))
            elif scene.frame.isDegree():
                offset_n = source.effective_lat - scene.frame.llLat
                offset_e = source.effective_lon - scene.frame.llLon

            im_extent = (
                scene.frame.E.min() - offset_e,
                scene.frame.E.max() - offset_e,
                scene.frame.N.min() - offset_n,
                scene.frame.N.max() - offset_n)

            if self.displacement_unit == 'rad':
                wavelength = scene.meta.wavelength
                if wavelength is None:
                    raise AttributeError(
                        'The satellite\'s wavelength is not set')

                stat_obs = displ2rad(stat_obs, wavelength)
                stat_syn = displ2rad(stat_syn, wavelength)
                res = displ2rad(res, wavelength)

                self.colormap = 'hsv'
                data_range = (0., num.pi)

            else:
                abs_displ = num.abs([stat_obs.min(), stat_obs.max(),
                                     stat_syn.min(), stat_syn.max(),
                                     res.min(), res.max()]).max()
                data_range = (-abs_displ, abs_displ)

            cmw = cm.ScalarMappable(cmap=self.colormap)
            cmw.set_clim(*data_range)
            cmw.set_array(stat_obs)

            axes = [fig.add_subplot(gs[0, 0]),
                    fig.add_subplot(gs[0, 1]),
                    fig.add_subplot(gs[0, 2])]

            ax = axes[0]
            ax.imshow(
                get_displacement_rgba(stat_obs, scene, cmw),
                extent=im_extent, origin='lower')
            draw_leaves(ax, scene, offset_e, offset_n)
            draw_source(ax, scene)
            add_arrow(ax, scene)
            init_axes(ax, scene, 'Observed')

            ax.text(.025, .025, 'Scene ID: %s' % scene.meta.scene_id,
                    fontsize=8, alpha=.7,
                    va='bottom', transform=ax.transAxes)
            if scene.frame.isMeter():
                ax.set_ylabel('Northing [km]', fontsize=self.font_size)

            ax = axes[1]
            ax.imshow(
                get_displacement_rgba(stat_syn, scene, cmw),
                extent=im_extent, origin='lower')
            draw_leaves(ax, scene, offset_e, offset_n)
            draw_source(ax, scene)
            add_arrow(ax, scene)
            init_axes(ax, scene, 'Model')
            ax.get_yaxis().set_visible(False)

            ax = axes[2]
            ax.imshow(
                get_displacement_rgba(res, scene, cmw),
                extent=im_extent, origin='lower')

            draw_leaves(ax, scene, offset_e, offset_n)
            draw_source(ax, scene)
            add_arrow(ax, scene)
            init_axes(ax, scene, 'Residual', last_axes=True)
            ax.get_yaxis().set_visible(False)

            for ax in axes:
                ax.set_xlim(*im_extent[:2])
                ax.set_ylim(*im_extent[2:])

            if closeup:
                if scene.frame.isMeter():
                    fn, fe = source.outline(cs='xy').T
                elif scene.frame.isDegree():
                    fn, fe = source.outline(cs='latlon').T
                    fn -= source.effective_lat
                    fe -= source.effective_lon

                if fn.size > 1:
                    off_n = (fn[0] + fn[1]) / 2
                    off_e = (fe[0] + fe[1]) / 2
                else:
                    off_n = fn[0]
                    off_e = fe[0]

                fault_size = 2*num.sqrt(max(abs(fn-off_n))**2
                                        + max(abs(fe-off_e))**2)
                fault_size *= self.map_scale
                if fault_size == 0.0:
                    extent = (scene.frame.N[-1] + scene.frame.E[-1]) / 2
                    fault_size = extent * .25

                for ax in axes:
                    ax.set_xlim(-fault_size/2 + off_e, fault_size/2 + off_e)
                    ax.set_ylim(-fault_size/2 + off_n, fault_size/2 + off_n)

            if self.map_limits is not None:
                xmin, xmax, ymin, ymax = self.map_limits
                assert xmin < xmax, 'bad map_limits xmin > xmax'
                assert ymin < ymax, 'bad map_limits ymin > ymax'

                for ax in axes:
                    ax.set_xlim(
                        xmin/ax.scale_x['scale'] - ax.scale_x['offset'],
                        xmax/ax.scale_x['scale'] - ax.scale_x['offset'],)
                    ax.set_ylim(
                        ymin/ax.scale_y['scale'] - ax.scale_y['offset'],
                        ymax/ax.scale_y['scale'] - ax.scale_y['offset'])

            if self.displacement_unit == 'm':
                def cfmt(x, p):
                    return '%.2f' % x
            elif self.displacement_unit == 'cm':
                def cfmt(x, p):
                    return '%.1f' % (x * 1e2)
            elif self.displacement_unit == 'mm':
                def cfmt(x, p):
                    return '%.0f' % (x * 1e3)
            elif self.displacement_unit == 'rad':
                def cfmt(x, p):
                    return '%.2f' % x
            else:
                raise AttributeError(
                    'unknown displacement unit %s' % self.displacement_unit)

            cbar_args = dict(
                orientation='horizontal',
                format=FuncFormatter(cfmt),
                use_gridspec=True)
            cbar_label = 'LOS Displacement [%s]' % self.displacement_unit

            if self.common_color_scale:
                cax = fig.add_subplot(gs[1, 1])
                cax.set_aspect(.05)
                cbar = fig.colorbar(cmw, cax=cax, **cbar_args)
                cbar.set_label(cbar_label)
            else:
                for idata, data in enumerate((stat_syn, stat_obs, res)):
                    cax = fig.add_subplot(gs[1, idata])
                    cax.set_aspect(.05)

                    if not self.displacement_unit == 'rad':
                        abs_displ = num.abs(data).max()
                        cmw.set_clim(-abs_displ, abs_displ)

                    cbar = fig.colorbar(cmw, cax=cax, **cbar_args)
                    cbar.set_label(cbar_label)

            return (item, fig)

        for ifig, (sat_target, result) in enumerate(zip(sat_targets, results)):
            yield generate_plot(sat_target, result, ifig)


class SatelliteTargetDisplacementCloseup(SatelliteTargetDisplacement):
    ''' Close-up of satellite surface displacements and modelled data. '''
    name = 'satellite_closeup'

    map_scale = Float.T(
        default=2.,
        help='Scale the map surroundings, larger value zooms out.')

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history(subset='harvest')
        optimiser = environ.get_optimiser()
        ds = environ.get_dataset()

        environ.setup_modelling()

        cm.create_group_mpl(
            self,
            self.draw_static_fits(ds, history, optimiser, closeup=True),
            title=u'InSAR Displacements (Closeup)',
            section='fits',
            feather_icon='zoom-in',
            description=u'''
Maps showing subsampled surface displacements as observed, modelled and the
residual (observed minus modelled).

The displacement values predicted by the orbit-ambiguity ramps are added to the
modelled displacements (middle panels). The color shows the LOS displacement
values associated with, and the extent of, every quadtree box. The light grey
dots show the focal point of pixels combined in the quadtree box. This point
corresponds to the position of the modelled data point.

The large dark grey dot shows the reference source position. The grey filled
box shows the surface projection of the modelled source, with the thick-lined
edge marking the upper fault edge. Map is focused around the fault's extent.
''')


def get_plot_classes():
    return [SatelliteTargetDisplacement, SatelliteTargetDisplacementCloseup]
