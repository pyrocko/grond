import numpy as num
from matplotlib import cm, gridspec

from grond.plot.config import PlotConfig
from grond.plot.collection import PlotItem

from matplotlib import pyplot as plt
from matplotlib import patches
from pyrocko.guts import Tuple, Float, String, Int, Bool

km = 1000.
guts_prefix = 'grond'


def scale_axes(axis, scale, offset=0.):
    from matplotlib.ticker import ScalarFormatter

    class FormatScaled(ScalarFormatter):

        @staticmethod
        def __call__(value, pos):
            return '{:,.1f}'.format((offset + value) * scale).replace(',', ' ')

    axis.set_major_formatter(FormatScaled())


class SatelliteTargetDisplacement(PlotConfig):
    ''' Maps showing surface displacements from satellite and modelled data '''

    name = 'satellite'
    dpi = Int.T(
        default=250)
    size_cm = Tuple.T(
        2, Float.T(),
        default=(22., 10.))
    colormap = String.T(
        default='RdBu',
        help='Colormap for the surface displacements')
    relative_coordinates = Bool.T(
        default=False,
        help='Show relative coordinates, initial location centered at 0N, 0E')

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history()
        optimiser = environ.get_optimiser()
        ds = environ.get_dataset()

        environ.setup_modelling()

        cm.create_group_mpl(
            self,
            self.draw_static_fits(ds, history, optimiser),
            title=u'Satellite Surface Displacements',
            section='fits',
            feather_icon='navigation',
            description=u' Maps showing subsampled surface displacements as '   
                        u' observed, modelled and the residual (observed minus'  
                        u' modelled).\n The displacement values predicted by' 
                        u' the orbit-ambiguity ramps are added to the modelled' 
                        u' displacements (middle panels). The color shows the' 
                        u' LOS displacement values associated with, and the' 
                        u' extent of, every quadtree box. The light grey dots' 
                        u' show the focal point of pixels combined in the' 
                        u' quadtree box. This point corresponds to the' 
                        u' position of the modelled data point.\n The large dark' 
                        u' grey dot shows the reference source position. The' 
                        u' grey filled box shows the surface projection of the' 
                        u' modelled source, with the thick-lined edge marking' 
                        u' the upper fault edge. '
                        u' Complete data extent is shown.')

    def draw_static_fits(self, ds, history, optimiser, closeup=False):
        from pyrocko.orthodrome import latlon_to_ne_numpy
        problem = history.problem

        sat_targets = problem.satellite_targets
        for target in sat_targets:
            target.set_dataset(ds)

        gms = problem.combine_misfits(history.misfits)
        isort = num.argsort(gms)
        gms = gms[isort]
        models = history.models[isort, :]
        xbest = models[0, :]

        source = problem.get_source(xbest)
        results = problem.evaluate(xbest, targets=sat_targets)

        def initAxes(ax, scene, title, last_axes=False):
            ax.set_title(title)
            ax.tick_params(length=2)

            if scene.frame.isMeter():
                ax.set_xlabel('Easting [km]')
                scale_x = {'scale': 1./km}
                scale_y = {'scale': 1./km}
                if not self.relative_coordinates:
                    import utm
                    utm_E, utm_N, utm_zone, utm_zone_letter =\
                        utm.from_latlon(source.lat, source.lon)
                    scale_x['offset'] = utm_E
                    scale_y['offset'] = utm_N

                    if last_axes:
                        ax.text(0.975, 0.025,
                                'UTM Zone %d%s' % (utm_zone, utm_zone_letter),
                                va='bottom', ha='right',
                                fontsize=8, alpha=.7,
                                transform=ax.transAxes)

            elif scene.frame.isDegree():
                ax.set_xlabel('Lon [°]')
                scale_x = {'scale': 1.}
                scale_y = {'scale': 1.}
                if not self.relative_coordinates:
                    scale_x['offset'] = source.lon
                    scale_y['offset'] = source.lat

            scale_axes(ax.get_xaxis(), **scale_x)
            scale_axes(ax.get_yaxis(), **scale_y)
            ax.set_aspect('equal')

        def drawSource(ax, scene):
            if scene.frame.isMeter():
                fn, fe = source.outline(cs='xy').T
            elif scene.frame.isDegree():
                fn, fe = source.outline(cs='latlon').T
                fn -= source.lat
                fe -= source.lon

            # source is centered
            ax.scatter(0., 0., color='black', s=3, alpha=.5, marker='o')
            ax.fill(fe, fn,
                    edgecolor=(0., 0., 0.),
                    facecolor=(.5, .5, .5), alpha=0.7)
            ax.plot(fe[0:2], fn[0:2], 'k', linewidth=1.3)

        def mapDisplacementGrid(displacements, scene):
            arr = num.full_like(scene.displacement, fill_value=num.nan)
            qt = scene.quadtree

            for syn_v, l in zip(displacements, qt.leaves):
                arr[l._slice_rows, l._slice_cols] = syn_v

            arr[scene.displacement_mask] = num.nan
            return arr

        def drawLeaves(ax, scene, offset_e=0, offset_n=0):
            rects = scene.quadtree.getMPLRectangles()
            for r in rects:
                r.set_edgecolor((.4, .4, .4))
                r.set_linewidth(.5)
                r.set_facecolor('none')
                r.set_x(r.get_x() - offset_e)
                r.set_y(r.get_y() - offset_n)
            map(ax.add_artist, rects)

            ax.scatter(scene.quadtree.leaf_coordinates[:, 0] - offset_e,
                       scene.quadtree.leaf_coordinates[:, 1] - offset_n,
                       s=.25, c='black', alpha=.1)

        def addArrow(ax, scene):
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
                    source.lat, source.lon))
            if target.scene.frame.isDegree():
                off_n = source.lat - target.scene.frame.llLat
                off_e = source.lon - target.scene.frame.llLon

            turE, turN, tllE, tllN = zip(
                *[(l.gridE.max()-off_e,
                   l.gridN.max()-off_n,
                   l.gridE.min()-off_e,
                   l.gridN.min()-off_n)
                  for l in target.scene.quadtree.leaves])

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
                wspace=.05, hspace=.2,
                left=.1, right=.975, top=.95,
                height_ratios=[12, 1])

            item = PlotItem(
                name='fig_%i' % ifig,
                attributes={'targets': [sat_target.path]},
                title=u'Satellite Surface Displacements - %s'
                      % scene.meta.scene_title,
                description=u'''Surface displacements derived from
satellite data, Scene {meta.scene_title} (id: {meta.scene_id}).
 (Left) the input data, (center) the
modelled data and (right) the model residual.'''.format(meta=scene.meta))

            stat_obs = result.statics_obs
            stat_syn = result.statics_syn['displacement.los']
            res = stat_obs - stat_syn

            if scene.frame.isMeter():
                offset_n, offset_e = map(float, latlon_to_ne_numpy(
                    scene.frame.llLat, scene.frame.llLon,
                    source.lat, source.lon))
            elif scene.frame.isDegree():
                offset_n = source.lat - scene.frame.llLat
                offset_e = source.lon - scene.frame.llLon

            im_extent = (scene.frame.E.min() - offset_e,
                         scene.frame.E.max() - offset_e,
                         scene.frame.N.min() - offset_n,
                         scene.frame.N.max() - offset_n)

            abs_displ = num.abs([stat_obs.min(), stat_obs.max(),
                                 stat_syn.min(), stat_syn.max(),
                                 res.min(), res.max()]).max()

            cmw = cm.ScalarMappable(cmap=self.colormap)
            cmw.set_clim(vmin=-abs_displ, vmax=abs_displ)
            cmw.set_array(stat_obs)

            axes = [fig.add_subplot(gs[0, 0]),
                    fig.add_subplot(gs[0, 1]),
                    fig.add_subplot(gs[0, 2])]

            ax = axes[0]
            ax.imshow(mapDisplacementGrid(stat_obs, scene),
                      extent=im_extent, cmap=self.colormap,
                      vmin=-abs_displ, vmax=abs_displ,
                      origin='lower')
            drawLeaves(ax, scene, offset_e, offset_n)
            drawSource(ax, scene)
            addArrow(ax, scene)
            initAxes(ax, scene, 'Observed')

            ax.text(.025, .025, 'Scene ID: %s' % scene.meta.scene_id,
                    fontsize=8, alpha=.7,
                    va='bottom', transform=ax.transAxes)
            if scene.frame.isDegree():
                ax.set_ylabel('Lat [°]')
            elif scene.frame.isMeter():
                ax.set_ylabel('Northing [km]')

            ax = axes[1]
            ax.imshow(mapDisplacementGrid(stat_syn, scene),
                      extent=im_extent, cmap=self.colormap,
                      vmin=-abs_displ, vmax=abs_displ,
                      origin='lower')
            drawLeaves(ax, scene, offset_e, offset_n)
            drawSource(ax, scene)
            addArrow(ax, scene)
            initAxes(ax, scene, 'Model')
            ax.get_yaxis().set_visible(False)

            ax = axes[2]
            ax.imshow(mapDisplacementGrid(res, scene),
                      extent=im_extent, cmap=self.colormap,
                      vmin=-abs_displ, vmax=abs_displ,
                      origin='lower')
            drawLeaves(ax, scene, offset_e, offset_n)
            drawSource(ax, scene)
            addArrow(ax, scene)
            initAxes(ax, scene, 'Residual', last_axes=True)
            ax.get_yaxis().set_visible(False)

            for ax in axes:
                ax.set_xlim(llE, urE)
                ax.set_ylim(llN, urN)

            if closeup:
                if scene.frame.isMeter():
                    fn, fe = source.outline(cs='xy').T
                elif scene.frame.isDegree():
                    fn, fe = source.outline(cs='latlon').T
                    fn -= source.lat
                    fe -= source.lon

                off_n = (fn[0] + fn[1]) / 2
                off_e = (fe[0] + fe[1]) / 2

                fault_size = 2*num.sqrt(max(abs(fn-off_n))**2
                                        + max(abs(fe-off_e))**2)
                fault_size *= self.map_scale

                for ax in axes:
                    ax.set_xlim(-fault_size/2 + off_e, fault_size/2 + off_e)
                    ax.set_ylim(-fault_size/2 + off_n, fault_size/2 + off_n)

            cax = fig.add_subplot(gs[1, :])
            cbar = fig.colorbar(cmw, cax=cax, orientation='horizontal',
                                use_gridspec=True)
            
            cbar.set_label('LOS Displacement [m]')

            return (item, fig)

        for ifig, (sat_target, result) in enumerate(zip(sat_targets, results)):
            yield generate_plot(sat_target, result, ifig)


class SatelliteTargetDisplacementCloseup(SatelliteTargetDisplacement):
    name = 'satellite_closeup'

    map_scale = Float.T(
        default=2.,
        help='Scale the map surroundings, larger value zooms out.')

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history()
        optimiser = environ.get_optimiser()
        ds = environ.get_dataset()

        environ.setup_modelling()

        cm.create_group_mpl(
            self,
            self.draw_static_fits(ds, history, optimiser, closeup=True),
            title=u'Satellite Displacements (Closeup)',
            section='fits',
            feather_icon='zoom-in',
            description=u' Maps showing subsampled surface displacements as '
                        u' observed, modelled and the residual (observed minus'
                        u' modelled).\n The displacement values predicted by'
                        u' the orbit-ambiguity ramps are added to the modelled'
                        u' displacements (middle panels). The color shows the'
                        u' LOS displacement values associated with, and the'
                        u' extent of, every quadtree box. The light grey dots'
                        u' show the focal point of pixels combined in the'
                        u' quadtree box. This point corresponds to the'
                        u' position of the modelled data point.'
                        u' \n The large dark'
                        u' grey dot shows the reference source position. The'
                        u' grey filled box shows the surface projection of the'
                        u' modelled source, with the thick-lined edge marking'
                        u' the upper fault edge. '
                        u' Map is focused around the fault\'s'
                        u' extent.')


def get_plot_classes():
    return [SatelliteTargetDisplacement, SatelliteTargetDisplacementCloseup]
