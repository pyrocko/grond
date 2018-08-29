import numpy as num
from matplotlib import cm, gridspec

from grond.plot.config import PlotConfig
from grond.plot.collection import PlotItem

from matplotlib import pyplot as plt
from pyrocko.guts import Tuple, Float, String, Int

km = 1000.
guts_prefix = 'grond'


def scale_axes(ax, scale):
    from matplotlib.ticker import ScalarFormatter

    class FormatScaled(ScalarFormatter):

        @staticmethod
        def __call__(value, pos):
            return '%d' % (value * scale)

    ax.get_xaxis().set_major_formatter(FormatScaled())
    ax.get_yaxis().set_major_formatter(FormatScaled())


class SatelliteTargetPlot(PlotConfig):
    ''' Maps showing surface displacements from satellite and modelled data '''

    name = 'fits_satellite'
    dpi = Int.T(
        default=250)
    size_cm = Tuple.T(
        2, Float.T(),
        default=(22., 10.))
    colormap = String.T(
        default='RdBu',
        help='Colormap for the surface displacements')

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
            description=u'Maps showing surface displacements'
                        u' from satellite and modelled data.')

    def draw_static_fits(self, ds, history, optimiser):
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

        def decorateAxes(ax, scene, title):
            ax.set_title(title)
            ax.tick_params(length=2)
            if scene.frame.isMeter():
                ax.set_xlabel('East [km]')
                scale_axes(ax, 1. / km)
            elif scene.frame.isDegree():
                ax.set_xlabel('Lon [°]')
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
                    facecolor=(.5, .5, .5), alpha=0.5)

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

        for ifig, (sat_target, result) in enumerate(
                zip(sat_targets, results)):

            scene = sat_target.scene

            fig = plt.figure()
            fig.set_size_inches(*self.size_inch)
            gs = gridspec.GridSpec(
                1, 3,
                wspace=.025, left=.1, bottom=.25,
                right=.95)

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

            axes = [plt.subplot(gs[0, i]) for i in range(3)]
            ax = axes[0]
            ax.imshow(mapDisplacementGrid(stat_obs, scene),
                      extent=im_extent, cmap=self.colormap,
                      vmin=-abs_displ, vmax=abs_displ,
                      origin='lower')
            drawLeaves(ax, scene, offset_e, offset_n)
            drawSource(ax, scene)
            decorateAxes(ax, scene, 'Data')

            if scene.frame.isDegree():
                ax.set_ylabel('Lat [°]')
            elif scene.frame.isMeter():
                ax.set_ylabel('North [km]')

            ax = axes[1]
            ax.imshow(mapDisplacementGrid(stat_syn, scene),
                      extent=im_extent, cmap=self.colormap,
                      vmin=-abs_displ, vmax=abs_displ,
                      origin='lower')
            drawLeaves(ax, scene, offset_e, offset_n)
            drawSource(ax, scene)
            decorateAxes(ax, scene, 'Model')
            ax.get_yaxis().set_visible(False)

            ax = axes[2]
            ax.imshow(mapDisplacementGrid(res, scene),
                      extent=im_extent, cmap=self.colormap,
                      vmin=-abs_displ, vmax=abs_displ,
                      origin='lower')
            drawLeaves(ax, scene, offset_e, offset_n)
            drawSource(ax, scene)
            decorateAxes(ax, scene, 'Residual')
            ax.get_yaxis().set_visible(False)

            for ax in axes:
                ax.set_xlim(llE, urE)
                ax.set_ylim(llN, urN)

            cax = fig.add_axes((.1, .125, .85, .025))
            cax.tick_params(length=2)
            cbar = fig.colorbar(cmw, cax=cax, orientation='horizontal')
            cbar.set_label('LOS Displacement [m]')

            fig.text(.02, .98, 'Scene ID: %s' % scene.meta.scene_id,
                     va='top', alpha=.3)

            yield (item, fig)


def get_plot_classes():
    return [SatelliteTargetPlot]
