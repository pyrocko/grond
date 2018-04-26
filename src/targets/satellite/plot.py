import numpy as num
from matplotlib import cm, gridspec

from grond.plot.config import PlotConfig
from grond.plot.collection import PlotItem

from matplotlib import pyplot as plt
from pyrocko.guts import Tuple, Float

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
    size_cm = Tuple.T(
        2, Float.T(),
        default=(20., 10.))

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history()
        optimiser = environ.get_optimiser()
        ds = environ.get_dataset()

        environ.setup_modelling()

        cm.create_group_mpl(
            self,
            self.draw_static_fits(ds, history, optimiser),
            title='Satellite Surface Displacements',
            section='results.satellite',
            description=self.__class__.__name__)

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

        def decorateAxes(ax, title):
            ax.set_title(title)
            ax.set_xlabel('[km]')
            scale_axes(ax, 1. / km)
            ax.set_aspect('equal')

        def drawRectangularOutline(ax):
            source.regularize()
            fn, fe = source.outline(cs='xy').T
            offset_n, offset_e = latlon_to_ne_numpy(
                sat_target.lats[0], sat_target.lons[0],
                source.lat, source.lon)
            fn += offset_n
            fe += offset_e
            ax.plot(offset_e, offset_n, marker='o')
            ax.plot(fe, fn, marker='o')
            # ax.fill(fe, fn, color=(0.5, 0.5, 0.5), alpha=0.5)
            # ax.plot(fe[:2], fn[:2], linewidth=2., color='black', alpha=0.5)

        def mapDisplacementGrid(displacements, scene):
            qt = scene.quadtree
            array = num.empty_like(scene.displacement)
            array.fill(num.nan)
            for syn_v, l in zip(displacements, qt.leaves):
                array[l._slice_rows, l._slice_cols] = syn_v

            array[scene.displacement_mask] = num.nan
            return array

        def drawTiles(ax, scene):
            rect = scene.quadtree.getMPLRectangles()
            for r in rect:
                r.set_edgecolor((.4, .4, .4))
                r.set_linewidth(.5)
                r.set_facecolor('none')
            map(ax.add_artist, rect)

            ax.scatter(scene.quadtree.leaf_coordinates[:, 0],
                       scene.quadtree.leaf_coordinates[:, 1],
                       s=.25, c='black', alpha=.1)

        for ifig, (sat_target, result) in enumerate(
                zip(sat_targets, results)):
            scene = target.scene

            item = PlotItem(
                name='fig_%i' % ifig,
                attributes={
                    'targets': sat_target.path
                },
                title='Satellite Surface Displacements - %s' % scene.title,
                description='''Surface displacements derived from
satellite data, Scene %s (id: %s). (Left) the input data, (center) the
modelled data and (right) the model residual.''')

            fig = plt.figure()
            fig.set_size_inches(*self.size_inch)
            gs = gridspec.GridSpec(
                1, 3,
                hspace=.0001, left=.06, bottom=.1,
                right=.9)

            axes = [plt.subplot(gs[0, i]) for i in range(3)]

            stat_obs = result.statics_obs
            cmw = cm.ScalarMappable(cmap='coolwarm')
            cmw.set_array(stat_obs)
            cmap = cmw.get_cmap()
            norm = cmw.norm

            stat_syn = result.statics_syn['displacement.los']
            res = (stat_obs - stat_syn)
            im_extent = (scene.frame.E.min(), scene.frame.E.max(),
                         scene.frame.N.min(), scene.frame.N.max())

            ax = axes[0]
            ax.imshow(mapDisplacementGrid(stat_obs, scene),
                      extent=im_extent, cmap=cmap,
                      origin='lower', norm=norm)
            drawTiles(ax, scene)
            drawRectangularOutline(ax)
            decorateAxes(ax, 'Data')
            ax.set_ylabel('[km]')

            ax = axes[1]
            ax.imshow(mapDisplacementGrid(stat_syn, scene),
                      extent=im_extent, cmap=cmap,
                      origin='lower', norm=norm)
            drawTiles(ax, scene)
            drawRectangularOutline(ax)
            decorateAxes(ax, 'Model')
            ax.get_yaxis().set_visible(False)

            ax = axes[2]
            ax.imshow(mapDisplacementGrid(res, scene),
                      extent=im_extent, cmap=cmap,
                      origin='lower', norm=norm)
            drawTiles(ax, scene)
            drawRectangularOutline(ax)
            decorateAxes(ax, 'Residual')
            ax.get_yaxis().set_visible(False)

            pos = ax.get_position()
            cax = fig.add_axes([pos.x1 + .01, pos.y0, 0.02, pos.y1 - pos.y0])
            cbar = fig.colorbar(cmw, cax=cax, orientation='vertical')
            cbar.set_label('[m]')

            yield (item, fig)


def get_plot_classes():
    return [SatelliteTargetPlot]
