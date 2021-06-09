import logging

import numpy as num
from matplotlib import cm, patches

from pyrocko import orthodrome as pod
from pyrocko.guts import Float, Bool, Tuple

from pyrocko.plot import mpl_color, mpl_init, automap

from grond.plot.section import SectionPlotConfig, SectionPlot
from grond.plot.collection import PlotItem
from grond.plot.config import PlotConfig
from grond.problems.plot import fixlim
from matplotlib import pyplot as plt

logger = logging.getLogger('grond.problem.singleforce.plot')

guts_prefix = 'grond'

km = 1e3


class SFLocationPlot(SectionPlotConfig):
    ''' MT location plot of the best solutions in three cross-sections. '''
    name = 'location_sf'
    normalisation_gamma = Float.T(
        default=3.,
        help='Normalisation of colors and alpha as :math:`x^\\gamma`.'
             'A linear colormap/alpha with :math:`\\gamma=1`.')

    def make(self, environ):
        environ.setup_modelling()
        cm = environ.get_plot_collection_manager()
        history = environ.get_history(subset='harvest')
        mpl_init(fontsize=self.font_size)
        self._to_be_closed = []
        cm.create_group_mpl(
            self,
            self.draw_figures(history),
            title=u'SF Location',
            section='solution',
            feather_icon='target',
            description=u'''
Location plot of the ensemble of best solutions in three cross-sections.

The coordinate range is defined by the search space given in the config file.
Dots show best single force mechanisms, and colors indicate low (red) and
high (blue) misfit.
''')
        for obj in self._to_be_closed:
            obj.close()

    def draw_figures(self, history, color_p_axis=False):
        from matplotlib import colors

        color = 'black'
        fontsize = self.font_size
        markersize = fontsize * 1.5

        problem = history.problem
        sp = SectionPlot(config=self)
        self._to_be_closed.append(sp)

        fig = sp.fig
        axes_en = sp.axes_xy
        axes_dn = sp.axes_zy
        axes_ed = sp.axes_xz

        bounds = problem.get_combined_bounds()

        models = history.get_sorted_primary_models()[::-1]

        iorder = num.arange(history.nmodels)

        for parname, set_label, set_lim in [
                ['east_shift', sp.set_xlabel, sp.set_xlim],
                ['north_shift', sp.set_ylabel, sp.set_ylim],
                ['depth', sp.set_zlabel, sp.set_zlim]]:

            ipar = problem.name_to_index(parname)
            par = problem.combined[ipar]
            set_label(par.get_label())
            xmin, xmax = fixlim(*par.scaled(bounds[ipar]))
            set_lim(xmin, xmax)

        def scale_size(source):
            return markersize * 1.5

        for axes, xparname, yparname in [
                (axes_en, 'east_shift', 'north_shift'),
                (axes_dn, 'depth', 'north_shift'),
                (axes_ed, 'east_shift', 'depth')]:

            ixpar = problem.name_to_index(xparname)
            iypar = problem.name_to_index(yparname)

            xpar = problem.combined[ixpar]
            ypar = problem.combined[iypar]

            xmin, xmax = fixlim(*xpar.scaled(bounds[ixpar]))
            ymin, ymax = fixlim(*ypar.scaled(bounds[iypar]))

            try:
                axes.set_facecolor(mpl_color('aluminium1'))
            except AttributeError:
                axes.patch.set_facecolor(mpl_color('aluminium1'))

            rect = patches.Rectangle(
                (xmin, ymin), xmax-xmin, ymax-ymin,
                facecolor=mpl_color('white'),
                edgecolor=mpl_color('aluminium2'),
                zorder=1)

            axes.add_patch(rect)

            cmap = cm.ScalarMappable(
                norm=colors.PowerNorm(
                    gamma=self.normalisation_gamma,
                    vmin=iorder.min(),
                    vmax=iorder.max()),

                cmap=plt.get_cmap('coolwarm'))

            for ix, x in enumerate(models):

                source = problem.get_source(x)
                fx = problem.extract(x, ixpar)
                fy = problem.extract(x, iypar)
                sx, sy = xpar.scaled(fx), ypar.scaled(fy)

                # TODO: Add rotation in cross-sections
                color = cmap.to_rgba(iorder[ix])

                alpha = (iorder[ix] - iorder.min()) / \
                    float(iorder.max() - iorder.min())
                alpha = alpha**self.normalisation_gamma

                axes.scatter(
                    [sx], [sy],
                    c=[color],
                    s=[scale_size(source)],
                    alpha=alpha,
                    zorder=2)

        item = PlotItem(name='main')
        return [[item, fig]]


class SFForcePlot(PlotConfig):
    ''' Maps showing horizontal and vertical force
        of the best single force model '''

    name = 'forces_singleforce'

    size_cm = Tuple.T(
        2, Float.T(),
        default=(15., 15.),
        help='width and length of the figure in cm')
    show_topo = Bool.T(
        default=False,
        help='show topography')
    show_grid = Bool.T(
        default=True,
        help='show the lat/lon grid')
    show_rivers = Bool.T(
        default=True,
        help='show rivers on the map')
    radius = Float.T(
        optional=True,
        help='radius of the map around campaign center lat/lon')

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history(subset='harvest')
        optimiser = environ.get_optimiser()
        ds = environ.get_dataset()

        environ.setup_modelling()

        cm.create_group_automap(
            self,
            self.draw_best_sf(ds, history, optimiser),
            title=u'Single Force Source Forces',
            section='solution',
            feather_icon='map',
            description=u'''
Maps showing location and force vectors of the best Single Force Source model.

Arrows show the modelled forces (red arrows). The top plot shows the horizontal
forces and the bottom plot the vertical force. The dot indicates the location
of the best single force source model.
''')

    def draw_best_sf(self, ds, history, optimiser, vertical=False):
        from grond.core import make_stats

        source = history.get_best_source()

        problem = history.problem
        models = history.models

        stats = make_stats(
            problem, models, history.get_primary_chain_misfits())

        def plot_sf(source, stats, ifig, vertical=False):
            orient = 'vertical' if vertical else 'horizontal'

            item = PlotItem(
                name='fig_%i' % ifig,
                attributes={},
                title=u'Best %s single force model force vector' % orient,
                description=u'''
Single force source %s force vector for the best model (red). The circle shows
the 95%% confidence ellipse.
''' % orient)

            event = source.pyrocko_event()

            radius = self.radius
            if radius is None or radius < 30.*km:
                logger.warn(
                    'Radius too small, defaulting to 30 km')
                radius = 30*km

            m = automap.Map(
                width=self.size_cm[0],
                height=self.size_cm[1],
                lat=event.lat,
                lon=event.lon,
                radius=radius,
                show_topo=self.show_topo,
                show_grid=self.show_grid,
                show_rivers=self.show_rivers,
                color_wet=(216, 242, 254),
                color_dry=(238, 236, 230))

            offset_scale = num.abs([source.fn, source.fe, source.fd]).sum()
            size = num.linalg.norm(self.size_cm)

            scale = (size / 5.) / offset_scale

            lat, lon = pod.ne_to_latlon(
                event.lat,
                event.lon,
                source.north_shift,
                source.east_shift)

            stats_dict = stats.get_values_dict()

            if vertical:
                rows = [lon, lat,
                        0., -source.fd * scale,
                        (stats_dict['fn.std'] + stats_dict['fe.std']) * scale,
                        stats_dict['fd.std'] * scale,
                        0.]

            else:
                rows = [lon, lat,
                        source.fe * scale, source.fn * scale,
                        stats_dict['fe.std'] * scale,
                        stats_dict['fn.std'] * scale,
                        0.]

            fontsize = 10.

            default_psxy_style = {
                'h': 0,
                'W': '2.0p,red',
                'A': '+p4p,black+e+a40',
                'G': 'red',
                't': 30,
                'L': True,
                'S': 'e1c/0.95/%d' % fontsize,
            }

            m.gmt.psvelo(
                in_rows=[rows],
                *m.jxyr,
                **default_psxy_style)

            m.gmt.psxy(
                S='c10p',
                in_rows=[[lon, lat]],
                W='1p,black',
                G='orange3',
                *m.jxyr)

            return (item, m)

        ifig = 0
        for vertical in (False, True):
            yield plot_sf(source, stats, ifig, vertical)
            ifig += 1

# TODO Fuzzy Single Force Plot
