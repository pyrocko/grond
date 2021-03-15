import logging

import numpy as num
from matplotlib import cm, patches

from pyrocko.guts import Float

from pyrocko.plot import mpl_color, mpl_init

from grond.plot.section import SectionPlotConfig, SectionPlot
from grond.plot.collection import PlotItem
from grond.problems.plot import fixlim
from matplotlib import pyplot as plt

logger = logging.getLogger('grond.problem.singleforce.plot')

guts_prefix = 'grond'


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
Symbols show best double-couple mechanisms, and colors indicate low (red) and
high (blue) misfit.
''')
        for obj in self._to_be_closed:
            obj.close()

    def draw_figures(self, history, color_p_axis=False):
        from matplotlib import colors

        color = 'black'
        fontsize = self.font_size
        markersize = fontsize * 1.5
        size_small = markersize * 0.5

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
            return size_small

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
                edgecolor=mpl_color('aluminium2'))

            axes.add_patch(rect)

            # fxs = xpar.scaled(problem.extract(models, ixpar))
            # fys = ypar.scaled(problem.extract(models, iypar))

            # axes.set_xlim(*fixlim(num.min(fxs), num.max(fxs)))
            # axes.set_ylim(*fixlim(num.min(fys), num.max(fys)))

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
                    colors=[color],
                    size=[scale_size(source)],
                    alpha=alpha,
                    zorder=1,
                    linewidth=0.25)

        item = PlotItem(name='main')
        return [[item, fig]]
