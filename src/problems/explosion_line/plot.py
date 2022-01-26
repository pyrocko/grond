import numpy as num

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import cm

from grond.plot.config import PlotConfig
from grond.plot.collection import PlotItem

from pyrocko.plot import mpl_init
from pyrocko.guts import Tuple, Float


km = 1e3

def km_formatter(x, pos):
    return "%.1f" % (x / km)


class ExplosionLineLocationPlot(PlotConfig):
    name = 'explosion_line_location'
    size_cm = Tuple.T(2, Float.T(), default=(20., 20.))

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history(subset='harvest')
        optimiser = environ.get_optimiser()

        mpl_init(fontsize=self.font_size)
        cm.create_group_mpl(
            self,
            self.draw_figures(history, optimiser),
            title=u'Explosion Line Location',
            section='solution',
            feather_icon='navigation',
            description=u'''
Geometry of the explosion line source.
''')

    def draw_figures(self, history, optimiser):
        problem = history.problem
        models = history.get_sorted_primary_models()[::-1]
        gms = history.get_primary_chain_misfits()[::-1]

        isort = num.argsort(gms)[::-1]
        iorder = num.arange(isort.size)

        cmap = cm.get_cmap('jet')

        fig = plt.figure(figsize=self.size_inch)
        ax = fig.add_subplot(1, 1, 1)

        for model, order in zip(models, iorder):
            source = problem.get_source(model)

            color = cmap(order/iorder.size)

            ax.arrow(
                source.east_shift, source.north_shift,
                source.end_east_shift - source.east_shift,
                source.end_north_shift - source.north_shift,
                length_includes_head=True,
                fc=color,
                ec='none',
                alpha=order/iorder.size * .1,
                width=.01*source.length
            )

        ax.set_xlabel('Easting [km]')
        ax.set_ylabel('Northing [km]')
        ax.set_aspect('equal')

        ax.xaxis.set_major_formatter(km_formatter)
        ax.yaxis.set_major_formatter(km_formatter)
        ax.grid(alpha=.3)

        item = PlotItem(name='explosion_line_loc')
        yield [item, fig]
