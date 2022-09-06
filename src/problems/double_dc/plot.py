import math
import logging

import numpy as num
from matplotlib import pyplot as plt

from pyrocko.guts import Float, Bool, Tuple
from pyrocko import gf
from pyrocko.plot import automap, mpl_init, beachball, mpl_color

from grond import stats
from grond.plot.collection import PlotItem
from grond.plot.config import PlotConfig

logger = logging.getLogger('grond.problem.double_dc.plot')

guts_prefix = 'grond'

km = 1e3


class DoubleDCDecompositionPlot(PlotConfig):
    '''
    Double DC decomposition plot.
    '''
    name = 'dc_decomposition'
    size_cm = Tuple.T(2, Float.T(), default=(15., 5.))

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history(subset='harvest')
        mpl_init(fontsize=self.font_size)
        cm.create_group_mpl(
            self,
            self.draw_figures(history),
            title=u'Double DC Decomposition',
            section='solution',
            feather_icon='sun',
            description=u'''
Double DC decomposition of the best and mean solution into its two contributing
focal mechanisms.
Shown are the ensemble best and the ensemble mean. The symbol size indicates
the relative strength of the mechanisms. The inversion result is consistent
and stable if ensemble mean and ensemble best have similar symbol size and
patterns.
''')

    def draw_figures(self, history):
        fontsize = self.font_size

        fig = plt.figure(figsize=self.size_inch)
        axes = fig.add_subplot(1, 1, 1, aspect=1.0)

        fig.subplots_adjust(left=0., right=1., bottom=0., top=1.)

        problem = history.problem
        models = history.models
        if models.size == 0:
            logger.warn('Empty models vector.')
            return []

        mean_source = stats.get_mean_source(
            problem, history.models)
        best_source = history.get_best_source()
        nlines_max = int(round(self.size_cm[1] / 5. * 4. - 1.0))

        def get_deco(source):
            if isinstance(source, gf.DoubleDCSource):
                return [source] + source.split()

        lines = []
        lines.append(
            ('Ensemble best', get_deco(best_source), mpl_color('aluminium5')))
        lines.append(
            ('Ensemble mean', get_deco(mean_source), mpl_color('aluminium5')))

        mag_max = max(dc.magnitude for (_, line, _) in lines for dc in line)

        for xpos, label in [
                (0., 'Double DC'),
                (2., 'DC 1'),
                (4., 'DC 2')]:

            axes.annotate(
                label,
                xy=(1 + xpos, nlines_max),
                xycoords='data',
                xytext=(0., 0.),
                textcoords='offset points',
                ha='center',
                va='center',
                color='black',
                fontsize=fontsize)

        for i, (label, deco, color_t) in enumerate(lines):
            ypos = nlines_max - i - 1.0
            [ddc, dc1, dc2] = deco
            size0 = ddc.magnitude / mag_max

            axes.annotate(
                label,
                xy=(-2., ypos),
                xycoords='data',
                xytext=(0., 0.),
                textcoords='offset points',
                ha='left',
                va='center',
                color='black',
                fontsize=fontsize)

            for xpos, dc_part, ratio, ops in [
                    (0., ddc, 1., '='),
                    (2., dc1, dc1.magnitude / mag_max, '+'),
                    (4., dc2, dc2.magnitude / mag_max, None)]:

                if ratio > 1e-4:
                    try:
                        beachball.plot_beachball_mpl(
                            dc_part.pyrocko_moment_tensor(), axes,
                            position=(1. + xpos, ypos),
                            size=0.9 * size0 * math.sqrt(ratio),
                            size_units='data',
                            color_t=color_t,
                            linewidth=1.0)
                    except beachball.BeachballError as e:
                        logger.warn(str(e))
                        axes.annotate(
                            'ERROR',
                            xy=(1. + xpos, ypos),
                            ha='center',
                            va='center',
                            color='red',
                            fontsize=fontsize)
                else:
                    axes.annotate(
                        'N/A',
                        xy=(1. + xpos, ypos),
                        ha='center',
                        va='center',
                        color='black',
                        fontsize=fontsize)
                if ops is not None:
                    axes.annotate(
                        ops,
                        xy=(2. + xpos, ypos),
                        ha='center',
                        va='center',
                        color='black',
                        fontsize=fontsize)
        axes.axison = False
        axes.set_xlim(-2.25, 9.75)
        axes.set_ylim(-0.5, nlines_max+0.5)
        item = PlotItem(name='main')
        return [[item, fig]]
