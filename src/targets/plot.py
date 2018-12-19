import numpy as num
import logging
import math

from pyrocko import plot
from pyrocko.guts import Tuple, Float

from matplotlib import pyplot as plt
from matplotlib import lines
from matplotlib.ticker import FuncFormatter

from grond.plot.config import PlotConfig

guts_prefix = 'grond'
km = 1e3
d2r = math.pi / 180.

logger = logging.getLogger('targets.plot')


class StationDistributionPlot(PlotConfig):
    ''' Plot showing all waveform fits for the ensemble of solutions'''

    name = 'station_distribution'
    size_cm = Tuple.T(
        2, Float.T(),
        default=(14., 13.),
        help='width and length of the figure in cm')
    font_size = Float.T(
        default=6,
        help='Font Size of all fonts, except title')
    font_size_title = Float.T(
        default=10,
        help='Font Size of title')

    def plot_station_distribution(
            self, azimuths, distances, weights, labels=None,
            scatter_kwargs=dict(), annotate_kwargs=dict(), maxsize=10**2):

        invalid_color = plot.mpl_color('scarletred2')

        scatter_default = {
            'alpha': .5,
            'zorder': 10,
            'c': plot.mpl_color('skyblue2'),
        }

        annotate_default = {
            'alpha': .8,
            'color': 'k',
            'fontsize': self.font_size,
            'ha': 'right',
            'va': 'top',
            'xytext': (-5, -5),
            'textcoords': 'offset points'

        }

        scatter_default.update(scatter_kwargs)
        annotate_default.update(annotate_kwargs)

        fig = plt.figure(figsize=self.size_inch)

        plot.mpl_margins(
            fig, nw=1, nh=1, w=5., h=5.,
            units=self.font_size)

        ax = fig.add_subplot(111, projection='polar')

        valid = ~num.isnan(weights)

        weights = weights.copy()
        weights[~valid] = weights[valid].min()
        colors = [scatter_default['c'] if s else invalid_color
                  for s in valid]

        scatter_default.pop('c')

        weights_ref = plot.nice_value(weights[valid].max())

        weights_scaled = (weights / weights_ref) * maxsize

        stations = ax.scatter(
            azimuths*d2r, distances, s=weights_scaled, c=colors,
            **scatter_default)

        annotations = []
        if labels is not None:
            for ilbl, label in enumerate(labels):
                ant = ax.annotate(label, (azimuths[ilbl]*d2r, distances[ilbl]),
                                  **annotate_default)
                annotations.append(ant)

        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.tick_params('y', labelsize=self.font_size, labelcolor='gray')
        ax.grid(alpha=.3)
        ax.set_ylim(0, distances.max()*1.1)
        ax.yaxis.set_major_formatter(
            FuncFormatter(lambda x, pos: '%d km' % (x/km)))

        # Legend
        entries = 4
        valid_marker = num.argmax(valid)
        fc = stations.get_facecolor()[valid_marker]
        ec = stations.get_edgecolor()[valid_marker]

        def get_min_precision(values):
            sig_prec = num.floor(
                num.isfinite(num.log10(weights)))
            return int(abs(sig_prec.min())) + 1

        legend_artists = [
            lines.Line2D(
                [0], [0], ls='none',
                marker='o', ms=num.sqrt(rad), mfc=fc, mec=ec)
            for rad in num.linspace(maxsize, .1*maxsize, entries)
        ]

        sig_prec = get_min_precision(weights)
        legend_annot = [
            '{value:.{prec}f}'.format(value=val, prec=sig_prec)
            for val in num.linspace(weights_ref, .1*weights_ref, entries)
        ]

        if not num.all(valid):
            legend_artists.append(
                lines.Line2D(
                    [0], [0], ls='none',
                    marker='o', ms=num.sqrt(maxsize),
                    mfc=invalid_color, mec=invalid_color))
            legend_annot.append('Excluded')

        legend = fig.legend(
            legend_artists, legend_annot,
            fontsize=self.font_size_title, loc=4,
            markerscale=1, numpoints=1,
            frameon=False)

        return fig, ax, legend
