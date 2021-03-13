import numpy as num
import logging
import math

from pyrocko import plot
from pyrocko.guts import Tuple, Float

from matplotlib import pyplot as plt, colors as mcolors, cm as mcm
from matplotlib import lines
from matplotlib.ticker import FuncFormatter

from grond.plot.config import PlotConfig

guts_prefix = 'grond'
km = 1e3
d2r = math.pi / 180.

logger = logging.getLogger('targets.plot')


def ndig_inc(vinc):
    assert vinc != 0.0
    ndig = -math.floor(math.log10(abs(vinc)))
    if ndig <= 0:
        return 0
    else:
        return ndig + len(('%5.3f' % (vinc * 10**ndig)).rstrip('0')) - 2


def make_scale(vmin, vmax, approx_ticks=3, mode='0-max', snap=True):
    cscale = plot.AutoScaler(approx_ticks=approx_ticks, mode=mode, snap=snap)
    vmin, vmax, vinc = cscale.make_scale((vmin, vmax))
    imin = math.ceil(vmin/vinc)
    imax = math.floor(vmax/vinc)
    vs = num.arange(imin, imax+1) * vinc
    vexp = cscale.make_exp(vinc)
    vexp_factor = 10**vexp
    vinc_base = vinc / vexp_factor
    ndig = ndig_inc(vinc_base)
    return vs, vexp, '%%.%if' % ndig


def darken(f, c):
    c = num.asarray(c)
    return f*c


class StationDistributionPlot(PlotConfig):
    ''' Plot showing all waveform fits for the ensemble of solutions'''

    name = 'seismic_stations_base'
    size_cm = Tuple.T(
        2, Float.T(),
        default=(16., 13.),
        help='width and length of the figure in cm')
    font_size = Float.T(
        default=10,
        help='font size of all text, except station labels')
    font_size_labels = Float.T(
        default=6,
        help='font size of station labels')

    def plot_station_distribution(
            self, azimuths, distances, weights, labels=None,
            colors=None, cmap=None, cnorm=None, clabel=None,
            scatter_kwargs=dict(), annotate_kwargs=dict(), maxsize=10**2,
            legend_title=''):

        invalid_color = plot.mpl_color('aluminium3')

        scatter_default = {
            'alpha': .5,
            'zorder': 10,
            'c': plot.mpl_color('skyblue2'),
        }

        annotate_default = {
            'alpha': .8,
            'color': 'k',
            'fontsize': self.font_size_labels,
            'ha': 'right',
            'va': 'top',
            'xytext': (-5, -5),
            'textcoords': 'offset points'
        }

        scatter_default.update(scatter_kwargs)
        annotate_default.update(annotate_kwargs)

        fig = plt.figure(figsize=self.size_inch)

        plot.mpl_margins(
            fig, nw=1, nh=1, left=3., right=10., top=3., bottom=3.,
            units=self.font_size)

        ax = fig.add_subplot(111, projection='polar')

        valid = num.isfinite(weights)
        valid[valid] = num.logical_and(valid[valid], weights[valid] > 0.0)

        weights = weights.copy()
        if num.sum(valid) == 0:
            weights[:] = 1.0
            weights_ref = 1.0
        else:
            weights[~valid] = weights[valid].min()
            weights_ref = plot.nice_value(weights[valid].max())

        if weights_ref == 0.:
            weights_ref = 1.0

        if colors is None:
            scolors = [
                scatter_default['c'] if s else invalid_color for s in valid]
        else:
            if cnorm is None:
                cnorm = mcolors.Normalize()
            elif isinstance(cnorm, tuple):
                cnorm = mcolors.Normalize(vmin=cnorm[0], vmax=cnorm[1])

            if cmap is None:
                cmap = plt.get_cmap('viridis')
            elif isinstance(cmap, str):
                cmap = plt.get_cmap(cmap)

            sm = mcm.ScalarMappable(norm=cnorm, cmap=cmap)
            sm.set_array(colors)

            scolors = [
                sm.to_rgba(value) for value in colors]

        scolors = num.array(scolors)

        scatter_default.pop('c')

        weights_scaled = (weights / weights_ref) * maxsize

        ws, exp, fmt = make_scale(0., weights_ref)
        ws = ws[1:]
        weight_clip_min = ws[0]
        weight_clip_min_scaled = (weight_clip_min / weights_ref) * maxsize
        weights_scaled = num.maximum(weight_clip_min_scaled, weights_scaled)

        stations = ax.scatter(
            azimuths*d2r, distances, s=weights_scaled, c=scolors,
            edgecolors=darken(0.5, scolors),
            linewidths=1.0,
            **scatter_default)

        if len(labels) < 30:  # TODO: remove after impl. of collision detection
            if labels is not None:
                for ilbl, label in enumerate(labels):
                    ax.annotate(
                        label, (azimuths[ilbl]*d2r, distances[ilbl]),
                        **annotate_default)

        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.tick_params('y', labelsize=self.font_size, labelcolor='gray')
        ax.grid(alpha=.2)
        ax.set_ylim(0, distances.max()*1.1)
        ax.yaxis.set_major_locator(plt.MaxNLocator(4))
        ax.yaxis.set_major_formatter(
            FuncFormatter(lambda x, pos: '%d km' % (x/km) if x != 0.0 else ''))

        # Legend
        valid_marker = num.argmax(valid)
        ecl = stations.get_edgecolor()
        fc = tuple(stations.get_facecolor()[valid_marker])
        ec = tuple(ecl[min(valid_marker, len(ecl)-1)])

        legend_data = []
        for w in ws[::-1]:
            legend_data.append((
                lines.Line2D(
                    [0], [0],
                    ls='none',
                    marker='o',
                    ms=num.sqrt(w/weights_ref*maxsize),
                    mfc=fc,
                    mec=ec),
                fmt % (w/10**exp) + (
                    '$\\times 10^{%i}$' % exp if exp != 0 else '')))

        if not num.all(valid):
            legend_data.append((
                lines.Line2D(
                    [0], [0],
                    marker='o',
                    ms=num.sqrt(maxsize),
                    mfc=invalid_color,
                    mec=darken(0.5, invalid_color),
                    mew=1.0),
                'Excluded'))

        legend = fig.legend(
            *zip(*legend_data),
            fontsize=self.font_size,
            loc='upper left', bbox_to_anchor=(0.77, 0.4),
            markerscale=1, numpoints=1,
            frameon=False)

        legend.set_title(
            legend_title,
            prop=dict(size=self.font_size))

        cb_axes = fig.add_axes([0.8, 0.6, 0.02, 0.3])

        if clabel is not None:
            fig.colorbar(sm, cax=cb_axes, label=clabel, extend='both')

        return fig, ax, legend
