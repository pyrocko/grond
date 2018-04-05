import logging

import numpy as num

from matplotlib import pyplot as plt

from pyrocko import plot

from grond.plot import Plotter
from .base import WaveformMisfitResult

logger = logging.getLogger('targets.waveform.plot')


class WaveformPlotter(Plotter):

    @classmethod
    def draw_check_figures(cls, sources, target, results):
        if not results:
            return []

        t0_mean = num.mean([s.time for s in sources])

        # distances = [
        #    s.distance_to(target) for s in sources]

        # distance_min = num.min(distances)
        # distance_max = num.max(distances)

        yabsmaxs = []

        for result in results:
            if isinstance(result, WaveformMisfitResult):
                yabsmaxs.append(
                    num.max(num.abs(
                        result.filtered_obs.get_ydata())))

        if yabsmaxs:
            yabsmax = max(yabsmaxs) or 1.0
        else:
            yabsmax = None

        fontsize = 10

        fig = None
        ii = 0
        for source, result in zip(sources, results):
            if not isinstance(result, WaveformMisfitResult):
                logger.warn(str(result))
                continue

            if result.tobs_shift != 0.0:
                t0 = result.tsyn_pick
            else:
                t0 = t0_mean

            if fig is None:
                fig = plt.figure(figsize=(4, 3))

                labelpos = plot.mpl_margins(
                    fig, nw=1, nh=1, w=1., h=5.,
                    units=fontsize)

                axes = fig.add_subplot(1, 1, 1)

                labelpos(axes, 2.5, 2.0)
                axes.set_frame_on(False)
                axes.set_ylim(1., 4.)
                axes.get_yaxis().set_visible(False)
                axes.set_title('%s' % target.string_id())
                axes.set_xlabel('Time [s]')

            t = result.filtered_obs.get_xdata()
            ydata = result.filtered_obs.get_ydata() / yabsmax
            axes.plot(
                t-t0, ydata*0.40 + 3.5, color='black', lw=1.0)

            color = plot.mpl_graph_color(ii)

            t = result.filtered_syn.get_xdata()
            ydata = result.filtered_syn.get_ydata()
            ydata = ydata / (num.max(num.abs(ydata)) or 1.0)

            axes.plot(t-t0, ydata*0.47 + 2.5, color=color, alpha=0.5, lw=1.0)

            t = result.processed_syn.get_xdata()
            ydata = result.processed_syn.get_ydata()
            ydata = ydata / (num.max(num.abs(ydata)) or 1.0)

            axes.plot(t-t0, ydata*0.47 + 1.5, color=color, alpha=0.5, lw=1.0)
            if result.tobs_shift != 0.0:
                axes.axvline(
                    result.tsyn_pick - t0,
                    color=(0.7, 0.7, 0.7),
                    zorder=2)

            t = result.processed_syn.get_xdata()
            taper = result.taper

            y = num.ones(t.size) * 0.9
            taper(y, t[0], t[1] - t[0])
            y2 = num.concatenate((y, -y[::-1]))
            t2 = num.concatenate((t, t[::-1]))
            axes.plot(t2-t0, y2 * 0.47 + 3.5, color=color, alpha=0.2, lw=1.0)
            ii += 1

        if fig is None:
            return []

        return [fig]
