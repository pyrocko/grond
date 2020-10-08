from __future__ import print_function
import logging
import numpy as num

from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter

from pyrocko.plot import mpl_init, mpl_margins, mpl_color
from pyrocko.guts import Tuple, Float
from pyrocko import trace

from grond.plot.config import PlotConfig
from grond.plot.collection import PlotItem

logger = logging.getLogger('grond.optimiser.highscore.plot')

guts_prefix = 'grond'


def _pcolormesh_same_dim(ax, x, y, v, **kwargs):
    # x, y, v must have the same dimension
    try:
        return ax.pcolormesh(x, y, v, shading='nearest', **kwargs)
    except TypeError:
        # matplotlib versions < 3.3
        return ax.pcolormesh(x, y, v[:-1, :-1], **kwargs)


class HighScoreOptimiserPlot(object):

    def __init__(
            self, optimiser, problem, history, xpar_name, ypar_name,
            movie_filename):

        self.optimiser = optimiser
        self.problem = problem
        self.chains = optimiser.chains(problem, history)
        self.history = history
        self.xpar_name = xpar_name
        self.ypar_name = ypar_name
        self.fontsize = 10.
        self.movie_filename = movie_filename
        self.show = False
        self.iiter = 0
        self.iiter_last_draw = 0
        self._volatile = []
        self._blocks_complete = set()

    def start(self):
        nfx = 1
        nfy = 1

        problem = self.problem

        ixpar = problem.name_to_index(self.xpar_name)
        iypar = problem.name_to_index(self.ypar_name)

        mpl_init(fontsize=self.fontsize)
        fig = plt.figure(figsize=(9.6, 5.4))
        labelpos = mpl_margins(fig, nw=nfx, nh=nfy, w=7., h=5., wspace=7.,
                               hspace=2., units=self.fontsize)

        xpar = problem.parameters[ixpar]
        ypar = problem.parameters[iypar]

        if xpar.unit == ypar.unit:
            axes = fig.add_subplot(nfy, nfx, 1, aspect=1.0)
        else:
            axes = fig.add_subplot(nfy, nfx, 1)

        labelpos(axes, 2.5, 2.0)

        axes.set_xlabel(xpar.get_label())
        axes.set_ylabel(ypar.get_label())

        axes.get_xaxis().set_major_locator(plt.MaxNLocator(4))
        axes.get_yaxis().set_major_locator(plt.MaxNLocator(4))

        xref = problem.get_reference_model()
        axes.axvline(xpar.scaled(xref[ixpar]), color='black', alpha=0.3)
        axes.axhline(ypar.scaled(xref[iypar]), color='black', alpha=0.3)

        self.fig = fig
        self.problem = problem
        self.xpar = xpar
        self.ypar = ypar
        self.axes = axes
        self.ixpar = ixpar
        self.iypar = iypar
        from matplotlib import colors
        n = self.optimiser.nbootstrap + 1
        hsv = num.vstack((
            num.random.uniform(0., 1., n),
            num.random.uniform(0.5, 0.9, n),
            num.repeat(0.7, n))).T

        self.bcolors = colors.hsv_to_rgb(hsv[num.newaxis, :, :])[0, :, :]
        self.bcolors[0, :] = [0., 0., 0.]

        bounds = self.problem.get_combined_bounds()

        from grond import plot
        self.xlim = plot.fixlim(*xpar.scaled(bounds[ixpar]))
        self.ylim = plot.fixlim(*ypar.scaled(bounds[iypar]))

        self.set_limits()

        from matplotlib.colors import LinearSegmentedColormap

        self.cmap = LinearSegmentedColormap.from_list('probability', [
            (1.0, 1.0, 1.0),
            (0.5, 0.9, 0.6)])

        self.writer = None
        if self.movie_filename:
            from matplotlib.animation import FFMpegWriter

            metadata = dict(title=problem.name, artist='Grond')

            self.writer = FFMpegWriter(
                fps=30,
                metadata=metadata,
                codec='libx264',
                bitrate=200000,
                extra_args=[
                    '-pix_fmt', 'yuv420p',
                    '-profile:v', 'baseline',
                    '-level', '3',
                    '-an'])

            self.writer.setup(self.fig, self.movie_filename, dpi=200)

        if self.show:
            plt.ion()
            plt.show()

    def set_limits(self):
        self.axes.autoscale(False)
        self.axes.set_xlim(*self.xlim)
        self.axes.set_ylim(*self.ylim)

    def draw_frame(self):

        self.chains.goto(self.iiter+1)
        msize = 15.

        for artist in self._volatile:
            artist.remove()

        self._volatile[:] = []

        nblocks = self.iiter // 100 + 1

        models = self.history.models[:self.iiter+1]

        for iblock in range(nblocks):
            if iblock in self._blocks_complete:
                continue

            models_add = self.history.models[
                iblock*100:min((iblock+1)*100, self.iiter+1)]

            fx = self.problem.extract(models_add, self.ixpar)
            fy = self.problem.extract(models_add, self.iypar)
            collection = self.axes.scatter(
                    self.xpar.scaled(fx),
                    self.ypar.scaled(fy),
                    color='black',
                    s=msize * 0.15, alpha=0.2, edgecolors='none')

            if models_add.shape[0] != 100:
                self._volatile.append(collection)
            else:
                self._blocks_complete.add(iblock)

        for ichain in range(self.chains.nchains):

            iiters = self.chains.indices(ichain)
            fx = self.problem.extract(models[iiters, :], self.ixpar)
            fy = self.problem.extract(models[iiters, :], self.iypar)

            nfade = 20
            t1 = num.maximum(0.0, iiters - (models.shape[0] - nfade)) / nfade
            factors = num.sqrt(1.0 - t1) * (1.0 + 15. * t1**2)

            msizes = msize * factors

            paths = self.axes.scatter(
                self.xpar.scaled(fx),
                self.ypar.scaled(fy),
                color=self.bcolors[ichain],
                s=msizes, alpha=0.5, edgecolors='none')

            self._volatile.append(paths)

        _, phase, iiter_phase = self.optimiser.get_sampler_phase(self.iiter)

        np = 1000
        models_prob = num.zeros((np, self.problem.nparameters))
        for ip in range(np):
            models_prob[ip, :] = phase.get_sample(
                self.problem, iiter_phase, self.chains)

        fx = self.problem.extract(models_prob, self.ixpar)
        fy = self.problem.extract(models_prob, self.iypar)

        if False:

            bounds = self.problem.get_combined_bounds()

            nx = 20
            ny = 20
            x_edges = num.linspace(
                bounds[self.ixpar][0], bounds[self.ixpar][1], nx)
            y_edges = num.linspace(
                bounds[self.iypar][0], bounds[self.iypar][1], ny)

            p, _, _ = num.histogram2d(fx, fy, bins=(x_edges, y_edges))
            x, y = num.meshgrid(x_edges, y_edges)

            artist = self.axes.pcolormesh(
                self.xpar.scaled(x),
                self.ypar.scaled(y),
                p, cmap=self.cmap, zorder=-1)

            self._volatile.append(artist)

        else:
            collection = self.axes.scatter(
                    self.xpar.scaled(fx),
                    self.ypar.scaled(fy),
                    color='green',
                    s=msize * 0.15, alpha=0.2, edgecolors='none')

            self._volatile.append(collection)

        if self.writer:
            self.writer.grab_frame()

        artist = self.axes.annotate(
            '%i (%s)' % (self.iiter+1, phase.__class__.__name__),
            xy=(0., 1.),
            xycoords='axes fraction',
            xytext=(self.fontsize/2., -self.fontsize/2.),
            textcoords='offset points',
            ha='left',
            va='top',
            fontsize=self.fontsize,
            fontstyle='normal')

        self._volatile.append(artist)

        if self.show:
            plt.draw()

        self.iiter_last_draw = self.iiter + 1

    def finish(self):
        if self.writer:
            self.writer.finish()

        if self.show:
            plt.show()
            plt.ioff()

    def render(self):
        self.start()

        while self.iiter < self.history.nmodels:
            logger.info('Rendering frame %i/%i.'
                        % (self.iiter+1, self.history.nmodels))
            self.draw_frame()
            self.iiter += 1

        self.finish()


def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return num.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


class HighScoreAcceptancePlot(PlotConfig):
    '''Model acceptance plot '''
    name = 'acceptance'
    size_cm = Tuple.T(2, Float.T(), default=(21., 14.9))

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        cm.create_group_mpl(
            self,
            self.draw_figures(environ),
            title=u'Acceptance',
            section='optimiser',
            description=u'''
Model acceptance and accepted model popularities.

The plots in this section can be used to investigate performance and
characteristics of the optimisation algorithm.
''',
            feather_icon='check')

    def draw_figures(self, environ):
        nwindow = 200
        show_raw_acceptance_rates = False
        optimiser = environ.get_optimiser()
        problem = environ.get_problem()
        history = environ.get_history()
        chains = optimiser.chains(problem, history)
        chains.load()

        acceptance = chains.acceptance_history

        nmodels_rate = history.nmodels - (nwindow - 1)
        if nmodels_rate < 1:
            logger.warning(
                'Cannot create plot acceptance: insufficient number of tested '
                'models.')

            return

        acceptance_rate = num.zeros((history.nchains, nmodels_rate))
        for ichain in range(history.nchains):
            acceptance_rate[ichain, :] = trace.moving_sum(
                acceptance[ichain, :], nwindow, mode='valid') / float(nwindow)

        acceptance_n = num.sum(acceptance, axis=0)

        acceptance_any = num.minimum(acceptance_n, 1)

        acceptance_any_rate = trace.moving_sum(
                acceptance_any, nwindow, mode='valid') / float(nwindow)

        acceptance_p = acceptance_n / float(history.nchains)

        popularity = trace.moving_sum(
            acceptance_p, nwindow, mode='valid') \
            / float(nwindow) / acceptance_any_rate

        mpl_init(fontsize=self.font_size)
        fig = plt.figure(figsize=self.size_inch)
        labelpos = mpl_margins(fig, w=7., h=5., units=self.font_size)

        axes = fig.add_subplot(1, 1, 1)
        labelpos(axes, 2.5, 2.0)

        imodels = num.arange(history.nmodels)

        imodels_rate = imodels[nwindow-1:]

        axes.plot(
            acceptance_n/history.nchains * 100.,
            '.',
            ms=2.0,
            color=mpl_color('skyblue2'),
            label='Popularity of Accepted Models',
            alpha=0.3)

        if show_raw_acceptance_rates:
            for ichain in range(chains.nchains):
                axes.plot(imodels_rate, acceptance_rate[ichain, :]*100.,
                          color=mpl_color('scarletred2'), alpha=0.2)

        axes.plot(
            imodels_rate,
            popularity * 100.,
            color=mpl_color('skyblue2'),
            label='Popularity (moving average)')
        axes.plot(
            imodels_rate,
            acceptance_any_rate*100.,
            color='black',
            label='Acceptance Rate (any chain)')

        axes.legend()

        axes.set_xlabel('Iteration')
        axes.set_ylabel('Acceptance Rate, Model Popularity')

        axes.set_ylim(0., 100.)
        axes.set_xlim(0., history.nmodels - 1)
        axes.grid(alpha=.2)
        axes.yaxis.set_major_formatter(FuncFormatter(lambda v, p: '%d%%' % v))

        iiter = 0
        bgcolors = [mpl_color('aluminium1'), mpl_color('aluminium2')]
        for iphase, phase in enumerate(optimiser.sampler_phases):
            axes.axvspan(
                iiter, iiter+phase.niterations,
                color=bgcolors[iphase % len(bgcolors)])

            iiter += phase.niterations

        yield (
            PlotItem(
                name='acceptance',
                description=u'''
Acceptance rate (black line) within a moving window of %d iterations.

A model is considered accepted, if it is accepted in at least one chain. The
popularity of accepted models is shown as blue dots. Popularity is defined as
the percentage of chains accepting the model (100%% meaning acceptance in all
chains). A moving average of the popularities is shown as blue line (same
averaging interval as for the acceptance rate). Different background colors
represent different sampler phases.
''' % nwindow),
            fig)

        mpl_init(fontsize=self.font_size)
        fig = plt.figure(figsize=self.size_inch)
        labelpos = mpl_margins(fig, w=7., h=5., units=self.font_size)

        axes = fig.add_subplot(1, 1, 1)
        labelpos(axes, 2.5, 2.0)

        nwindow2 = max(1, int(history.nmodels / (self.size_inch[1] * 100)))
        nmodels_rate2 = history.nmodels - (nwindow2 - 1)
        acceptance_rate2 = num.zeros((history.nchains, nmodels_rate2))
        for ichain in range(history.nchains):
            acceptance_rate2[ichain, :] = trace.moving_sum(
                acceptance[ichain, :], nwindow2, mode='valid') \
                / float(nwindow2)

        imodels_rate2 = imodels[nwindow2-1:]

        _pcolormesh_same_dim(
            axes,
            imodels_rate2,
            num.arange(history.nchains),
            num.log(0.01+acceptance_rate2),
            cmap='GnBu')

        if history.sampler_contexts is not None:
            axes.plot(
                imodels,
                history.sampler_contexts[:, 1],
                '.',
                ms=2.0,
                color='black',
                label='Breeding Chain',
                alpha=0.3)

        axes.set_xlabel('Iteration')
        axes.set_ylabel('Bootstrap Chain')
        axes.set_xlim(0, history.nmodels - 1)
        axes.set_ylim(0, history.nchains - 1)

        axes.xaxis.grid(alpha=.4)

        yield (
            PlotItem(
                name='acceptance_img',
                description=u'''
Model acceptance per bootstrap chain averaged over %d models (background color,
low to high acceptance as light to dark colors).

Black dots mark the base chains used when sampling new models (directed sampler
phases only).
''' % nwindow2),
            fig)


__all__ = [
    'HighScoreOptimiserPlot', 'HighScoreAcceptancePlot']
