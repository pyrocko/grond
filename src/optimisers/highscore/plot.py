from __future__ import print_function
import logging
import numpy as num

from matplotlib import pyplot as plt

from pyrocko.plot import mpl_init, mpl_margins

logger = logging.getLogger('grond.optimiser.highscore.plot')


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

        phase, iiter_phase = self.optimiser.get_sampler_phase(self.iiter)

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
            logger.info('rendering frame %i/%i' % (
                self.iiter+1, self.history.nmodels))
            self.draw_frame()
            self.iiter += 1

        self.finish()


__all__ = [
    'HighScoreOptimiserPlot']
