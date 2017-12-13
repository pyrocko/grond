import math
import numpy as num

from matplotlib import pyplot as plt

from pyrocko.plot import mpl_init, mpl_margins
from grond import plot


class HighScoreOptimizerPlot(object):

    def __init__(self, optimizer, problem, history, xpar_name, ypar_name):
        self.optimizer = optimizer
        self.problem = problem
        self.history = history
        self.xpar_name = xpar_name
        self.ypar_name = ypar_name
        self.fontsize = 10.
        self.movie_filename = None
        self.show = True
        self.iiter = None

    def start(self):
        self.iiter = 0
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

        xref = problem.get_xref()
        axes.axhline(xpar.scaled(xref[ixpar]), color='black', alpha=0.3)
        axes.axvline(ypar.scaled(xref[iypar]), color='black', alpha=0.3)

        self.fig = fig
        self.problem = problem
        self.xpar = xpar
        self.ypar = ypar
        self.axes = axes
        self.ixpar = ixpar
        self.iypar = iypar
        from matplotlib import colors
        n = self.optimizer.nbootstrap + 1
        hsv = num.vstack((
            num.random.uniform(0., 1., n),
            num.random.uniform(0.5, 0.9, n),
            num.repeat(0.7, n))).T

        self.bcolors = colors.hsv_to_rgb(hsv[num.newaxis, :, :])[0, :, :]
        self.bcolors[0, :] = [0., 0., 0.]

        bounds = self.problem.get_combined_bounds()

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
                bitrate=200000)

            self.writer.setup(self.fig, self.movie_filename, dpi=200)

        #if self.show:
            #plt.ion()
            #plt.show()

    def set_limits(self):
        self.axes.set_xlim(*self.xlim)
        self.axes.set_ylim(*self.ylim)

    def draw_frame(self, iiter):
        msize = 15.

        self.axes.cla()

        fx = self.problem.extract(self.history.models[:iiter], self.ixpar)
        fy = self.problem.extract(self.history.models[:iiter], self.iypar)

        self.axes.scatter(
            self.xpar.scaled(fx),
            self.ypar.scaled(fy),
            color='black',
            s=msize * 0.15, alpha=0.2, edgecolors='none')

    def update(self, xhist, chains_i, ibase, jchoice, local_sxs, factor, phase,
               compensate_excentricity):

        msize = 15.

        iiter_frame = len(xhist)

        self.axes.cla()

        if jchoice is not None and local_sxs is not None:

            nx = 100
            ny = 100

            sx = local_sxs[jchoice][self.ixpar] * factor
            sy = local_sxs[jchoice][self.iypar] * factor

            p = num.zeros((ny, nx))

            for j in [jchoice]:  # range(self.problem.nbootstrap+1):

                if compensate_excentricity:
                    ps = core.excentricity_compensated_probabilities(
                            xhist[chains_i[j, :], :], local_sxs[jchoice], 2.)
                else:
                    ps = num.ones(chains_i.shape[1])

                bounds = self.get_combined_bounds()

                x = num.linspace(
                    bounds[self.ixpar][0], bounds[self.ixpar][1], nx)
                y = num.linspace(
                    bounds[self.iypar][0], bounds[self.iypar][1], ny)

                for ichoice in range(chains_i.shape[1]):
                    iiter = chains_i[j, ichoice]
                    vx = xhist[iiter, self.ixpar]
                    vy = xhist[iiter, self.iypar]
                    pdfx = 1.0 / math.sqrt(2.0 * sx**2 * math.pi) * num.exp(
                        -(x - vx)**2 / (2.0 * sx**2))

                    pdfy = 1.0 / math.sqrt(2.0 * sy**2 * math.pi) * num.exp(
                        -(y - vy)**2 / (2.0 * sy**2))

                    p += ps[ichoice] * pdfx[num.newaxis, :] * \
                        pdfy[:, num.newaxis]

            self.axes.pcolormesh(x, y, p, cmap=self.cmap)

        fx = self.problem.extract(xhist, self.ixpar)
        fy = self.problem.extract(xhist, self.iypar)

        self.axes.scatter(
            self.xpar.scaled(fx),
            self.ypar.scaled(fy),
            color='black',
            s=msize * 0.15, alpha=0.2, edgecolors='none')

        for ibootstrap in range(self.optimizer.nbootstrap + 1):

            iiters = chains_i[ibootstrap, :]
            fx = self.problem.extract(xhist[iiters, :], self.ixpar)
            fy = self.problem.extract(xhist[iiters, :], self.iypar)

            nfade = 20
            factors = 1.0 + 5.0 * (num.maximum(
                0.0, iiters - (xhist.shape[0] - nfade)) / nfade)**2

            msizes = msize * factors

            self.axes.scatter(
                self.xpar.scaled(fx),
                self.ypar.scaled(fy),
                color=self.bcolors[ibootstrap],
                s=msizes, alpha=0.5, edgecolors='none')

        # if ibase is not None:
        #     fx = self.problem.extract(
        #         xhist[(ibase, -1), :], self.ixpar)
        #     fy = self.problem.extract(
        #         xhist[(ibase, -1), :], self.iypar)

        #     self.axes.plot(
        #         self.xpar.scaled(fx),
        #         self.ypar.scaled(fy),
        #         color='black')

        fx = self.problem.extract(xhist[-1:, :], self.ixpar)
        fy = self.problem.extract(xhist[-1:, :], self.iypar)

        self.axes.scatter(
            self.xpar.scaled(fx),
            self.ypar.scaled(fy),
            s=msize * 5.0,
            color='none',
            edgecolors=self.bcolors[ibootstrap])

        self.axes.annotate(
            '%i (%s)' % (iiter_frame, phase),
            xy=(0., 1.),
            xycoords='axes fraction',
            xytext=(self.fontsize/2., -self.fontsize/2.),
            textcoords='offset points',
            ha='left',
            va='top',
            fontsize=self.fontsize,
            fontstyle='normal')

        self.set_limits()

        self.post_update()

        if self.writer:
            self.writer.grab_frame()

        if self.show:
            plt.draw()

    def post_update(self):
        pass

    def finish(self):
        if self.writer:
            self.writer.finish()

        if self.show:
            plt.show()
            #plt.ioff()


    def render(self):
        self.start()
        
        self.draw_frame(100)
        self.finish()
