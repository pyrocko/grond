import logging
import numpy as num
from scipy import signal

from matplotlib import cm, pyplot as plt

from pyrocko.guts import Tuple, Float, Int, StringChoice
from pyrocko.plot import mpl_papersize, mpl_margins, mpl_graph_color, mpl_init

from grond.plot.config import PlotConfig
from grond.plot.collection import PlotItem

from .highscore.plot import *  # noqa


logger = logging.getLogger('grond.problem.plot')


def fixlim(lo, hi):
    if lo == hi:
        return lo - 1.0, hi + 1.0
    else:
        return lo, hi


class SequencePlot(PlotConfig):
    ''' Draws single parameter values as a function of optimization progress'''
    name = 'sequence'
    size_cm = Tuple.T(2, Float.T(), default=(21., 14.9))
    misfit_cutoff = Float.T(optional=True)
    ibootstrap = Int.T(optional=True)
    sort_by = StringChoice.T(
        choices=['iteration', 'misfit'],
        default='iteration')
    subplot_layout = Tuple.T(2, Int.T(), default=(2, 3))
    marker_size = Float.T(default=1.5)

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history()
        mpl_init(fontsize=self.font_size)
        cm.create_group_mpl(self, self.draw_figures(history))

    def draw_figures(self, history):
        misfit_cutoff = self.misfit_cutoff
        sort_by = self.sort_by

        problem = history.problem
        models = history.models

        npar = problem.nparameters
        ndep = problem.ndependants
        fontsize = self.font_size
        nfx, nfy = self.subplot_layout

        imodels = num.arange(history.nmodels)
        bounds = problem.get_combined_bounds()

        xref = problem.get_reference_model()

        gms = problem.combine_misfits(history.misfits)
        gms_softclip = num.where(gms > 1.0, 0.2 * num.log10(gms) + 1.0, gms)

        isort = num.argsort(gms)[::-1]

        if sort_by == 'iteration':
            imodels = imodels[isort]
        elif sort_by == 'misfit':
            imodels = num.arange(imodels.size)
        else:
            assert False

        gms = gms[isort]
        gms_softclip = gms_softclip[isort]
        models = models[isort, :]

        iorder = num.empty_like(isort)
        iorder = num.arange(iorder.size)

        if misfit_cutoff is None:
            ibest = num.ones(gms.size, dtype=num.bool)
        else:
            ibest = gms < misfit_cutoff

        def config_axes(axes, nfx, nfy, impl, iplot, nplots):
            if (impl - 1) % nfx != nfx - 1:
                axes.get_yaxis().tick_left()

            if (impl - 1) >= (nfx * (nfy - 1)) or iplot >= nplots - nfx:
                axes.set_xlabel('Iteration')
                if not (impl - 1) // nfx == 0:
                    axes.get_xaxis().tick_bottom()
            elif (impl - 1) // nfx == 0:
                axes.get_xaxis().tick_top()
                axes.set_xticklabels([])
            else:
                axes.get_xaxis().set_visible(False)

        # nfz = (npar + ndep + 1 - 1) / (nfx*nfy) + 1
        cmap = cm.YlOrRd
        cmap = cm.jet
        msize = self.marker_size
        axes = None
        figs = []
        fig = None
        alpha = 0.5
        for ipar in range(npar):
            impl = ipar % (nfx * nfy) + 1

            if impl == 1:
                fig = plt.figure(figsize=self.size_inch)
                labelpos = mpl_margins(
                    fig, nw=nfx, nh=nfy,
                    w=7., h=5.,
                    wspace=7., hspace=2., units=fontsize)

                item = PlotItem(name='fig_%i' % (len(figs)+1))
                item.attributes['parameters'] = []
                figs.append((item, fig))

            par = problem.parameters[ipar]

            figs[-1][0].attributes['parameters'].append(par.name)

            axes = fig.add_subplot(nfy, nfx, impl)
            labelpos(axes, 2.5, 2.0)

            axes.set_ylabel(par.get_label())
            axes.get_yaxis().set_major_locator(plt.MaxNLocator(4))

            config_axes(axes, nfx, nfy, impl, ipar, npar + ndep + 1)

            axes.set_ylim(*fixlim(*par.scaled(bounds[ipar])))
            axes.set_xlim(0, history.nmodels)

            axes.scatter(
                imodels[ibest], par.scaled(models[ibest, ipar]), s=msize,
                c=iorder[ibest], edgecolors='none', cmap=cmap, alpha=alpha)

            axes.axhline(par.scaled(xref[ipar]), color='black', alpha=0.3)

        for idep in range(ndep):
            # ifz, ify, ifx = num.unravel_index(ipar, (nfz, nfy, nfx))
            impl = (npar + idep) % (nfx * nfy) + 1

            if impl == 1:
                fig = plt.figure(figsize=mpl_papersize('a5', 'landscape'))
                labelpos = mpl_margins(
                    fig, nw=nfx, nh=nfy,
                    w=7., h=5.,
                    wspace=7., hspace=2., units=fontsize)

                item = PlotItem(name='fig_%i' % (len(figs)+1))
                item.attributes['parameters'] = []

                figs.append((item, fig))

            par = problem.dependants[idep]
            figs[-1][0].attributes['parameters'].append(par.name)

            axes = fig.add_subplot(nfy, nfx, impl)
            labelpos(axes, 2.5, 2.0)

            axes.set_ylabel(par.get_label())
            axes.get_yaxis().set_major_locator(plt.MaxNLocator(4))

            config_axes(axes, nfx, nfy, impl, npar + idep, npar + ndep + 1)

            axes.set_ylim(*fixlim(*par.scaled(bounds[npar + idep])))
            axes.set_xlim(0, history.nmodels)

            ys = problem.make_dependant(models[ibest, :], par.name)
            axes.scatter(
                imodels[ibest], par.scaled(ys), s=msize, c=iorder[ibest],
                edgecolors='none', cmap=cmap, alpha=alpha)

            y = problem.make_dependant(xref, par.name)
            axes.axhline(par.scaled(y), color='black', alpha=0.3)

        impl = (npar + ndep) % (nfx * nfy) + 1
        if impl == 1:
            fig = plt.figure(figsize=mpl_papersize('a5', 'landscape'))
            labelpos = mpl_margins(fig, nw=nfx, nh=nfy, w=7., h=5., wspace=7.,
                                   hspace=2., units=fontsize)

            item = PlotItem(name='fig_%i' % (len(figs)+1))
            item.attributes['parameters'] = []

            figs.append(fig)

        axes = fig.add_subplot(nfy, nfx, impl)
        labelpos(axes, 2.5, 2.0)

        config_axes(axes, nfx, nfy, impl, npar + ndep, npar + ndep + 1)

        axes.set_ylim(0., 1.5)
        axes.set_yticks([0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
        axes.set_yticklabels(
            ['0.0', '0.2', '0.4', '0.6', '0.8', '1', '10', '100'])

        axes.scatter(
            imodels[ibest], gms_softclip[ibest], c=iorder[ibest],
            s=msize, edgecolors='none', cmap=cmap, alpha=alpha)

        axes.axhspan(1.0, 1.5, color=(0.8, 0.8, 0.8), alpha=0.2)
        axes.axhline(1.0, color=(0.5, 0.5, 0.5), zorder=2)

        axes.set_xlim(0, history.nmodels)
        axes.set_xlabel('Iteration')

        axes.set_ylabel('Misfit')

        return figs


class ContributionsPlot(PlotConfig):
    '''Relative contribution of single targets to the global misfit'''
    name = 'contributions'
    size_cm = Tuple.T(2, Float.T(), default=(21., 14.9))

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history()
        mpl_init(fontsize=self.font_size)
        cm.create_group_mpl(self, self.draw_figures(history))

    def draw_figures(self, history):

        fontsize = self.font_size

        fig = plt.figure(figsize=self.size_inch)
        labelpos = mpl_margins(fig, nw=2, nh=2, w=7., h=5., wspace=2.,
                               hspace=5., units=fontsize)

        problem = history.problem
        if not problem:
            logger.warn('problem not set')
            return []

        models = history.models

        if models.size == 0:
            logger.warn('empty models vector')
            return []

        imodels = num.arange(history.nmodels)

        gms = problem.combine_misfits(history.misfits)**problem.norm_exponent

        isort = num.argsort(gms)[::-1]

        gms = gms[isort]

        gms_softclip = num.where(gms > 1.0, 0.1 * num.log10(gms) + 1.0, gms)

        gcms = problem.combine_misfits(
            history.misfits, get_contributions=True)

        gcms = gcms[isort, :]

        jsort = num.argsort(gcms[-1, :])[::-1]

        # ncols = 4
        # nrows = ((problem.ntargets + 1) - 1) / ncols + 1

        axes = fig.add_subplot(2, 2, 1)
        labelpos(axes, 2.5, 2.0)

        axes.set_ylabel('Relative contribution (smoothed)')
        axes.set_ylim(0.0, 1.0)

        axes2 = fig.add_subplot(2, 2, 3, sharex=axes)
        labelpos(axes2, 2.5, 2.0)

        axes2.set_xlabel(
            'Tested model, sorted descending by global misfit value')

        axes2.set_ylabel('Square of misfit')

        axes2.set_ylim(0., 1.5)
        axes2.axhspan(1.0, 1.5, color=(0.8, 0.8, 0.8))
        axes2.set_yticks(
            [0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5])
        axes2.set_yticklabels(
            ['0.0', '0.2', '0.4', '0.6', '0.8', '1', '10', '100', '1000',
             '10000', '100000'])

        axes2.set_xlim(imodels[0], imodels[-1])

        rel_ms_sum = num.zeros(history.nmodels)
        rel_ms_smooth_sum = num.zeros(history.nmodels)
        ms_smooth_sum = num.zeros(history.nmodels)
        b = num.hanning(100)
        b /= num.sum(b)
        a = [1]
        ii = 0
        for itarget in jsort:
            target = problem.targets[itarget]
            ms = gcms[:, itarget]
            ms = num.where(num.isfinite(ms), ms, 0.0)
            if num.all(ms == 0.0):
                continue

            rel_ms = ms / gms

            rel_ms_smooth = signal.filtfilt(b, a, rel_ms)

            ms_smooth = rel_ms_smooth * gms_softclip

            rel_poly_y = num.concatenate(
                [rel_ms_smooth_sum[::-1], rel_ms_smooth_sum + rel_ms_smooth])
            poly_x = num.concatenate([imodels[::-1], imodels])

            add_args = {}
            if ii < 20:
                add_args['label'] = '%s (%.2g)' % (
                    target.string_id(), num.mean(rel_ms[-1]))

            axes.fill(
                poly_x, rel_poly_y,
                alpha=0.5,
                color=mpl_graph_color(ii),
                **add_args)

            poly_y = num.concatenate(
                [ms_smooth_sum[::-1], ms_smooth_sum + ms_smooth])

            axes2.fill(poly_x, poly_y, alpha=0.5, color=mpl_graph_color(ii))

            rel_ms_sum += rel_ms

            # axes.plot(
            #    imodels, rel_ms_sum, color='black', alpha=0.1, zorder=-1)

            ms_smooth_sum += ms_smooth
            rel_ms_smooth_sum += rel_ms_smooth
            ii += 1

        axes.legend(
            title='Contributions (top twenty)',
            bbox_to_anchor=(1.05, 0.0, 1.0, 1.0),
            loc='upper left',
            ncol=1, borderaxespad=0., prop={'size': 9})

        axes2.plot(imodels, gms_softclip, color='black')
        axes2.axhline(1.0, color=(0.5, 0.5, 0.5))

        return [[PlotItem(name='main'), fig]]


class BootstrapPlot(PlotConfig):
    ''' Sorted misfit (descending) of single bootstrap chains'''
    name = 'bootstrap'
    size_cm = Tuple.T(2, Float.T(), default=(21., 14.9))

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history()
        optimiser = environ.get_optimiser()
        mpl_init(fontsize=self.font_size)
        cm.create_group_mpl(self, self.draw_figures(history, optimiser))

    def draw_figures(self, history, optimiser):

        fig = plt.figure()

        problem = history.problem
        gms = problem.combine_misfits(history.misfits)

        imodels = num.arange(history.nmodels)

        axes = fig.add_subplot(1, 1, 1)

        gms_softclip = num.where(gms > 1.0, 0.1 * num.log10(gms) + 1.0, gms)

        ibests = []
        for ibootstrap in range(optimiser.nbootstrap):
            bms = optimiser.bootstrap_misfits(
                problem, history.misfits, ibootstrap)

            isort_bms = num.argsort(bms)[::-1]

            ibests.append(isort_bms[-1])

            bms_softclip = num.where(
                bms > 1.0, 0.1 * num.log10(bms) + 1.0, bms)
            axes.plot(imodels, bms_softclip[isort_bms], color='red', alpha=0.2)

        isort = num.argsort(gms)[::-1]
        iorder = num.empty(isort.size)
        iorder[isort] = imodels

        axes.plot(iorder[ibests], gms_softclip[ibests], 'x', color='black')

        m = num.median(gms[ibests])
        s = num.std(gms[ibests])

        axes.axhline(m + s, color='black', alpha=0.5)
        axes.axhline(m, color='black')
        axes.axhline(m - s, color='black', alpha=0.5)

        axes.plot(imodels, gms_softclip[isort], color='black')

        axes.set_xlim(imodels[0], imodels[-1])
        axes.set_xlabel(
            'Tested model, sorted descending by global misfit value')

        return [[PlotItem(name='main'), fig]]


def get_plot_classes():
    return [SequencePlot, ContributionsPlot, BootstrapPlot]
