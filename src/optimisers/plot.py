import logging
import numpy as num
from scipy import signal

from matplotlib import cm, pyplot as plt

from pyrocko.guts import Tuple, Float, Int, StringChoice, Bool
from pyrocko.plot import mpl_margins, mpl_graph_color, mpl_init

from grond.plot.config import PlotConfig
from grond.plot.collection import PlotItem

logger = logging.getLogger('grond.problem.plot')

guts_prefix = 'grond'


def fixlim(lo, hi):
    if lo == hi:
        return lo - 1.0, hi + 1.0
    else:
        return lo, hi


class SequencePlot(PlotConfig):
    '''
    Draws all parameter values evaluated during the optimisation

    The sequence of all the parameter values is either a function of the
    optimisation in progress or of the misfit from high to low. This plot can
    be used to check on convergence or see if model parameters push the given
    bounds. The color always shows the relative misfit. Relatively high misfits
    are in cold blue colors and relatively low misfits in red. The last panel
    gives the corresponding misfit values.
    '''

    name = 'sequence'
    size_cm = Tuple.T(2, Float.T(), default=(14., 6.))
    misfit_cutoff = Float.T(optional=True)
    ibootstrap = Int.T(optional=True)
    sort_by = StringChoice.T(
        choices=['iteration', 'misfit'],
        default='iteration')
    subplot_layout = Tuple.T(2, Int.T(), default=(1, 1))
    marker_size = Float.T(default=1.5)
    show_reference = Bool.T(default=True)

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history()
        optimiser = environ.get_optimiser()

        mpl_init(fontsize=self.font_size)
        cm.create_group_mpl(
            self,
            self.draw_figures(history, optimiser),
            title=u'Sequence Plots',
            section='optimiser',
            description=u'''
Sequence plots for all parameters of the optimisation.

The sequence of all the parameter values is either a function of the
optimisation in progress or of the misfit from high to low. This plot can be
used to check on convergence or to see if model parameters push the given
bounds.

The color always shows the relative misfit. Relatively high misfits are in
cold blue colors and relatively low misfits in red. The last panel gives the
corresponding misfit values.
''',
            feather_icon='fast-forward')

    def draw_figures(self, history, optimiser):
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

        gms = history.get_primary_chain_misfits()
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
        fig = None
        item_fig = None
        nfigs = 0
        alpha = 0.5
        for ipar in range(npar):
            impl = ipar % (nfx * nfy) + 1

            if impl == 1:
                if item_fig:
                    yield item_fig
                    nfigs += 1

                fig = plt.figure(figsize=self.size_inch)
                labelpos = mpl_margins(
                    fig, nw=nfx, nh=nfy,
                    left=7.,
                    right=2.,
                    top=1.,
                    bottom=5.,
                    wspace=7., hspace=2., units=fontsize)

                item = PlotItem(name='fig_%i' % (nfigs+1))
                item.attributes['parameters'] = []
                item_fig = (item, fig)

            par = problem.parameters[ipar]

            item_fig[0].attributes['parameters'].append(par.name)

            axes = fig.add_subplot(nfy, nfx, impl)
            labelpos(axes, 2.5, 2.0)

            axes.set_ylabel(par.get_label())
            axes.get_yaxis().set_major_locator(plt.MaxNLocator(4))

            config_axes(axes, nfx, nfy, impl, ipar, npar + ndep + 1)

            axes.set_ylim(*fixlim(*par.scaled(bounds[ipar])))
            axes.set_xlim(0, history.nmodels)

            axes.scatter(
                imodels[ibest], par.scaled(models[ibest, ipar]), s=msize,
                c=iorder[ibest], edgecolors='none', cmap=cmap, alpha=alpha,
                rasterized=True)

            if self.show_reference:
                axes.axhline(par.scaled(xref[ipar]), color='black', alpha=0.3)

        for idep in range(ndep):
            # ifz, ify, ifx = num.unravel_index(ipar, (nfz, nfy, nfx))
            impl = (npar + idep) % (nfx * nfy) + 1

            if impl == 1:
                if item_fig:
                    yield item_fig
                    nfigs += 1

                fig = plt.figure(figsize=self.size_inch)
                labelpos = mpl_margins(
                    fig, nw=nfx, nh=nfy,
                    left=7.,
                    right=2.,
                    top=1.,
                    bottom=5.,
                    wspace=7., hspace=2., units=fontsize)

                item = PlotItem(name='fig_%i' % (nfigs+1))
                item.attributes['parameters'] = []

                item_fig = (item, fig)

            par = problem.dependants[idep]
            item_fig[0].attributes['parameters'].append(par.name)

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
                edgecolors='none', cmap=cmap, alpha=alpha, rasterized=True)

            if self.show_reference:
                y = problem.make_dependant(xref, par.name)
                axes.axhline(par.scaled(y), color='black', alpha=0.3)

        impl = (npar + ndep) % (nfx * nfy) + 1
        if impl == 1:
            if item_fig:
                yield item_fig
                nfigs += 1

            fig = plt.figure(figsize=self.size_inch)
            labelpos = mpl_margins(
                fig, nw=nfx, nh=nfy,
                left=7.,
                right=2.,
                top=1.,
                bottom=5.,
                wspace=7., hspace=2., units=fontsize)

            item = PlotItem(name='fig_%i' % (nfigs+1))
            item.attributes['parameters'] = []

            item_fig = (item, fig)

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

        yield item_fig
        nfigs += 1


class ContributionsPlot(PlotConfig):
    ''' Relative contribution of single targets to the global misfit
    '''

    name = 'contributions'
    size_cm = Tuple.T(2, Float.T(), default=(21., 14.9))

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history()
        optimiser = environ.get_optimiser()
        dataset = environ.get_dataset()

        environ.setup_modelling()

        mpl_init(fontsize=self.font_size)
        cm.create_group_mpl(
            self,
            self.draw_figures(dataset, history, optimiser),
            title=u'Target Contributions',
            section='solution',
            feather_icon='thermometer',
            description=u'''
Contributions of the targets to the total misfit.

The relative contribution that each single target has in the global misfit
result is plotted relative and unscales as a function of global misfit
(descending).

The target contribution is shown in color-filled curves with the bottom curve
on the bottom and the best-fit target on top. This plot can be used to analyse
the balance of targets in the optimisations. For ideal configurations, the
target contributions are of similar size. If the contribution of a single
target is much larger than those of all others, the weighting should be
modified.
''')

    def draw_figures(self, dataset, history, optimiser):

        fontsize = self.font_size

        fig = plt.figure(figsize=self.size_inch)
        labelpos = mpl_margins(fig, nw=2, nh=2, w=7., h=5., wspace=2.,
                               hspace=5., units=fontsize)

        problem = history.problem
        if not problem:
            logger.warn('Problem not set.')
            return []

        models = history.models

        if models.size == 0:
            logger.warn('Empty models vector.')
            return []

        for target in problem.targets:
            target.set_dataset(dataset)

        imodels = num.arange(history.nmodels)

        gms = history.get_sorted_primary_misfits()[::-1]
        isort = history.get_sorted_misfits_idx(chain=0)[::-1]

        gms **= problem.norm_exponent
        gms_softclip = num.where(gms > 1.0, 0.1 * num.log10(gms) + 1.0, gms)

        gcms = problem.combine_misfits(
            history.misfits,
            extra_correlated_weights=optimiser.get_correlated_weights(problem),
            get_contributions=True)

        gcms = gcms[isort, :]
        nmisfits = gcms.shape[1]  # noqa

        ncontributions = sum([1 if t.plot_misfits_cumulative else t.nmisfits
                              for t in problem.targets])
        cum_gcms = num.zeros((history.nmodels, ncontributions))

        # Squash matrix and sum large targets.nmisifts, eg SatelliteTarget
        plot_target_labels = []
        idx = 0
        idx_cum = 0
        for itarget, target in enumerate(problem.targets):
            target_gcms = gcms[:, idx:idx+target.nmisfits]
            if target.plot_misfits_cumulative:
                cum_gcms[:, idx_cum] = target_gcms.sum(axis=1)
                plot_target_labels.append(target.string_id())
                idx_cum += 1
            else:
                cum_gcms[:, idx_cum:idx_cum+target.nmisfits] = target_gcms
                plot_target_labels.extend(target.misfits_string_ids())
                idx_cum += target.nmisfits
            idx += target.nmisfits

        jsort = num.argsort(cum_gcms[-1, :])[::-1]

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
        b = num.hanning(min(100, history.nmodels//5))
        b /= num.sum(b)
        a = [1]
        ii = 0

        for idx in jsort:
            target_label = plot_target_labels[idx]
            ms = cum_gcms[:, idx]

            ms = num.where(num.isfinite(ms), ms, 0.0)
            if num.all(ms == 0.0):
                continue

            rel_ms = ms / gms

            if b.shape[0] > 5:
                rel_ms_smooth = signal.filtfilt(b, a, rel_ms)
            else:
                rel_ms_smooth = rel_ms

            ms_smooth = rel_ms_smooth * gms_softclip

            rel_poly_y = num.concatenate(
                [rel_ms_smooth_sum[::-1], rel_ms_smooth_sum + rel_ms_smooth])
            poly_x = num.concatenate([imodels[::-1], imodels])

            add_args = {}
            if ii < 20:
                add_args['label'] = '%s (%.2g)' % (
                    target_label, num.mean(rel_ms[-1]))

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
    '''
    Sorted misfit (descending) of single bootstrap chains

    For each bootstrap configuration, all models are sorted according to their
    misfit value (red lines) and their global misfit value (black line). (They
    are sorted individually for each line). The best model of every bootstrap
    configuration (right end model of red lines) is marked as a cross in the
    global misfit configuration. The horizontal black lines indicate mean and
    +- standard deviation of the y-axis values of these crosses. If the
    bootstrap configurations converge to the same region in model-space, all
    crosses should be close to the right end of the plot. If this is not the
    case, some bootstrap configurations have converged to very different places
    in model-space. This would be an indicator that there might be
    inconsistencies in the observations (maybe due to faulty or noisy or
    misoriented data). Also the shape of the curve in general can give
    information. A well-behaved optimisation run has approximately linear
    functions in this plot. Only at the end they should have a higher downward
    gradient. This would be the place where the objective functions of the
    bootstrap start to disagree.
    '''

    name = 'bootstrap'
    size_cm = Tuple.T(2, Float.T(), default=(21., 14.9))
    show_ticks = Bool.T(default=False)

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history()
        optimiser = environ.get_optimiser()
        mpl_init(fontsize=self.font_size)
        cm.create_group_mpl(
            self,
            self.draw_figures(history, optimiser),
            title=u'Bootstrap Misfit',
            section='optimiser',
            feather_icon='trending-down',
            description=u'''
Sorted misfit (descending) of single bootstrap chains.

For each bootstrap configuration, all models are sorted according to their
misfit value (red lines) and their global misfit value (black line). (They are
sorted individually for each line). The best model of every bootstrap
configuration (right end model of red lines) is marked as a cross in the global
misfit configuration. The horizontal black lines indicate mean and +- standard
deviation of the y-axis values of these crosses.

If the bootstrap configurations converge to the same region in model-space, all
crosses should be close to the right end of the plot. If this is not the case,
some bootstrap configurations have converged to very different places in
model-space. This would indicate that there might be inconsistencies in the
observations (maybe due to faulty or noisy or misoriented data). Also the shape
of the curve in general can give information. A well-behaved optimisation run
has approximately linear functions in this plot. Only at the end they should
have a higher downward gradient. This would be the place where the objective
functions of the bootstrap start to disagree.
''')

    def draw_figures(self, history, optimiser):

        fig = plt.figure()

        imodels = num.arange(history.nmodels)
        gms = history.bootstrap_misfits[:, 0]

        gms_softclip = num.where(gms > 1.0,
                                 0.1 * num.log10(gms) + 1.0,
                                 gms)
        axes = fig.add_subplot(1, 1, 1)

        ibests = []
        for ibootstrap in range(history.nchains):
            # if ibootstrap ==0:
            #    global, no-bootstrapping misfits, chain
            #    gms = history.bootstrap_misfits[:, ibootstrap]
            #    gms_softclip = num.where(gms > 1.0,
            #                            0.1 * num.log10(gms) + 1.0,
            #                            gms)

            bms = history.bootstrap_misfits[:, ibootstrap]
            isort_bms = num.argsort(bms)[::-1]
            ibests.append(isort_bms[-1])

            bms_softclip = num.where(
                bms > 1.0, 0.1 * num.log10(bms) + 1.0, bms)
            axes.plot(imodels, bms_softclip[isort_bms], color='red', alpha=0.2)

        isort = num.argsort(gms)[::-1]
        iorder = num.empty(isort.size)
        iorder[isort] = imodels

        axes.plot(iorder[ibests], gms_softclip[ibests], 'x', color='black')

        m = num.median(gms_softclip[ibests])
        s = num.std(gms_softclip[ibests])
        axes.axhline(m + s, color='black', alpha=0.5)
        axes.axhline(m, color='black')
        axes.axhline(m - s, color='black', alpha=0.5)

        axes.plot(imodels, gms_softclip[isort], color='black')

        axes.set_xlim(imodels[0], imodels[-1])
        axes.set_xlabel(
            'Tested model, sorted descending by global misfit value')

        return [(PlotItem(name='main'), fig)]


def get_plot_classes():
    return [SequencePlot, ContributionsPlot, BootstrapPlot]
