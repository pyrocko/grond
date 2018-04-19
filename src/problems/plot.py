import math
import logging
import random

import numpy as num
from matplotlib import cm, patches

from pyrocko.plot import mpl_margins, mpl_color
from pyrocko.plot import beachball, hudson

from grond.plot import Plotter
from grond import meta, core

logger = logging.getLogger('grond.problem.plot')


def fixlim(lo, hi):
    if lo == hi:
        return lo - 1.0, hi + 1.0
    else:
        return lo, hi


def eigh_sorted(mat):
    evals, evecs = num.linalg.eigh(mat)
    iorder = num.argsort(evals)
    return evals[iorder], evecs[:, iorder]


class ProblemPlotter(Plotter):

    @classmethod
    @Plotter.register('results', 'jointpar')
    def plot_jointpar(cls):
        pass

    @classmethod
    @Plotter.register('results', 'histogram')
    def plot_histogram(cls):
        pass

    @classmethod
    @Plotter.register('results', 'location')
    def plot_location(cls):
        pass

    @classmethod
    @Plotter.register('results', 'solution')
    def plot_solution(cls):
        pass

    @classmethod
    @Plotter.register('results', 'hudson')
    def plot_hudson(cls):
        pass

    @classmethod
    def draw_jointpar_figures(
            cls, history, optimizer, plt,
            misfit_cutoff=None,
            ibootstrap=None,
            color=None,
            exclude=None,
            include=None,
            draw_ellipses=False):

        color = 'misfit'
        # exclude = ['duration']
        # include = ['magnitude', 'rel_moment_iso', 'rel_moment_clvd', 'depth']
        neach = 6
        figsize = (8, 8)
        # cmap = cm.YlOrRd
        # cmap = cm.jet
        msize = 1.5

        problem = history.problem
        if not problem:
            return []

        models = history.models

        bounds = problem.get_combined_bounds()
        for ipar in range(problem.ncombined):
            par = problem.combined[ipar]
            lo, hi = bounds[ipar]
            if lo == hi:
                if exclude is None:
                    exclude = []

                exclude.append(par.name)

        xref = problem.get_xref()

        if ibootstrap is not None:
            gms = optimizer.bootstrap_misfits(
                problem, history.misfits, ibootstrap)
        else:
            gms = problem.combine_misfits(history.misfits)

        isort = num.argsort(gms)[::-1]

        gms = gms[isort]
        models = models[isort, :]

        if misfit_cutoff is not None:
            ibest = gms < misfit_cutoff
            gms = gms[ibest]
            models = models[ibest]

        nmodels = models.shape[0]

        if color == 'dist':
            mx = num.mean(models, axis=0)
            cov = num.cov(models.T)
            mdists = core.mahalanobis_distance(models, mx, cov)
            icolor = meta.ordersort(mdists)

        elif color == 'misfit':
            iorder = num.arange(nmodels)
            icolor = iorder

        elif color in problem.parameter_names:
            ind = problem.name_to_index(color)
            icolor = problem.extract(models, ind)

        smap = {}
        iselected = 0
        for ipar in range(problem.ncombined):
            par = problem.combined[ipar]
            if exclude and par.name in exclude or \
                    include and par.name not in include:
                continue

            smap[iselected] = ipar
            iselected += 1

        nselected = iselected

        if nselected < 2:
            logger.warn('cannot draw joinpar figures with less than two '
                        'parameters selected')
            return []

        nfig = (nselected - 2) // neach + 1

        figs = []
        for ifig in range(nfig):
            figs_row = []
            for jfig in range(nfig):
                if ifig >= jfig:
                    figs_row.append(plt.figure(figsize=figsize))
                else:
                    figs_row.append(None)

            figs.append(figs_row)

        for iselected in range(nselected):
            ipar = smap[iselected]
            ypar = problem.combined[ipar]
            for jselected in range(iselected):
                jpar = smap[jselected]
                xpar = problem.combined[jpar]

                ixg = (iselected - 1)
                iyg = jselected

                ix = ixg % neach
                iy = iyg % neach

                ifig = ixg // neach
                jfig = iyg // neach

                aind = (neach, neach, (ix * neach) + iy + 1)

                fig = figs[ifig][jfig]

                axes = fig.add_subplot(*aind)

                axes.axvline(0., color=mpl_color('aluminium3'), lw=0.5)
                axes.axhline(0., color=mpl_color('aluminium3'), lw=0.5)
                for spine in axes.spines.values():
                    spine.set_edgecolor(mpl_color('aluminium5'))
                    spine.set_linewidth(0.5)

                xmin, xmax = fixlim(*xpar.scaled(bounds[jpar]))
                ymin, ymax = fixlim(*ypar.scaled(bounds[ipar]))

                if ix == 0 or jselected + 1 == iselected:
                    for (xpos, xoff, x) in [
                            (0.0, 10., xmin),
                            (1.0, -10., xmax)]:

                        axes.annotate(
                            '%.3g%s' % (x, xpar.get_unit_suffix()),
                            xy=(xpos, 1.05),
                            xycoords='axes fraction',
                            xytext=(xoff, 5.),
                            textcoords='offset points',
                            verticalalignment='bottom',
                            horizontalalignment='left',
                            rotation=45.)

                if iy == neach - 1 or jselected + 1 == iselected:
                    for (ypos, yoff, y) in [
                            (0., 10., ymin),
                            (1.0, -10., ymax)]:

                        axes.annotate(
                            '%.3g%s' % (y, ypar.get_unit_suffix()),
                            xy=(1.0, ypos),
                            xycoords='axes fraction',
                            xytext=(5., yoff),
                            textcoords='offset points',
                            verticalalignment='bottom',
                            horizontalalignment='left',
                            rotation=45.)

                axes.set_xlim(xmin, xmax)
                axes.set_ylim(ymin, ymax)

                axes.get_xaxis().set_ticks([])
                axes.get_yaxis().set_ticks([])

                if iselected == nselected - 1 or ix == neach - 1:
                    axes.annotate(
                        xpar.get_label(with_unit=False),
                        xy=(0.5, -0.05),
                        xycoords='axes fraction',
                        verticalalignment='top',
                        horizontalalignment='right',
                        rotation=45.)

                if iy == 0:
                    axes.annotate(
                        ypar.get_label(with_unit=False),
                        xy=(-0.05, 0.5),
                        xycoords='axes fraction',
                        verticalalignment='top',
                        horizontalalignment='right',
                        rotation=45.)

                fx = problem.extract(models, jpar)
                fy = problem.extract(models, ipar)

                axes.scatter(
                    xpar.scaled(fx),
                    ypar.scaled(fy),
                    c=icolor,
                    s=msize, alpha=0.5, cmap='coolwarm', edgecolors='none')

                if draw_ellipses:
                    cov = num.cov((xpar.scaled(fx), ypar.scaled(fy)))
                    evals, evecs = eigh_sorted(cov)
                    evals = num.sqrt(evals)
                    ell = patches.Ellipse(
                        xy=(
                            num.mean(xpar.scaled(fx)),
                            num.mean(ypar.scaled(fy))),
                        width=evals[0] * 2,
                        height=evals[1] * 2,
                        angle=num.rad2deg(
                            num.arctan2(evecs[1][0], evecs[0][0])))

                    ell.set_facecolor('none')
                    axes.add_artist(ell)

                fx = problem.extract(xref, jpar)
                fy = problem.extract(xref, ipar)

                ref_color = mpl_color('aluminium6')
                ref_color_light = 'none'
                axes.plot(
                    xpar.scaled(fx), ypar.scaled(fy), 's',
                    mew=1.5, ms=5, mfc=ref_color_light, mec=ref_color)

        figs_flat = []
        for figs_row in figs:
            figs_flat.extend(fig for fig in figs_row if fig is not None)

        return figs_flat

    @classmethod
    def draw_histogram_figures(cls, history, optimizer, plt):

        import scipy.stats
        from grond.core import make_stats

        exclude = None
        include = None
        figsize = (5, 3)
        fontsize = 10
        method = 'gaussian_kde'
        ref_color = mpl_color('aluminium6')
        stats_color = mpl_color('scarletred2')

        problem = history.problem
        misfits = history.misfits
        if not problem:
            return []

        models = history.models

        bounds = problem.get_combined_bounds()
        for ipar in range(problem.ncombined):
            par = problem.combined[ipar]
            vmin, vmax = bounds[ipar]
            if vmin == vmax:
                if exclude is None:
                    exclude = []

                exclude.append(par.name)

        xref = problem.get_xref()

        smap = {}
        iselected = 0
        for ipar in range(problem.ncombined):
            par = problem.combined[ipar]
            if exclude and par.name in exclude or \
                    include and par.name not in include:
                continue

            smap[iselected] = ipar
            iselected += 1

        nselected = iselected

        pnames = [
            problem.combined[smap[iselected]].name
            for iselected in range(nselected)]

        rstats = make_stats(problem, models, misfits, pnames=pnames)

        figs = []
        for iselected in range(nselected):
            ipar = smap[iselected]
            par = problem.combined[ipar]
            vs = problem.extract(models, ipar)
            vmin, vmax = bounds[ipar]

            fig = plt.figure(figsize=figsize)
            labelpos = mpl_margins(
                fig, nw=1, nh=1, w=7., bottom=5., top=1, units=fontsize)

            axes = fig.add_subplot(1, 1, 1)
            labelpos(axes, 2.5, 2.0)
            axes.set_xlabel(par.get_label())
            axes.set_ylabel('PDF')
            axes.set_xlim(*fixlim(*par.scaled((vmin, vmax))))

            for method in 'gaussian_kde', 'histogram':
                if method == 'gaussian_kde':
                    kde = scipy.stats.gaussian_kde(vs)
                    vps = num.linspace(vmin, vmax, 600)
                    pps = kde(vps)

                elif method == 'histogram':
                    pps, edges = num.histogram(vs, density=True)
                    vps = 0.5 * (edges[:-1] + edges[1:])

                axes.plot(
                    par.scaled(vps), par.inv_scaled(pps), color=stats_color)

            pstats = rstats.parameter_stats_list[iselected]

            axes.axvspan(
                par.scaled(pstats.minimum),
                par.scaled(pstats.maximum),
                color=stats_color, alpha=0.1)
            axes.axvspan(
                par.scaled(pstats.percentile16),
                par.scaled(pstats.percentile84),
                color=stats_color, alpha=0.1)
            axes.axvspan(
                par.scaled(pstats.percentile5),
                par.scaled(pstats.percentile95),
                color=stats_color, alpha=0.1)

            axes.axvline(
                par.scaled(pstats.median),
                color=stats_color, alpha=0.5)
            axes.axvline(
                par.scaled(pstats.mean),
                color=stats_color, ls=':', alpha=0.5)

            axes.axvline(
                par.scaled(problem.extract(xref, ipar)),
                color=ref_color)

            figs.append(fig)

        return figs

    @classmethod
    def draw_solution_figure(
            cls, history, optimizer, plt,
            misfit_cutoff=None,
            beachball_type='full'):

        fontsize = 10.

        fig = plt.figure(figsize=(6, 2))
        axes = fig.add_subplot(1, 1, 1, aspect=1.0)
        fig.subplots_adjust(left=0., right=1., bottom=0., top=1.)

        problem = history.problem
        if not problem:
            logger.warn('problem not set')
            return []

        models = history.models

        if models.size == 0:
            logger.warn('empty models vector')
            return []

        gms = problem.combine_misfits(history.misfits)
        isort = num.argsort(gms)
        iorder = num.empty_like(isort)
        iorder[isort] = num.arange(iorder.size)[::-1]

        mean_source = core.get_mean_source(
            problem, history.models)

        best_source = core.get_best_source(
            problem, history.models, history.misfits)

        ref_source = problem.base_source

        for xpos, label in [
                (0., 'Full'),
                (2., 'Isotropic'),
                (4., 'Deviatoric'),
                (6., 'CLVD'),
                (8., 'DC')]:

            axes.annotate(
                label,
                xy=(1 + xpos, 3),
                xycoords='data',
                xytext=(0., 0.),
                textcoords='offset points',
                ha='center',
                va='center',
                color='black',
                fontsize=fontsize)

        decos = []
        for source in [best_source, mean_source, ref_source]:
            mt = source.pyrocko_moment_tensor()
            deco = mt.standard_decomposition()
            decos.append(deco)

        moment_full_max = max(deco[-1][0] for deco in decos)

        for ypos, label, deco, color_t in [
                (2., 'Ensemble best', decos[0], mpl_color('aluminium5')),
                (1., 'Ensemble mean', decos[1], mpl_color('scarletred1')),
                (0., 'Reference', decos[2], mpl_color('aluminium3'))]:

            [(moment_iso, ratio_iso, m_iso),
             (moment_dc, ratio_dc, m_dc),
             (moment_clvd, ratio_clvd, m_clvd),
             (moment_devi, ratio_devi, m_devi),
             (moment_full, ratio_full, m_full)] = deco

            size0 = moment_full / moment_full_max

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

            for xpos, mt_part, ratio, ops in [
                    (0., m_full, ratio_full, '-'),
                    (2., m_iso, ratio_iso, '='),
                    (4., m_devi, ratio_devi, '='),
                    (6., m_clvd, ratio_clvd, '+'),
                    (8., m_dc, ratio_dc, None)]:

                if ratio > 1e-4:
                    try:
                        beachball.plot_beachball_mpl(
                            mt_part, axes,
                            beachball_type='full',
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
        axes.set_ylim(-0.5, 3.5)

        return [fig]

    @classmethod
    def draw_location_figure(cls, history, optimizer, plt):
        from matplotlib import colors

        color = 'black'
        fontsize = 10.
        markersize = fontsize * 1.5
        # markersize_small = markersize * 0.2
        beachballsize = markersize
        beachballsize_small = beachballsize * 0.5
        beachball_type = 'dc'
        width = 7.
        figsize = (width, width / (4. / 3.))

        problem = history.problem

        fig = plt.figure(figsize=figsize)
        axes_en = fig.add_subplot(2, 2, 1)
        axes_dn = fig.add_subplot(2, 2, 2)
        axes_ed = fig.add_subplot(2, 2, 3)

        bounds = problem.get_combined_bounds()

        gms = problem.combine_misfits(history.misfits)

        isort = num.argsort(gms)[::-1]

        gms = gms[isort]
        models = history.models[isort, :]

        iorder = num.arange(history.nmodels)

        for axes, xparname, yparname in [
                (axes_en, 'east_shift', 'north_shift'),
                (axes_dn, 'depth', 'north_shift'),
                (axes_ed, 'east_shift', 'depth')]:

            ixpar = problem.name_to_index(xparname)
            iypar = problem.name_to_index(yparname)

            xpar = problem.combined[ixpar]
            ypar = problem.combined[iypar]

            axes.set_xlabel(xpar.get_label())
            axes.set_ylabel(ypar.get_label())

            xmin, xmax = fixlim(*xpar.scaled(bounds[ixpar]))
            ymin, ymax = fixlim(*ypar.scaled(bounds[iypar]))

            axes.set_aspect(1.0)
            axes.set_xlim(xmin, xmax)
            axes.set_ylim(ymin, ymax)

            # fxs = xpar.scaled(problem.extract(models, ixpar))
            # fys = ypar.scaled(problem.extract(models, iypar))

            # axes.set_xlim(*fixlim(num.min(fxs), num.max(fxs)))
            # axes.set_ylim(*fixlim(num.min(fys), num.max(fys)))

            cmap = cm.ScalarMappable(
                norm=colors.Normalize(
                    vmin=num.min(iorder),
                    vmax=num.max(iorder)),

                cmap=plt.get_cmap('coolwarm'))

            for ix, x in enumerate(models):
                source = problem.get_source(x)
                mt = source.pyrocko_moment_tensor()
                fx = problem.extract(x, ixpar)
                fy = problem.extract(x, iypar)
                sx, sy = xpar.scaled(fx), ypar.scaled(fy)

                color = cmap.to_rgba(iorder[ix])

                alpha = (iorder[ix] - num.min(iorder)) / \
                    float(num.max(iorder) - num.min(iorder))

                try:
                    beachball.plot_beachball_mpl(
                        mt, axes,
                        beachball_type=beachball_type,
                        position=(sx, sy),
                        size=beachballsize_small,
                        color_t=color,
                        alpha=alpha,
                        zorder=1,
                        linewidth=0.25)

                except beachball.BeachballError as e:
                    logger.warn(str(e))

        return [fig]

    @classmethod
    def draw_hudson_figure(cls, history, optimizer, plt):

        color = 'black'
        fontsize = 10.
        markersize = fontsize * 1.5
        markersize_small = markersize * 0.2
        beachballsize = markersize
        beachballsize_small = beachballsize * 0.5
        width = 7.
        beachball_type = 'dc'
        figsize = (width, width / (4. / 3.))

        problem = history.problem
        mean_source = core.get_mean_source(problem, history.models)
        best_source = core.get_best_source(
            problem, history.models, history.misfits)

        fig = plt.figure(figsize=figsize)
        axes = fig.add_subplot(1, 1, 1)

        data = []
        for ix, x in enumerate(history.models):
            source = problem.get_source(x)
            mt = source.pyrocko_moment_tensor()
            u, v = hudson.project(mt)

            if random.random() < 0.1:
                try:
                    beachball.plot_beachball_mpl(
                        mt, axes,
                        beachball_type=beachball_type,
                        position=(u, v),
                        size=beachballsize_small,
                        color_t=color,
                        alpha=0.5,
                        zorder=1,
                        linewidth=0.25)
                except beachball.BeachballError as e:
                    logger.warn(str(e))

            else:
                data.append((u, v))

        if data:
            u, v = num.array(data).T
            axes.plot(
                u, v, 'o',
                color=color,
                ms=markersize_small,
                mec='none',
                mew=0,
                alpha=0.25,
                zorder=0)

        hudson.draw_axes(axes)

        mt = mean_source.pyrocko_moment_tensor()
        u, v = hudson.project(mt)

        try:
            beachball.plot_beachball_mpl(
                mt, axes,
                beachball_type=beachball_type,
                position=(u, v),
                size=beachballsize,
                color_t=color,
                zorder=2,
                linewidth=0.5)
        except beachball.BeachballError as e:
            logger.warn(str(e))

        mt = best_source.pyrocko_moment_tensor()
        u, v = hudson.project(mt)

        axes.plot(
            u, v, 's',
            markersize=markersize,
            mew=1,
            mec='black',
            mfc='none',
            zorder=-2)

        mt = problem.base_source.pyrocko_moment_tensor()
        u, v = hudson.project(mt)

        try:
            beachball.plot_beachball_mpl(
                mt, axes,
                beachball_type=beachball_type,
                position=(u, v),
                size=beachballsize,
                color_t='red',
                zorder=2,
                linewidth=0.5)
        except beachball.BeachballError as e:
            logger.warn(str(e))

        return [fig]
