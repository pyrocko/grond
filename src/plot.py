from __future__ import print_function
import math
import re
import random
import logging
import os
import os.path as op
from collections import defaultdict
import numpy as num
from scipy import signal

from matplotlib import pyplot as plt
from matplotlib import cm, patches

from pyrocko import beachball, guts, util, gf
from pyrocko import hudson

from pyrocko.guts import Object, Float, Int, List, Tuple, String, Unicode, Dict

from pyrocko.plot import mpl_init, mpl_papersize, mpl_margins, \
    mpl_graph_color, mpl_color

from grond import meta, core
from grond.problems.base import ModelHistory, load_problem_info


logger = logging.getLogger('grond.plot')


class StringID(gf.StringID):
    pass


class PlotFormat(Object):

    @property
    def extension(self):
        return self.name

    def get_dpi(self, size_cm):
        return None


class PNG(PlotFormat):
    name = 'png'

    dpi = Float.T(optional=True)
    size_pixels = Int.T(optional=True)
    width_pixels = Int.T(optional=True)
    height_pixels = Int.T(optional=True)

    @property
    def extension(self):
        if self.dpi is not None:
            return 'd%i.png' % self.dpi
        elif self.size_pixels is not None:
            return 's%i.png' % self.size_pixels
        elif self.width_pixels is not None:
            return 'w%i.png' % self.width_pixels
        elif self.height_pixels is not None:
            return 'h%i.png' % self.height_pixels

    def get_dpi(self, size_cm):
        inch = 2.54
        w_cm, h_cm = size_cm
        w_inch, h_inch = w_cm/inch, h_cm/inch
        if self.dpi:
            return self.dpi
        elif self.size_pixels is not None:
            return min(self.size_pixels/w_inch, self.size_pixels/h_inch)
        elif self.width_pixels is not None:
            return self.width_pixels/w_inch
        elif self.height_pixels is not None:
            return self.height_pixels/h_inch


class PDF(PlotFormat):
    name = 'pdf'


class PlotConfig(Object):
    name = 'undefined'
    formats = List.T(PlotFormat.T())
    size_cm = Tuple.T(2, Float.T())


class PlotPage(Object):
    name = StringID.T()
    attributes = Dict.T(StringID.T(), List.T(String.T()))


class PlotBook(Object):
    name = StringID.T()
    variant = StringID.T()
    description = Unicode.T(optional=True)
    formats = List.T(PlotFormat.T())
    size_cm = Tuple.T(2, Float.T())
    pages = List.T(PlotPage.T())

    def filename_image(self, page, format):
        return '%s.%s.%s' % (
            self.name,
            self.variant,
            page.name,
            format.extension)


class BookShelve(Object):
    book_refs = List.T(Tuple.T(StringID.T(), StringID.T()))


class PlotShelveManager(object):

    def __init__(self, path):
        self._path = path
        self.load_shelve()

    def load_shelve(self):
        self._shelve = guts.load(filename=self.path_index())

    def dump_shelve(self):
        guts.dump(filename=self.path_index())

    def path_shelve(self):
        return op.join(self._path, 'plot_shelve.yaml')

    def path_image(self, book, page, format):
        return op.join(self._path, book.filename_image(page, format))

    def path_book(self, book_ref=None, book=None):
        if book_ref is not None:
            book_name, book_variant = book_ref
        else:
            book_name = book.name
            book_variant = book.variant

        return op.join(
            self._path, book_name, book_variant)

    def create_book(self, config, iter_page_figure, **kwargs):
        book = PlotBook(
            formats=guts.clone(config.formats),
            size_cm=config.size_cm,
            name=config.plotname,
            **kwargs)

        path_book = self.path_book(book=book)
        if os.path.exists(path_book):
            self.remove_book(book.name, book.variant)

        for page, fig in iter_page_figure:
            book.pages.append(page)
            for format in book.formats:
                path = self.path_image(book, page, format)
                util.ensuredirs(path)
                fig.savefig(
                    path,
                    format=format.name,
                    dpi=format.get_dpi(book.size_cm))

                logger.info('figure saved: %s' % path)

        book.dump(filename=path_book)
        self.add(book)

    def remove_book(self, book_name, book_variant):
        book = guts.load(filename=self.path_book(book_name, book_variant))
        for page in book.pages:
            for format in book.formats:
                path = self.path_image(book, page, format)
                os.unlink(path)

        path_book = self.path_book(book)
        os.unlink(path_book)


def save_figs(figs, plot_dirname, plotname, formats, dpi):
    for fmt in formats:
        if fmt not in ['pdf', 'png']:
            raise core.GrondError('unavailable output format: %s' % fmt)

    assert re.match(r'^[a-zA-Z0-9_.]+$', plotname)

    # remove files from previous runs
    pat = re.compile(r'^%s-[0-9]+\.(%s)$' % (
        re.escape(plotname), '|'.join(formats)))
    if op.exists(plot_dirname):
        for entry in os.listdir(plot_dirname):
            if pat.match(entry):
                os.unlink(op.join(plot_dirname, entry))

    fns = []
    for ifig, fig in enumerate(figs):
        for format in formats:
            fn = op.join(plot_dirname, '%s-%02i.%s' % (plotname, ifig, format))
            util.ensuredirs(fn)

            fig.savefig(fn, format=format, dpi=dpi)
            logger.info('figure saved: %s' % fn)
            fns.append(fn)

    return fns


class Plotter(object):

    @classmethod
    def draw_summary_figures(cls, sources, target, results, config):
        raise NotImplementedError('to be implemented in subclass')

    @classmethod
    def draw_check_figures(cls, sources, target, results, config):
        raise NotImplementedError('to be implemented in subclass')

    @classmethod
    def draw_result_figures(cls, sources, target, results, config):
        raise NotImplementedError('to be implemented in subclass')


class NoPlotterClassAvailable(Exception):
    pass


def light(color, factor=0.5):
    return tuple(1-(1-c)*factor for c in color)


def dark(color, factor=0.5):
    return tuple(c*factor for c in color)


def fixlim(lo, hi):
    if lo == hi:
        return lo - 1.0, hi + 1.0
    else:
        return lo, hi


def eigh_sorted(mat):
    evals, evecs = num.linalg.eigh(mat)
    iorder = num.argsort(evals)
    return evals[iorder], evecs[:, iorder]


def draw_sequence_figures(
        history, optimizer, plt, misfit_cutoff=None, sort_by='iteration'):

    problem = history.problem
    npar = problem.nparameters
    ndep = problem.ndependants

    imodels = num.arange(history.nmodels)
    bounds = problem.get_combined_bounds()

    xref = problem.get_xref()

    models = history.models

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

    fontsize = 10.0

    nfx = 2
    nfy = 3
    # nfz = (npar + ndep + 1 - 1) / (nfx*nfy) + 1
    cmap = cm.YlOrRd
    cmap = cm.jet
    msize = 1.5
    axes = None
    figs = []
    fig = None
    alpha = 0.5
    for ipar in range(npar):
        impl = ipar % (nfx * nfy) + 1

        if impl == 1:
            fig = plt.figure(figsize=mpl_papersize('a5', 'landscape'))
            labelpos = mpl_margins(fig, nw=nfx, nh=nfy, w=7., h=5., wspace=7.,
                                   hspace=2., units=fontsize)
            figs.append(fig)

        par = problem.parameters[ipar]

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
            labelpos = mpl_margins(fig, nw=nfx, nh=nfy, w=7., h=5., wspace=7.,
                                   hspace=2., units=fontsize)
            figs.append(fig)

        par = problem.dependants[idep]

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
        figs.append(fig)

    axes = fig.add_subplot(nfy, nfx, impl)
    labelpos(axes, 2.5, 2.0)

    config_axes(axes, nfx, nfy, impl, npar + ndep, npar + ndep + 1)

    axes.set_ylim(0., 1.5)
    axes.set_yticks([0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
    axes.set_yticklabels(['0.0', '0.2', '0.4', '0.6', '0.8', '1', '10', '100'])

    axes.scatter(
        imodels[ibest], gms_softclip[ibest], c=iorder[ibest],
        s=msize, edgecolors='none', cmap=cmap, alpha=alpha)

    axes.axhspan(1.0, 1.5, color=(0.8, 0.8, 0.8), alpha=0.2)
    axes.axhline(1.0, color=(0.5, 0.5, 0.5), zorder=2)

    axes.set_xlim(0, history.nmodels)
    axes.set_xlabel('Iteration')

    axes.set_ylabel('Misfit')

    return figs


def draw_jointpar_figures(
        history, optimizer, plt, misfit_cutoff=None, ibootstrap=None,
        color=None, exclude=None, include=None, draw_ellipses=False):

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
                for (xpos, xoff, x) in [(0.0, 10., xmin), (1.0, -10., xmax)]:
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
                for (ypos, yoff, y) in [(0., 10., ymin), (1.0, -10., ymax)]:
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
                    xy=(num.mean(xpar.scaled(fx)), num.mean(ypar.scaled(fy))),
                    width=evals[0] * 2,
                    height=evals[1] * 2,
                    angle=num.rad2deg(num.arctan2(evecs[1][0], evecs[0][0])))

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


def draw_histogram_figures(history, optimizer, plt):

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

            axes.plot(par.scaled(vps), par.inv_scaled(pps), color=stats_color)

        pstats = rstats.parameter_stats_list[iselected]

        axes.axvspan(
            par.scaled(pstats.minimum), par.scaled(pstats.maximum),
            color=stats_color, alpha=0.1)
        axes.axvspan(
            par.scaled(pstats.percentile16), par.scaled(pstats.percentile84),
            color=stats_color, alpha=0.1)
        axes.axvspan(
            par.scaled(pstats.percentile5), par.scaled(pstats.percentile95),
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


def draw_solution_figure(
        history, optimizer, plt, misfit_cutoff=None, beachball_type='full'):

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


def draw_contributions_figure(history, optimizer, plt):

    fontsize = 10.

    fig = plt.figure(figsize=mpl_papersize('a5', 'landscape'))
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

    axes2.set_xlabel('Tested model, sorted descending by global misfit value')

    axes2.set_ylabel('Square of misfit')

    axes2.set_ylim(0., 1.5)
    axes2.axhspan(1.0, 1.5, color=(0.8, 0.8, 0.8))
    axes2.set_yticks([0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5])
    axes2.set_yticklabels(
        ['0.0', '0.2', '0.4', '0.6', '0.8', '1', '10', '100', '1000', '10000',
         '100000'])

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

        # axes.plot(imodels, rel_ms_sum, color='black', alpha=0.1, zorder=-1)

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

    return [fig]


def draw_bootstrap_figure(history, optimizer, plt):

    fig = plt.figure()

    problem = history.problem
    gms = problem.combine_misfits(history.misfits)

    imodels = num.arange(history.nmodels)

    axes = fig.add_subplot(1, 1, 1)

    gms_softclip = num.where(gms > 1.0, 0.1 * num.log10(gms) + 1.0, gms)

    ibests = []
    for ibootstrap in range(optimizer.nbootstrap):
        bms = optimizer.bootstrap_misfits(
            problem, history.misfits, ibootstrap)

        isort_bms = num.argsort(bms)[::-1]

        ibests.append(isort_bms[-1])

        bms_softclip = num.where(bms > 1.0, 0.1 * num.log10(bms) + 1.0, bms)
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
    axes.set_xlabel('Tested model, sorted descending by global misfit value')

    return [fig]


def draw_hudson_figure(history, optimizer, plt):

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


def draw_location_figure(history, optimizer, plt):
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
            norm=colors.Normalize(vmin=num.min(iorder), vmax=num.max(iorder)),
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


def xpop(s, k):
    try:
        s.remove(k)
        return k

    except KeyError:
        return None


plot_dispatch = {
    'bootstrap': draw_bootstrap_figure,
    'sequence': draw_sequence_figures,
    'contributions': draw_contributions_figure,
    'jointpar': draw_jointpar_figures,
    'histogram': draw_histogram_figures,
    'hudson': draw_hudson_figure,
    'fits': draw_fits_figures,
    'fits_ensemble': draw_fits_ensemble_figures,
    'fits_statics': draw_fits_figures_statics,
    'solution': draw_solution_figure,
    'location': draw_location_figure}


def available_plotnames():
    return list(plot_dispatch.keys())


def plot_result(dirname, plotnames_want,
                save=False, save_path=None, formats=('pdf',), dpi=None):

    if isinstance(formats, str):
        formats = formats.split(',')

    plotnames_want = set(plotnames_want)
    plotnames_avail = set(plot_dispatch.keys())

    if save_path is None:
        plot_dirname = op.join(dirname, 'plots')
    else:
        plot_dirname = save_path
        save = True

    unavailable = plotnames_want - plotnames_avail
    if unavailable:
        raise core.GrondError(
            'unavailable plotname: %s' % ', '.join(unavailable))

    fontsize = 10.0

    mpl_init(fontsize=fontsize)
    fns = defaultdict(list)
    config = None

    optimizer_fn = op.join(dirname, 'optimizer.yaml')
    optimizer = guts.load(filename=optimizer_fn)

    if 3 != len({'bootstrap', 'sequence', 'contributions'} - plotnames_want):
        problem = load_problem_info(dirname)
        history = ModelHistory(problem, path=dirname)

        for plotname in ['bootstrap', 'sequence', 'contributions']:
            if plotname in plotnames_want:
                figs = plot_dispatch[plotname](history, optimizer, plt)
                if save:
                    fns[plotname].extend(
                        save_figs(figs, plot_dirname, plotname, formats, dpi))

                    for fig in figs:
                        plt.close(fig)

    if 8 != len({
            'fits',
            'fits_statics',
            'fits_ensemble',
            'jointpar',
            'histogram',
            'hudson',
            'solution',
            'location'} - plotnames_want):

        problem = load_problem_info(dirname)
        history = ModelHistory(problem, path=meta.xjoin(dirname, 'harvest'))

        for plotname in ['fits', 'fits_ensemble', 'fits_statics']:
            if plotname in plotnames_want:
                event_name = problem.base_source.name

                if config is None:
                    config = guts.load(
                        filename=op.join(dirname, 'config.yaml'))
                    config.set_basepath(dirname)
                    config.setup_modelling_environment(problem)

                ds = config.get_dataset(event_name)
                figs = plot_dispatch[plotname](ds, history, optimizer, plt)
                if save:
                    fns[plotname].extend(
                        save_figs(figs, plot_dirname, plotname, formats, dpi))

                    for fig in figs:
                        plt.close(fig)

        for plotname in [
                'jointpar',
                'histogram',
                'hudson',
                'solution',
                'location']:

            if plotname in plotnames_want:
                figs = plot_dispatch[plotname](history, optimizer, plt)
                if save:
                    fns[plotname].extend(
                        save_figs(figs, plot_dirname, plotname, formats, dpi))

                    for fig in figs:
                        plt.close(fig)

    if not save:
        plt.show()

    return fns


def make_movie(dirname, xpar_name, ypar_name, movie_filename):
    optimizer_fn = op.join(dirname, 'optimizer.yaml')
    optimizer = guts.load(filename=optimizer_fn)
    problem = load_problem_info(dirname)
    history = ModelHistory(problem, path=dirname)
    movie_maker = optimizer.get_movie_maker(
        problem, history, xpar_name, ypar_name, movie_filename)

    movie_maker.render()


def draw_target_check_figures(sources, target, results):
    try:
        plotter_class = target.get_plotter_class()
        return plotter_class.draw_check_figures(sources, target, results)
    except NoPlotterClassAvailable:
        return []
