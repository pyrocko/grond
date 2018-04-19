from __future__ import print_function
import re
import logging
import os
import os.path as op
from collections import defaultdict

from matplotlib import pyplot as plt

from pyrocko import guts, util, gf

from pyrocko.guts import Object, Float, Int, List, Tuple, String, Unicode, Dict

from pyrocko.plot import mpl_init

from grond import meta, core
from grond.problems.base import ModelHistory, load_problem_info
from grond.plot import plotter


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
    book_refs = List.T(Tuple.T(2, StringID.T()))


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


def light(color, factor=0.5):
    return tuple(1-(1-c)*factor for c in color)


def dark(color, factor=0.5):
    return tuple(c*factor for c in color)


def xpop(s, k):
    try:
        s.remove(k)
        return k

    except KeyError:
        return None


plot_dispatch = {
    # 'bootstrap': draw_bootstrap_figure,
    # 'sequence': draw_sequence_figures,
    # 'contributions': draw_contributions_figure,
    # 'jointpar': draw_jointpar_figures,
    # 'histogram': draw_histogram_figures,
    # 'hudson': draw_hudson_figure,
    # 'fits': draw_fits_figures,
    # 'fits_ensemble': draw_fits_ensemble_figures,
    # 'fits_statics': draw_fits_figures_statics,
    # 'solution': draw_solution_figure,
    # 'location': draw_location_figure
}


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

    except plotter.NoPlotterClassAvailable:
        return []
