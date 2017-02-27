import tornado.ioloop
import grond
import os.path as op
import logging
import numpy as num

from collections import OrderedDict

from pyrocko.guts import Object, Bool, String  # noqa

from tornado.web import RequestHandler, StaticFileHandler
from tornado import gen

from bokeh.embed import autoload_server
from bokeh.application import Application
from bokeh.server.server import Server as BokehServer
from bokeh.application.handlers import Handler as BokehHandler

from bokeh.models import ColumnDataSource
from bokeh import layouts
from bokeh.plotting import figure


def makeColorGradient(misfits, fr=1., fg=.5, fb=1.,
                      pr=0, pg=2.5, pb=4):
    misfits /= misfits.max()
    r = num.sin(fr * misfits + pr) * 127 + 128
    g = num.sin(fg * misfits + pg) * 127 + 128
    b = num.sin(fb * misfits + pb) * 127 + 128
    return ['#%02x%02x%02x' % (r[i], g[i], b[i]) for i in xrange(misfits.size)]


class BaraddurRequestHandler(RequestHandler):
    def initialize(self, config):
        self.config = config


class BaraddurBokehHandler(BokehHandler):
    def __init__(self, config, *args, **kwargs):
        BokehHandler.__init__(self, *args, **kwargs)
        self.config = config


class Status(BaraddurRequestHandler):

    class MisfitsPlot(BaraddurBokehHandler):

        def modify_document(self, doc):
            self.nmodels = 0
            self.source = ColumnDataSource(
                data={'n': [],
                      'gm': []})
            self.update_misfits()

            plot = figure(webgl=True,
                          x_axis_label='Iteration #',
                          y_axis_label='Misfit')
            plot.scatter('n', 'gm',
                         source=self.source, alpha=.4)

            doc.add_root(plot)
            doc.add_periodic_callback(self.update_misfits, 1e3)

        @gen.coroutine
        def update_misfits(self):
            mx, misfits = grond.core.load_problem_data(
                self.config.rundir, self.config.problem,
                skip_models=self.nmodels)
            new_nmodels = mx.shape[0]

            fits = num.mean(misfits, axis=1)
            self.source.stream(dict(gm=fits[:, 0],
                                    n=num.arange(new_nmodels,
                                                 dtype=num.int) +
                                    self.nmodels + 1))
            self.nmodels += new_nmodels

    bokeh_handlers = {'misfit_plot': MisfitsPlot}

    @gen.coroutine
    def get(self):
        self.render('status.html',
                    pages=pages,
                    misfit_plot=autoload_server(None, url='/misfit_plot'),
                    problem=self.config.problem)


class Parameters(BaraddurRequestHandler):

    class ParameterPlots(BaraddurBokehHandler):

        ncols = 4

        def modify_document(self, doc):
            self.nmodels = 0
            problem = self.config.problem

            self.source = ColumnDataSource()
            for p in ['n'] + [p.name for p in problem.parameters]:
                self.source.add([], p)
            self.update_parameters()

            plots = []
            for par in problem.parameters:
                fig = figure(webgl=True,
                             x_axis_label='Iteration #',
                             y_axis_label='%s [%s]' % (par.label, par.unit))
                fig.scatter('n', par.name,
                            source=self.source, alpha=.4)
                plots.append(fig)
            plots += ([None] * (self.ncols - (len(plots) % self.ncols)))

            grid = layouts.gridplot(
                plots,
                responsive=True,
                ncols=self.ncols)

            doc.add_root(grid)
            doc.add_periodic_callback(self.update_parameters, 2.5*1e3)

        @gen.coroutine
        def update_parameters(self):
            problem = self.config.problem

            try:
                mx, misfits = grond.core.load_problem_data(
                    self.config.rundir, problem, skip_models=self.nmodels)
            except IOError:
                return

            new_nmodels = mx.shape[0]

            new_data = {}
            for ip, par in enumerate(problem.parameters):
                new_data[par.name] = mx[:, ip]
            new_data['n'] = num.arange(new_nmodels, dtype=num.int) +\
                self.nmodels + 1

            self.source.stream(new_data)
            self.nmodels += new_nmodels

    bokeh_handlers = {'parameter_plot': ParameterPlots}

    @gen.coroutine
    def get(self):

        self.render('parameter_plots.html',
                    pages=pages,
                    parameter_plot=autoload_server(
                        None,
                        url='/parameter_plot'),
                    problem=self.config.problem)


class Summary(BaraddurRequestHandler):

    @gen.coroutine
    def get(self):
        self.render('summary.html',
                    pages=pages,
                    problem=self.config.problem)


class MainHandler(BaraddurRequestHandler):

    @gen.coroutine
    def get(self):
        self.render('index.html',
                    pages=pages)


pages = OrderedDict([
    ('Status', Status),
    ('Parameters', Parameters),
    ('Summary', Summary),
    ('Sequences', MainHandler),
])


class BaraddurConfig(Object):
    rundir = String.T()
    template_path = String.T(
        default='tpl',
        optional=True)
    debug = Bool.T(
        default=True,
        optional=True)

    @property
    def problem(self):
        return grond.core.load_problem_info(self.rundir)


class Baraddur(BokehServer):
    def __init__(self, rundir=None, *args, **kwargs):
        self.config = BaraddurConfig(rundir=rundir)

        BokehServer.__init__(
            self,
            self.get_bokeh_apps(),
            io_loop=tornado.ioloop.IOLoop.current(),
            extra_patterns=self.get_tornado_handlers(),
            hosts='0.0.0.0')

        tornado_app = self._tornado
        tornado_app.listen(8080)
        tornado_app.settings['template_path'] = self.config.template_path

        if self.config.debug:
            tornado_app.settings.setdefault('autoreload', True)
            tornado_app.settings.setdefault('compiled_template_cache', False)
            tornado_app.settings.setdefault('static_hash_cache', False)
            tornado_app.settings.setdefault('serve_traceback', True)
            # Automatically reload modified modules
            from tornado import autoreload
            autoreload.start()

            console = logging.StreamHandler()
            console.setLevel(logging.DEBUG)
            formatter = logging.Formatter(
                '%(name)-12s: %(levelname)-8s %(message)s')
            console.setFormatter(formatter)
            rootLogger = logging.getLogger('')
            rootLogger.addHandler(console)
            rootLogger.setLevel(logging.DEBUG)

    def get_bokeh_apps(self):
        bokeh_apps = {}
        for tornado_handler in pages.itervalues():
            handler_docs = getattr(tornado_handler, 'bokeh_handlers',
                                   None)
            if handler_docs is None:
                continue

            for url, bokeh_handler in handler_docs.iteritems():
                bokeh_apps['/%s' % url] = Application(bokeh_handler(
                    self.config))
        return bokeh_apps

    def get_tornado_handlers(self):
        return [(r'/', pages.values()[0],
                 {'config': self.config})] +\
               [(r'/%s' % title, handler,
                 {'config': self.config})
                for title, handler in pages.iteritems()] +\
               [(r'/css/(.*)', StaticFileHandler,
                {'path': op.join(op.dirname(__file__), 'css')})]

    def start(self):
        BokehServer.start(self)
        tornado.ioloop.IOLoop.current().start()


if __name__ == '__main__':
    baraddur = Baraddur(
        rundir='/home/marius/Development/testing/grond/rundir')
    baraddur.start()
