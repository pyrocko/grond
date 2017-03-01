import tornado.ioloop
import grond
import os
import os.path as op
import logging
import numpy as num
import socket

from collections import OrderedDict

from pyrocko.guts import Object, Bool, String, Int, List

from tornado.web import RequestHandler, StaticFileHandler
from tornado import gen

from bokeh.embed import autoload_server
from bokeh.application import Application
from bokeh.server.server import Server as BokehServer
from bokeh.application.handlers import Handler as BokehHandler

from bokeh.models import ColumnDataSource
from bokeh import layouts
from bokeh.plotting import figure

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('grond.baraddur')


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


class GrondBokehModel(object):
    def __init__(self, config):
        self.config = config
        self.set_rundir(self.config.rundir)

    def set_rundir(self, rundir):
        logger.debug('Loading problem from %s' % rundir)
        self.rundir = rundir
        self.problem = grond.core.load_problem_info(self.rundir)
        self.parameters = self.problem.parameters
        self.nparameters = self.problem.nparameters
        self.ntargets = self.problem.ntargets

    def get_models(self, skip_nmodels=0):
        fn = op.join(self.rundir, 'models')
        with open(fn, 'r') as f:
            nmodels = os.fstat(f.fileno()).st_size / (self.nparameters * 8)
            nmodels -= skip_nmodels
            f.seek(skip_nmodels * self.nparameters * 8)
            data = num.fromfile(
                f, dtype='<f8', count=nmodels * self.nparameters)\
                .astype(num.float)

        nmodels = data.size/self.nparameters - skip_nmodels
        models = data.reshape((nmodels, self.nparameters))

        mods_dict = {}
        for ip, par in enumerate(self.parameters):
            mods_dict[par.name] = models[:, ip]
        mods_dict['niter'] = num.arange(nmodels, dtype=num.int) + (nmodels+1)
        return nmodels, mods_dict

    def get_misfits(self, skip_nmodels=0):
        fn = op.join(self.rundir, 'misfits')

        with open(fn, 'r') as f:
            nmodels = os.fstat(f.fileno()).st_size / (self.nparameters * 8)
            nmodels -= skip_nmodels
            f.seek(skip_nmodels * self.ntargets * 2 * 8)
            data = num.fromfile(
                f, dtype='<f8', count=nmodels*self.ntargets*2)\
                .astype(num.float)

        data = data.reshape((nmodels, self.ntargets*2))

        combi = num.empty_like(data)
        combi[:, 0::2] = data[:, :self.ntargets]
        combi[:, 1::2] = data[:, self.ntargets:]

        assert(data.size/self.nparameters - skip_nmodels == nmodels)
        misfits = combi.reshape((nmodels, self.ntargets, 2))

        mf_dict = {}
        for it in xrange(self.ntargets):
            mf_dict['target_%03d' % it] = misfits[:, it, 0]
        mf_dict['target_mean'] = num.mean(misfits[:, :, 0])

        mf_dict['niter'] = num.arange(nmodels, dtype=num.int) + (nmodels+1)
        return nmodels, misfits


class Status(BaraddurRequestHandler):

    class MisfitsPlot(BaraddurBokehHandler):

        def modify_document(self, doc):
            self.nmodels = 0
            self.source = ColumnDataSource(
                data={'n': num.ndarray(0),
                      'gm': num.ndarray(0)})
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
                self.source.add(num.ndarray(0), p)
            self.update_parameters()

            plots = []
            for par in problem.parameters:
                fig = figure(webgl=True,
                             x_axis_label='Iteration #',
                             y_axis_label='%s [%s]' % (par.label, par.unit))
                fig.scatter('n', par.name,
                            source=self.source, alpha=.4)
                plots.append(fig)
            plots += [None] * (self.ncols - (len(plots) % self.ncols))
            print plots

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


class Targets(BaraddurRequestHandler):

    class TargetContributionPlot(BaraddurBokehHandler):

        def modify_document(self, doc):
            self.nmodels = 0
            self.source = ColumnDataSource()
            self.update_contributions()

            plot = figure(webgl=True,
                          x_axis_label='Iteration #',
                          y_axis_label='Misfit')

            doc.add_root(plot)
            doc.add_periodic_callback(self.update_contributions, 1e3)

        @gen.coroutine
        def update_contributions(self):
            mx, misfits = grond.core.load_problem_data(
                self.config.rundir, self.config.problem,
                skip_models=self.nmodels)

            print misfits
            new_nmodels = mx.shape[0]

            # self.source.stream(dict(m=misfits[:, :, 0],
            #                         n=num.arange(new_nmodels,
            #                                      dtype=num.int) +
            #                         self.nmodels + 1))
            self.nmodels += new_nmodels

    bokeh_handlers = {'contribution_plot': TargetContributionPlot}

    @gen.coroutine
    def get(self):
        self.render('targets.html',
                    contribution_plot=autoload_server(
                        None,
                        url='/contribution_plot'),
                    pages=pages,
                    problem=self.config.problem)


pages = OrderedDict([
    ('Summary', Summary),
    ('Status', Status),
    ('Parameters', Parameters),
    ('Targets', Targets),
])


class BaraddurConfig(Object):
    rundir = String.T(
        help='Grond rundir.')
    debug = Bool.T(
        default=False,
        optional=True)
    hosts = List.T(
        String.T(),
        default=['*'],
        optional=True,
        help='List of allowed hosts, default is all \'*\'.')
    port = Int.T(
        default=8080,
        optional=True,
        help='Port to listen on.')

    @property
    def problem(self):
        return grond.core.load_problem_info(self.rundir)


class Baraddur(BokehServer):
    def __init__(self, rundir=None, *args, **kwargs):
        self.config = BaraddurConfig(rundir=rundir)
        print self.config
        self.ioloop = tornado.ioloop.IOLoop.current()
        port_offset = 0

        while True:
            try:
                BokehServer.__init__(
                    self,
                    self.get_bokeh_apps(),
                    io_loop=self.ioloop,
                    extra_patterns=self.get_tornado_handlers(),
                    port=self.config.port + port_offset,
                    host=self.config.hosts)
                break
            except socket.error as se:
                if se.errno == 98 and port_offset < 50:  # Port in use
                    port_offset += 1
                    logger.info('Port %d in use, bailing to %d'
                                % (self.config.port + port_offset - 1,
                                   self.config.port + port_offset))
                else:
                    raise se
        tornado_app = self._tornado
        tornado_app.settings['template_path'] = op.join(
            op.dirname(op.abspath(__file__)), 'templates')

        if self.config.debug:
            tornado_app.settings.setdefault('autoreload', True)
            tornado_app.settings.setdefault('compiled_template_cache', False)
            tornado_app.settings.setdefault('static_hash_cache', False)
            tornado_app.settings.setdefault('serve_traceback', True)
            # Automatically reload modified modules
            from tornado import autoreload
            autoreload.start()

            logging.getLogger('').setLevel(logging.DEBUG)

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

    def start(self, signal=None):
        logger.info('Starting Baraddur server on http://localhost:%d'
                    % (self.port))

        if signal is not None:
            def shutdown():
                if not signal.empty():
                    self.stop()
            tornado.ioloop.PeriodicCallback(shutdown, 2000).start()

        BokehServer.start(self)
        self.ioloop.start()

    def stop(self, *args, **kwargs):
        print args, kwargs
        logger.info('Stopping Baraddur server...')
        BokehServer.stop(self)
        self.ioloop.stop()


if __name__ == '_123_main__':
    baraddur = Baraddur(
        rundir='/home/marius/Development/testing/grond/rundir')
    baraddur.start()
    print 'here!'
