import tornado.ioloop
import grond
import os
import glob
import os.path as op
import logging
import numpy as num
import socket

from scipy import ndimage

from collections import OrderedDict
from datetime import datetime

from pyrocko.guts import Object, Bool, String, Int, List, Float

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


class BaraddurModel(object):
    def __init__(self, rundir):
        logger.debug('Loading problem from %s' % rundir)
        self.rundir = op.abspath(rundir)
        self.problem = grond.core.load_problem_info(self.rundir)
        self.parameters = self.problem.parameters
        self.nparameters = self.problem.nparameters
        self.targets = self.problem.targets
        self.ntargets = self.problem.ntargets

    def _jp(self, file):
        return op.join(self.rundir, file)

    @property
    def start_time(self):
        return datetime.fromtimestamp(
                os.stat(self._jp('problem.yaml')).st_mtime)

    @property
    def duration(self):
        return datetime.fromtimestamp(
            os.stat(self._jp('models')).st_mtime) - self.start_time

    @property
    def name(self):
        return op.split(self.rundir)[-1]

    @property
    def best_misfit(self):
        _, data = self.get_misfits(keys='misfit_mean')
        return data['misfit_mean'].min()

    @property
    def niterations(self):
        return op.getsize(self._jp('models')) / (self.nparameters * 8)

    @property
    def iter_per_second(self):
        return self.niterations / self.duration.seconds

    @staticmethod
    def validate_rundir(rundir):

        def jp(fn):
            return op.join(rundir, fn)

        for fn in ['models', 'misfits', 'problem.yaml']:
            if op.lexists(jp(fn)):
                if op.getsize(jp(fn)) == 0:
                    return False
            else:
                return False
        return True

    def get_models(self, skip_nmodels=0, keys=None):
        fn = self._jp('models')

        with open(fn, 'r') as f:
            nmodels = os.fstat(f.fileno()).st_size / (self.nparameters * 8)
            nmodels -= skip_nmodels
            f.seek(skip_nmodels * self.nparameters * 8)
            data = num.fromfile(
                f, dtype='<f8', count=nmodels * self.nparameters)\
                .astype(num.float32)

        nmodels = data.size/self.nparameters
        models = data.reshape((nmodels, self.nparameters))

        models_dict = {}
        for ip, par in enumerate(self.parameters):
            models_dict[par.name] = models[:, ip]
        models_dict['niter'] = num.arange(nmodels, dtype=num.int)\
            + skip_nmodels + 1

        if keys is not None:
            for k in models_dict.keys():
                if k not in keys:
                    models_dict.pop(k)

        return nmodels, models_dict

    def get_misfits(self, skip_nmodels=0, keys=None):
        fn = self._jp('misfits')

        with open(fn, 'r') as f:
            nmodels = os.fstat(f.fileno()).st_size / (self.ntargets * 2 * 8)
            nmodels -= skip_nmodels
            f.seek(skip_nmodels * self.ntargets * 2 * 8)
            data = num.fromfile(
                f, dtype='<f8', count=nmodels*self.ntargets*2)\
                .astype(num.float32)

        nmodels = data.size/(self.ntargets * 2)

        data = data.reshape((nmodels, self.ntargets*2))
        combi = num.empty_like(data)
        combi[:, 0::2] = data[:, :self.ntargets]
        combi[:, 1::2] = data[:, self.ntargets:]
        misfits = combi.reshape((nmodels, self.ntargets, 2))

        mf_dict = {}
        for it, target in enumerate(self.targets):
            mf_dict[target.id] = misfits[:, it, 0]
        mf_dict['misfit_mean'] = num.mean(misfits, axis=1)[:, 0]

        mf_dict['niter'] = num.arange(nmodels, dtype=num.int)\
            + skip_nmodels + 1

        if keys is not None:
            for k in mf_dict.keys():
                if k not in keys:
                    mf_dict.pop(k)

        return nmodels, mf_dict


class SmoothColumnDataSource(ColumnDataSource):
    '''Gauss smoothed ColumnDataSource '''
    def __init__(self, *args, **kwargs):
        self._gaussian_kw = {
            'sigma': 3,
        }
        for kw in ['sigma', 'order', 'output', 'mode', 'cval', 'truncate']:
            if kw in kwargs.keys():
                self._gaussian_kw[kw] = kwargs.pop(kw)
        ColumnDataSource.__init__(self, *args, **kwargs)

    def stream(self, data, **kwargs):
        sigma = self._gaussian_kw['sigma']
        if data.keys() == self.data.keys():
            raise AttributeError('streaming data must represent'
                                 ' existing column')
        for key, arr in data.iteritems():
            if arr.ndim > 1 or not isinstance(arr, num.ndarray):
                raise AttributeError('data is not numpy.ndarray of 1d')

            nvalues = arr.size
            if nvalues < sigma:
                self.data[key][-(sigma-nvalues):]
                d = num.empty(sigma)
                d[:(sigma-nvalues)] = self.data[key][-(sigma-nvalues):]
                d[(sigma-nvalues):] = arr[:]
                arr = ndimage.gaussian_filter1d(d, **self._gaussian_kw)
                arr = arr[-nvalues:]
            else:
                arr = ndimage.gaussian_filter1d(arr, **self._gaussian_kw)
            data[key] = arr

        ColumnDataSource.stream(self, data, **kwargs)


class BaraddurRequestHandler(RequestHandler):
    def initialize(self, config):
        self.config = config
        self.model = config.model


class BaraddurBokehHandler(BokehHandler):
    def __init__(self, config, *args, **kwargs):
        BokehHandler.__init__(self, *args, **kwargs)
        self.config = config
        self.model = config.model
        self.nmodels = 0

        self.source = ColumnDataSource()

    def update_source(self, func):
        new_nmodels, new_data = func(
            skip_nmodels=self.nmodels,
            keys=self.source.data.keys())
        self.source.stream(new_data)
        self.nmodels += new_nmodels

    @gen.coroutine
    def update_models(self):
        return self.update_source(self.model.get_models)

    @gen.coroutine
    def update_misfits(self):
        return self.update_source(self.model.get_misfits)


class Rundirs(BaraddurRequestHandler):

    def get(self):
        models = [BaraddurModel(rundir=rd)
                  for rd in self.config.available_rundirs]
        self.render('rundirs.html',
                    pages=pages,
                    project_dir=self.config.project_dir,
                    models=models,
                    active_rundir=self.config.rundir)

    def post(self):
        rundir = self.get_argument('rundir', None)
        if rundir is not None:
            self.config.rundir = rundir
        self.get()


class Summary(BaraddurRequestHandler):

    @gen.coroutine
    def get(self):
        self.render('summary.html',
                    pages=pages,
                    problem=self.model.problem)


class Misfit(BaraddurRequestHandler):

    class MisfitsPlot(BaraddurBokehHandler):

        def modify_document(self, doc):
            self.source.add(num.ndarray(0, dtype=num.float32),
                            'misfit_mean')
            self.source.add(num.ndarray(0, dtype=num.float32),
                            'niter')
            self.update_misfits()

            plot = figure(webgl=True,
                          responsive=True,
                          plot_height=300,
                          x_axis_label='Iteration #',
                          y_axis_label='Misfit')
            plot.scatter('niter', 'misfit_mean',
                         source=self.source, alpha=.4)

            doc.add_root(plot)
            doc.add_periodic_callback(self.update_misfits,
                                      self.config.plot_update_small)

    bokeh_handlers = {'misfit': MisfitsPlot}

    @gen.coroutine
    def get(self):
        self.render('misfit.html',
                    pages=pages,
                    misfit_plot=autoload_server(
                        None,
                        app_path='/misfit',
                        url='plots'),
                    problem=self.model.problem)


class Parameters(BaraddurRequestHandler):

    class ParameterPlots(BaraddurBokehHandler):

        ncols = 4

        def modify_document(self, doc):
            problem = self.model.problem

            for param in ['niter'] + [p.name for p in problem.parameters]:
                self.source.add(
                    num.ndarray(0, dtype=num.float32),
                    param)
            self.model = BaraddurModel(self.config.rundir)
            self.update_models()

            plots = []
            for par in problem.parameters:
                fig = figure(webgl=True,
                             responsive=True,
                             x_axis_label='Iteration #',
                             y_axis_label='%s [%s]' % (par.label, par.unit))
                fig.scatter('niter', par.name,
                            source=self.source, alpha=.4)
                plots.append(fig)
            if (len(plots) % self.ncols) != 0:
                plots += [None] * (self.ncols - (len(plots) % self.ncols))

            grid = layouts.gridplot(
                plots,
                responsive=True,
                ncols=self.ncols)

            doc.add_root(grid)
            doc.add_periodic_callback(self.update_models,
                                      self.config.plot_update_small)

    bokeh_handlers = {'parameters': ParameterPlots}

    @gen.coroutine
    def get(self):

        self.render('parameters.html',
                    pages=pages,
                    parameter_plot=autoload_server(
                        None,
                        app_path='/parameters',
                        url='plots'),
                    problem=self.model.problem)


class Targets(BaraddurRequestHandler):

    class TargetContributionPlot(BaraddurBokehHandler):

        def __init__(self, *args, **kwargs):
            BaraddurBokehHandler.__init__(self, *args, **kwargs)
            self.source = SmoothColumnDataSource(sigma=20)

        def modify_document(self, doc):
            plot = figure(webgl=True,
                          responsive=True,
                          plot_height=300,
                          x_axis_label='Iteration #',
                          y_axis_label='Misfit')

            target_ids = [t.id for t in self.model.targets]
            for param in ['niter'] + target_ids:
                self.source.add(
                    num.ndarray(0, dtype=num.float32),
                    param)
            self.update_misfits()
            colors = ['blue', 'red']
            for it, target in enumerate(self.model.targets):
                plot.scatter('niter', target.id,
                             source=self.source,
                             color=colors[it])

            # plot.patches(
            #     ['niter'] * self.model.ntargets,
            #     target_ids,
            #     source=self.source)

            doc.add_root(plot)
            doc.add_periodic_callback(self.update_misfits,
                                      self.config.plot_update_small)

    bokeh_handlers = {'targets': TargetContributionPlot}

    @gen.coroutine
    def get(self):
        self.render('targets.html',
                    target_misfit_plot=autoload_server(
                        None,
                        app_path='/targets',
                        url='plots'),
                    pages=pages,
                    problem=self.model.problem)


pages = OrderedDict([
    ('Rundirs', Rundirs),
    ('Summary', Summary),
    ('Misfit', Misfit),
    ('Parameters', Parameters),
    ('Targets', Targets),
])


class BaraddurConfig(Object):
    rundir = String.T()
    project_dir = String.T(
        help='Grond project dir.',
        optional=True)
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
    plot_update_small = Float.T(
        default=1e3,
        optional=True,
        help='Update rate for small, performant plots. E.g. misfits')
    plot_update_big = Float.T(
        default=2.5e3,
        optional=True,
        help='Update rate for big, heavy plots. E.g. parameters.')

    def __init__(self, *args, **kwargs):
        Object.__init__(self, *args, **kwargs)
        self._rundir = None
        self._model = None

    @property
    def rundir(self):
        return self._rundir

    @rundir.getter
    def rundir(self):
        if self._rundir is None:
            if len(self.available_rundirs) > 0:
                self.rundir = self.available_rundirs[0]
        return self._rundir

    @rundir.setter
    def rundir(self, value):
        self._rundir = value
        self._model = BaraddurModel(rundir=self._rundir)

    @property
    def model(self):
        if self._model is None:
            self._model = BaraddurModel(rundir=self.rundir)
        return self._model

    @property
    def available_rundirs(self):
        rundirs = []
        for d in [d for d in glob.glob(op.join(self.project_dir, '*'))
                  if op.isdir(d)]:
            if BaraddurModel.validate_rundir(d):
                rundirs.append(d)
        return rundirs


class Baraddur(BokehServer):

    def __init__(self, config=None, project_dir=None, *args, **kwargs):
        if config is not None:
            self.config = config
        elif project_dir is not None:
            self.config = BaraddurConfig(project_dir=project_dir)
        else:
            raise AttributeError('No project dir set!')
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
                    logger.warning('Port %d in use, bailing to %d'
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
                bokeh_apps['/plots/%s' % url] = Application(bokeh_handler(
                    self.config))
        return bokeh_apps

    def get_tornado_handlers(self):
        return [(r'/', pages.values()[0],
                 {'config': self.config})] +\
               [(r'/%s' % title, handler,
                 {'config': self.config})
                for title, handler in pages.iteritems()] +\
               [(r'/res/(.*)', StaticFileHandler,
                {'path': op.join(op.dirname(__file__), 'res')})]

    def start(self, signal=None):
        logger.info('Starting Baraddur server on http://%s:%d'
                    % ('localhost', self.port))

        if signal is not None:
            def shutdown():
                if not signal.empty():
                    self.stop()
            tornado.ioloop.PeriodicCallback(shutdown, 2000).start()

        BokehServer.start(self)
        self.ioloop.start()

    def stop(self, *args, **kwargs):
        logger.info('Stopping Baraddur server...')
        BokehServer.stop(self)
        self.ioloop.stop()


if __name__ == '__main__':
    baraddur = Baraddur(
        rundir='/home/marius/Development/testing/grond/rundir')
    baraddur.start()
