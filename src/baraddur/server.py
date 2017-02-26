import tornado.ioloop
import grond
import os.path as op
import numpy as np

from collections import OrderedDict

from tornado.web import RequestHandler, Application
from tornado import gen

from bokeh.application import Application as BokehApplication
from bokeh.embed import autoload_server
from bokeh.server.server import Server

from bokeh.layouts import column
from bokeh.models import ColumnDataSource, Slider
from bokeh.application.handlers import FunctionHandler
from bokeh.plotting import figure

problem = grond.core.load_problem_info('.')


class Summary(RequestHandler):
    from bokeh.models import map_plots

    @gen.coroutine
    def get(self):
        p = figure(plot_width=400, plot_height=400)
        p.circle([1, 2, 3, 4, 5], [6, 7, 2, 4, 5], size=20, color="navy",
                 alpha=0.5)

        script = autoload_server(p, url='http://localhost:5001')
        self.render('summary.html',
                    pages=pages,
                    problem=problem,
                    script=script)


class MainHandler(RequestHandler):
    @gen.coroutine
    def get(self):
        self.render('index.html',
                    pages=pages)


pages = OrderedDict([
    ('Summary', Summary),
    ('Sequences', MainHandler),
    ('Tradeoffs', MainHandler),
])

if __name__ == '__main__':
    bokeh_app = BokehApplication()
    bokeh_server = Server({'/bokeh': bokeh_app},
                          ioloop=tornado.ioloop.IOLoop.instance(),
                          num_procs=1)

    tornado_app = Application(
        [(r'/', Summary)] +
        [(r'/%s' % title, handler) for title, handler in pages.iteritems()],
        debug=False,
        autoreload=False,
        compiled_template_cache=False,
        static_path=op.join(op.dirname(__file__), 'static'),
        template_path=op.join(op.dirname(__file__), 'tpl'),
        )
    tornado_app.listen(8888)

    tornado.ioloop.IOLoop.current().start()
