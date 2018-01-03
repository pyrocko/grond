import time
import logging
import os.path as op
import numpy as num
from datetime import timedelta

from pyrocko import util, guts
from grond.problems import ModelHistory
from .base import Listener


logger = logging.getLogger('TerminalListener')


class RingBuffer(num.ndarray):
    def __new__(cls, *args, **kwargs):
        cls = num.ndarray.__new__(cls, *args, **kwargs)
        cls.fill(0.)
        return cls

    def __init__(self, *args, **kwargs):
        self.pos = 0

    def put(self, value):
        self[self.pos] = value
        self.pos += 1
        self.pos %= self.size


class color:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


class TerminalListener(Listener):

    col_width = 15
    row_name = color.BOLD + '{:<{col_param_width}s}' + color.END
    parameter_fmt = '{:{col_width}s}'

    def __init__(self, rundir):
        Listener.__init__(self)
        self.rundir = rundir
        self._iiter = 0
        self.iter_per_second = 0
        self._iter_buffer = RingBuffer(20)

    def run(self):
        logger.info('Waiting to follow %s' % self.rundir)

        self.history = ModelHistory.follow(self.rundir)

        optimizer_fn = op.join(self.rundir, 'optimizer.yaml')
        self.optimizer = guts.load(filename=optimizer_fn)

        self.problem = self.history.problem
        self.niter = self.optimizer.niterations

        self.starttime = time.time()
        self.last_update = self.starttime

        self.history.add_listener(self)
        self.start_watch()

    def start_watch(self):
        while True:
            self.history.update()
            time.sleep(.1)

    @property
    def runtime(self):
        return timedelta(seconds=time.time() - self.starttime)

    @property
    def iiter(self):
        return self._iiter

    @iiter.setter
    def iiter(self, iiter):
        dt = time.time() - self.last_update
        self._iter_buffer.put(float((iiter - self.iiter) / dt))
        self.iter_per_second = float(self._iter_buffer.mean())
        self._iiter = iiter
        self.last_update = time.time()

    @property
    def runtime_remaining(self):
        if self.iter_per_second == 0.:
            return timedelta()
        return timedelta(seconds=(self.niter - self.iiter)
                         / self.iter_per_second)

    def extend(self, *args):
        self.iiter = self.history.nmodels
        problem = self.history.problem

        lines = []
        ladd = lines.append

        def fmt(s):
            return util.gform(s, significant_digits=(self.col_width-1-6)//2)

        ladd('Problem name: {p.name}'
             '\t({s.runtime} - remaining {s.runtime_remaining}'
             ' @ {s.iter_per_second:.1f} iter/s)'
             .format(s=self, p=problem))
        ladd('Iteration {s.iiter} / {s.niter}'
             .format(s=self))

        out_ln = self.row_name +\
            ''.join([self.parameter_fmt] * len(problem.parameter_names))
        col_param_width = max([len(p) for p in problem.parameter_names]) + 2

        ladd(out_ln.format(
            *['Parameter'] + list(problem.parameter_names),
            col_param_width=col_param_width,
            col_width=self.col_width,
            type='s'))

        # for ip, parameter_name in enumerate(problem.parameter_names):
        #     ladd(out_ln.format(
        #          parameter_name,
        #          *[fmt(v[ip]) for v in problem.parameter_sets.values()],
        #          col_param_width=col_param_width,
        #          col_width=self.col_width))

        # ladd(problem.extra_text.format(
        #     col_param_width=col_param_width,
        #     col_width=self.col_width,))

        lines[0:0] = ['\033[2J']
        ladd('')
        print('\n'.join(lines))
