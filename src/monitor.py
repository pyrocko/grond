import time
import logging
import os.path as op
import numpy as num
from datetime import timedelta

from pyrocko import util, guts
from grond.problems import ModelHistory


logger = logging.getLogger('grond.monit')


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


class GrondMonitor(object):

    col_width = 15
    row_name = color.BOLD + '{:<{col_param_width}s}' + color.END
    parameter_fmt = '{:{col_width}s}'

    def __init__(self, rundir):
        self.rundir = rundir

        self.iter_per_second = 0
        self._iiter = 0
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
        return timedelta(seconds=round(time.time() - self.starttime))

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
        return timedelta(seconds=round((self.niter - self.iiter)
                         / self.iter_per_second))

    def extend(self, *args):
        ''' Connected and called through the self.history.add_listener '''
        self.iiter = self.history.nmodels
        problem = self.history.problem
        optimizer_status = self.optimizer.get_status(self.history)

        lines = []
        lnadd = lines.append

        parameter_names = [p.name_nogroups for p in problem.parameters]

        def fmt(s):
            return util.gform(s, significant_digits=(self.col_width-1-6)//2)

        lnadd('Problem name: {p.name}'
              '\t({s.runtime} - remaining {s.runtime_remaining}'
              ' @ {s.iter_per_second:.1f} iter/s)'
              .format(s=self, p=problem))
        lnadd('Iteration {s.iiter} / {s.niter}'
              .format(s=self))
        lnadd(optimizer_status.extra_text)

        col_param_width = max([len(p) for p in parameter_names]) + 2

        out_ln = self.row_name +\
            ''.join([self.parameter_fmt] * optimizer_status.ncolumns)

        lnadd(out_ln.format(
            *['Parameter'] + list(optimizer_status.column_names),
            col_param_width=col_param_width,
            col_width=self.col_width,
            type='s'))

        for ip, parameter_name in enumerate(parameter_names):
            lnadd(out_ln.format(
                 parameter_name,
                 *[fmt(v[ip]) for v in optimizer_status.values],
                 col_param_width=col_param_width,
                 col_width=self.col_width))

        lines[0:0] = ['\033[2J']
        lnadd('')
        print('\n'.join(lines))

    @classmethod
    def watch(cls, rundir):
        import threading
        import signal

        monitor = cls(rundir)

        monitor_thread = threading.Thread(target=monitor.run)
        monitor_thread.start()

        return monitor_thread
