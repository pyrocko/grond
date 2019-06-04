from __future__ import print_function, absolute_import

import time
import logging
import threading
import os.path as op
import numpy as num
from datetime import timedelta

from pyrocko import util, guts
from grond.environment import Environment


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


class TerminalMonitor(object):

    def __init__(self, nlines_scroll=10):
        self.nlines_scroll = nlines_scroll

    def __enter__(self):
        print('\033[2J\033[1;%ir\033[%i;1H'
              % (self.nlines_scroll, self.nlines_scroll-1), end=None)
        return self

    def __exit__(self, type, value, traceback):
        print('\033[r\033[%i;1H\033[0J\033[%i;1H'
              % (self.nlines_scroll+1, self.nlines_scroll-1))

    def _start_show(self):
        print('\033[%i;1H\033[0J' % (self.nlines_scroll+1), end=None)

    def _end_show(self):
        print('\033[%i;1H' % (self.nlines_scroll-1), end=None)

    def show(self, s):
        self._start_show()
        print(s, end=None)
        self._end_show()


class GrondMonitor(threading.Thread):

    col_width = 15
    row_name = color.BOLD + '{:<{col_param_width}s}' + color.END
    parameter_fmt = '{:{col_width}s}'

    def __init__(self, rundir):
        threading.Thread.__init__(self)
        self.rundir = rundir
        self.daemon = True

        self.sig_terminate = threading.Event()
        self.iter_per_second = 0
        self._iiter = 0
        self._iter_buffer = RingBuffer(20)
        self._tm = None

    def run(self):
        logger.info('Waiting to follow environment %s...' % self.rundir)
        env = Environment.discover(self.rundir)
        if env is None:
            logger.error('Could not attach to Grond environment.')
            return

        self.environment = env
        self.history = self.environment.get_history()

        optimiser_fn = op.join(self.rundir, 'optimiser.yaml')
        self.optimiser = guts.load(filename=optimiser_fn)

        self.problem = self.history.problem
        self.niter = self.optimiser.niterations

        self.starttime = time.time()
        self.last_update = self.starttime

        self.history.add_listener(self)

        with TerminalMonitor(10) as tm:

            self._tm = tm

            ii = 0
            while True:
                ii += 1
                self.history.update()
                time.sleep(0.1)
                if self.sig_terminate.is_set():
                    break

        logger.debug('Monitor thread exiting.')

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
        optimiser_status = self.optimiser.get_status(self.history)
        row_names = optimiser_status.row_names

        lines = []
        lnadd = lines.append

        def fmt(s):
            return util.gform(s, significant_digits=(self.col_width-1-6)//2)

        lnadd('Problem:   {p.name}'
              '\t({s.runtime} - remaining {s.runtime_remaining}'
              ' @ {s.iter_per_second:.1f} iter/s)'
              .format(s=self, p=problem))
        lnadd('Iteration: {s.iiter} / {s.niter}'
              .format(s=self))
        if optimiser_status.extra_header is not None:
            lnadd(optimiser_status.extra_header)

        col_param_width = max([len(p) for p in row_names]) + 2

        out_ln = self.row_name +\
            ''.join([self.parameter_fmt] * optimiser_status.ncolumns)

        lnadd(out_ln.format(
            *['Parameter'] + list(optimiser_status.column_names),
            col_param_width=col_param_width,
            col_width=self.col_width,
            type='s'))

        for ip, parameter_name in enumerate(row_names):
            lnadd(out_ln.format(
                 parameter_name,
                 *[fmt(v[ip]) for v in optimiser_status.values],
                 col_param_width=col_param_width,
                 col_width=self.col_width))

        if optimiser_status.extra_footer is not None:
            lnadd(optimiser_status.extra_footer)

        self._tm.show('\n'.join(lines))

    def terminate(self):
        logger.debug('Setting thread termination flag.')
        self.sig_terminate.set()
        self.join()

    @classmethod
    def watch(cls, rundir):
        monitor = cls(rundir)
        monitor.start()
        return monitor
