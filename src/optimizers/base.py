import logging
import time
import numpy as num
from datetime import timedelta
from collections import OrderedDict

from pyrocko.guts import Object

guts_prefix = 'grond'

logger = logging.getLogger('grond.solver')


class SimpleTimedelta(timedelta):
    def __str__(self):
        return timedelta.__str__(self).split('.')[0]


class RingBuffer(num.ndarray):
    def __init__(self, *args, **kwargs):
        num.ndarray.__init__(self, *args, **kwargs)
        self.fill(0.)
        self.pos = 0

    def put(self, value):
        self[self.pos] = value
        self.pos += 1
        self.pos %= self.size


class SolverState(object):
    problem_name = ''
    parameter_sets = OrderedDict()
    parameter_names = []

    starttime = time.time()

    niter = 0
    iter_per_second = 0.

    extra_text = ''

    _iiter = 0
    _iter_buffer = RingBuffer(25)
    _last_update = time.time()

    @property
    def iiter(self):
        return self._iiter

    @iiter.setter
    def iiter(self, value):
        dt = time.time() - self._last_update
        self._iter_buffer.put(float((value - self._iiter) / dt))
        self.iter_per_second = float(self._iter_buffer.mean())
        self._iiter = value
        self._last_update = time.time()

    @property
    def runtime(self):
        return SimpleTimedelta(seconds=time.time() - self.starttime)

    @property
    def runtime_remaining(self):
        if self.iter_per_second == 0.:
            return SimpleTimedelta()
        return SimpleTimedelta(seconds=(self.niter - self.iiter)
                               / self.iter_per_second)


class Solver(object):
    state = SolverState()

    def solve(
            self, problem, rundir=None, status=(), plot=None, xs_inject=None,
            notifier=None):
        raise NotImplemented()


class SolverConfig(Object):

    def get_solver(self):
        return Solver()


__all__ = '''
    Solver
    SolverState
    SolverConfig
'''.split()
