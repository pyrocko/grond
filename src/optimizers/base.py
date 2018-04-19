import logging

from pyrocko.guts import Object
from grond.meta import GrondError

guts_prefix = 'grond'

logger = logging.getLogger('grond.optimizers.base')


class BadProblem(GrondError):
    pass


class Optimizer(Object):

    def optimize(self, problem):
        raise NotImplemented()

    @property
    def niterations(self):
        raise NotImplementedError()

    def get_status(self, problem):
        pass


class OptimizerConfig(Object):
    pass


__all__ = '''
    BadProblem
    Optimizer
    OptimizerConfig
'''.split()
