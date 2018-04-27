import logging

from pyrocko.guts import Object
from grond.meta import GrondError

guts_prefix = 'grond'

logger = logging.getLogger('grond.optimisers.base')


class BadProblem(GrondError):
    pass


class Optimiser(Object):

    def optimise(self, problem):
        raise NotImplemented()

    @property
    def niterations(self):
        raise NotImplementedError()

    def get_status(self, history):
        pass

    @classmethod
    def get_plot_classes(cls):
        from . import plot
        return plot.get_plot_classes()


class OptimiserConfig(Object):
    pass


class OptimiserStatus(object):
    __slots__ = ['row_names', 'columns', 'extra_text']

    def __init__(self, row_names, columns, extra_text):
        self.row_names = row_names
        self.columns = columns
        self.extra_text = extra_text

    @property
    def column_names(self):
        return self.columns.keys()

    @property
    def ncolumns(self):
        return len(self.columns)

    @property
    def values(self):
        return self.columns.values()


__all__ = '''
    BadProblem
    Optimiser
    OptimiserConfig
'''.split()
