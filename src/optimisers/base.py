import logging

from pyrocko.guts import Object
from grond.meta import GrondError, has_get_plot_classes

guts_prefix = 'grond'

logger = logging.getLogger('grond.optimisers.base')


class BadProblem(GrondError):
    pass


@has_get_plot_classes
class Optimiser(Object):

    def __init__(self, **kwargs):
        Object.__init__(self, **kwargs)
        self._nthreads = 1

    def set_nthreads(self, nthreads):
        logger.debug('Setting nthreads to %d', nthreads)
        self._nthreads = nthreads

    def optimise(self, problem):
        raise NotImplementedError

    @property
    def niterations(self):
        raise NotImplementedError

    def get_status(self, history):
        pass

    def init_bootstraps(self, problem):
        raise NotImplementedError

    @classmethod
    def get_plot_classes(cls):
        from . import plot
        return plot.get_plot_classes()


class OptimiserConfig(Object):
    pass


class OptimiserStatus(object):
    __slots__ = ['row_names', 'column_data', 'extra_header', 'extra_footer']

    def __init__(self, row_names, column_data,
                 extra_header=None, extra_footer=None):
        self.row_names = row_names
        self.column_data = column_data
        self.extra_header = extra_header
        self.extra_footer = extra_footer

    @property
    def column_names(self):
        return self.column_data.keys()

    @property
    def ncolumns(self):
        return len(self.column_data)

    @property
    def values(self):
        return self.column_data.values()


__all__ = '''
    BadProblem
    Optimiser
    OptimiserConfig
'''.split()
