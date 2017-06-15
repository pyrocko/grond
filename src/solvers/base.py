import logging

from pyrocko.guts import Object

guts_prefix = 'grond'

logger = logging.getLogger('grond.solver')


class Solver(object):
    def solve(
            self, problem, rundir=None, status=(), plot=None, xs_inject=None):

        raise NotImplemented()


class SolverConfig(Object):

    def get_solver(self):
        return Solver()


__all__ = '''
    Solver
    SolverConfig
'''.split()
