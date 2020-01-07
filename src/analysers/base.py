from pyrocko.guts import Object

from grond.meta import GrondError

guts_prefix = 'grond'


class Analyser(object):

    def analyse(self, problem, ds):
        pass

    def get_rstate(self, problem, seed=None):
        return problem.get_rstate_manager().get_rstate(
            self.__class__.__name__, seed)


class AnalyserConfig(Object):

    def get_analyser(self):
        return Analyser


class AnalyserResult(Object):

    class NoResults(GrondError):
        pass

    pass


__all__ = '''
    Analyser
    AnalyserConfig
'''.split()
