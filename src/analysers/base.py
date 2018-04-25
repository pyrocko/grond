from pyrocko.guts import Object, Bool

guts_prefix = 'grond'


class Analyser(object):

    def analyse(self, problem, ds):
        pass


class AnalyserConfig(Object):
        
    def get_analyser(self):
        return Analyser


__all__ = '''
    Analyser
    AnalyserConfig
'''.split()
