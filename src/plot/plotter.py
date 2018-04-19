
class Plotter(object):

    def __init__(self):
        raise NotImplementedError(
            'Plotter should not be instantiated, use its classmethods')

    @classmethod
    def draw_check_figures(cls, sources, target, results):
        raise NotImplementedError('to be implemented in subclass')

    @staticmethod
    def register(section, name):
        def decorate(meth):
            Plotter.plots[section][name].append(meth)
            return meth

        return decorate

    plots = defaultdict(lambda: defaultdict(list))


class NoPlotterClassAvailable(Exception):
    pass
