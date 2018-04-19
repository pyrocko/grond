

class Plotter(object):
    def __init__(self):
        raise NotImplementedError(
            'Plotter should not be instantiated, use its classmethods')


class NoPlotterClassAvailable(Exception):
    pass


__all__ = [
    'Plotter',
    'NoPlotterClassAvailable',
]
