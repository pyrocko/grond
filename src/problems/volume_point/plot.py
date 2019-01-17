from ..plot import MTLocationPlot
from pyrocko.plot import mpl_init


class VolumePointLocationPlot(MTLocationPlot):
    name = 'location_volume'
    beachball_type = 'full'

    def make(self, environ):
        environ.setup_modelling()
        cm = environ.get_plot_collection_manager()
        history = environ.get_history(subset='harvest')
        mpl_init(fontsize=self.font_size)
        self._to_be_closed = []
        cm.create_group_mpl(
            self,
            self.draw_figures(history),
            title=u'Volume Location',
            section='solution',
            feather_icon='target',
            description=u'''
Location plot of the ensemble of best solutions in three cross-sections.

The coordinate range is defined by the search space given in the config file.
Symbols show best volume locations, and colors indicate low (red) and
high (blue) misfit.
''')
        for obj in self._to_be_closed:
            obj.close()
