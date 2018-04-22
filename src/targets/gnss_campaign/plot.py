from pyrocko.guts import Tuple, Float
from grond.plot.config import PlotConfig


class GNSSTargetMisfitPlot(PlotConfig):

    name = 'fits_gnss'
    size_cm = Tuple.T(2, Float.T(), default=(20., 20.))

    def make(self, environ):
        return

        cm = environ.get_plot_collection_manager()
        history = environ.get_history()
        optimizer = environ.get_optimizer()
        ds = environ.get_dataset()

        cm.create_group_gmt(
            self,
            self.draw_static_fits(ds, history, optimizer))

    def plot_misfit(self):
        return self
