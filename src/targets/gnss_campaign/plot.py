from pyrocko.guts import Tuple, Float
from grond.plot.config import PlotConfig
from pyrocko.plot import automap

guts_prefix = 'grond'


class GNSSTargetMisfitPlot(PlotConfig):

    name = 'fits_gnss'
    size_cm = Tuple.T(2, Float.T(), default=(20., 20.))

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history()
        optimizer = environ.get_optimizer()
        ds = environ.get_dataset()

        cm.create_group_pygmt(
            self,
            self.draw_gnss_fits(ds, history, optimizer))

    def plot_gnss_fits(self, ds, history, optimizer):
        lat, lon = self.get_center_latlon()
        radius = self.get_radius()

        m = automap.Map(
            width=30.,
            height=30.,
            lat=lat,
            lon=lon,
            radius=radius,
            show_topo=True,
            show_grid=True,
            show_rivers=True,
            color_wet=(216, 242, 254),
            color_dry=(238, 236, 230)
            )
