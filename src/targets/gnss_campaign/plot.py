from pyrocko.guts import Tuple, Float, Bool
from grond.plot.config import PlotConfig, PlotItem
from pyrocko.plot import automap

guts_prefix = 'grond'


class GNSSTargetMisfitPlot(PlotConfig):

    name = 'fits_gnss'
    size_cm = Tuple.T(
        2, Float.T(),
        default=(20., 20.),
        help='size of the figure in cm')
    show_topo = Bool.T(
        default=True,
        help='show topography')
    show_grid = Bool.T(
        default=True,
        help='show the lat/lon grid')
    show_rivers = Bool.T(
        default=True,
        help='show rivers on the map')

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history()
        optimiser = environ.get_optimiser()
        ds = environ.get_dataset()

        cm.create_group_automap(
            self,
            self.draw_gnss_fits(ds, history, optimiser))

    def plot_gnss_fits(self, ds, history, optimiser):
        problem = history.problem

        for itarget, gnss_target in enumerate(problem.gnss_targets):
            camp = gnss_target.campaign
            item = PlotItem(
                name='fig_%i' % itarget,
                attributes={
                    'targets': gnss_target.path
                })

            lat, lon = camp.get_center_latlon()
            radius = camp.get_radius()

            m = automap.Map(
                width=self.size[0],
                height=self.size[1],
                lat=lat,
                lon=lon,
                radius=radius,
                show_topo=self.show_topo,
                show_grid=self.show_grid,
                show_rivers=self.show_rivers,
                color_wet=(216, 242, 254),
                color_dry=(238, 236, 230)
                )

            m.add_gnss_campaign(camp)
            yield item, m
