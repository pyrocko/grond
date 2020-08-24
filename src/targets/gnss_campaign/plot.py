import logging
import numpy as num
from pyrocko.model import gnss

from pyrocko.plot import automap, mpl_init
from pyrocko import orthodrome as od

from grond.plot.config import PlotConfig
from grond.plot.collection import PlotItem
from grond.problems import CMTProblem, RectangularProblem, \
    VLVDProblem

from ..plot import StationDistributionPlot

import copy
from pyrocko.guts import Tuple, Float, Bool

guts_prefix = 'grond'
km = 1e3

logger = logging.getLogger('grond.targets.gnss_campaign.plot')


class GNSSTargetMisfitPlot(PlotConfig):
    ''' Maps showing horizontal surface displacements
        of a GNSS campaign and model '''

    name = 'gnss'

    size_cm = Tuple.T(
        2, Float.T(),
        default=(30., 30.),
        help='width and length of the figure in cm')
    show_topo = Bool.T(
        default=False,
        help='show topography')
    show_grid = Bool.T(
        default=True,
        help='show the lat/lon grid')
    show_rivers = Bool.T(
        default=True,
        help='show rivers on the map')
    radius = Float.T(
        optional=True,
        help='radius of the map around campaign center lat/lon')

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history(subset='harvest')
        optimiser = environ.get_optimiser()
        ds = environ.get_dataset()

        environ.setup_modelling()

        cm.create_group_automap(
            self,
            self.draw_gnss_fits(ds, history, optimiser),
            title=u'GNSS Displacements',
            section='fits',
            feather_icon='map',
            description=u'''
Maps showing station positions and statiom names of the GNSS targets.

Arrows the observed surface displacements (black arrows) and synthetic
displacements (red arrows). The top plot shows the horizontal displacements and
the bottom plot the vertical displacements. The grey filled box shows the
surface projection of the modelled source, with the thick-lined edge marking
the upper fault edge.
''')

    def draw_gnss_fits(self, ds, history, optimiser, vertical=False):
        problem = history.problem

        gnss_targets = problem.gnss_targets
        for target in gnss_targets:
            target.set_dataset(ds)

        xbest = history.get_best_model()
        source = history.get_best_source()

        results = problem.evaluate(
            xbest, result_mode='full', targets=gnss_targets)

        def plot_gnss(gnss_target, result, ifig, vertical=False):
            campaign = gnss_target.campaign
            item = PlotItem(
                name='fig_%i' % ifig,
                attributes={
                    'targets': gnss_target.path
                },
                title=u'Static GNSS Surface Displacements - Campaign %s'
                      % campaign.name,
                description=u'''
Static surface displacement from GNSS campaign %s (black vectors) and
displacements derived from best model (red).
''' % campaign.name)

            event = source.pyrocko_event()
            locations = campaign.stations + [event]

            lat, lon = od.geographic_midpoint_locations(locations)

            if self.radius is None:
                coords = num.array([loc.effective_latlon for loc in locations])
                radius = od.distance_accurate50m_numpy(
                            lat[num.newaxis], lon[num.newaxis],
                            coords[:, 0].max(), coords[:, 1]).max()
                radius *= 1.1

            if radius < 30.*km:
                logger.warn(
                    'Radius of GNSS campaign %s too small, defaulting'
                    ' to 30 km' % campaign.name)
                radius = 30*km

            model_camp = gnss.GNSSCampaign(
                stations=copy.deepcopy(campaign.stations),
                name='grond model')
            for ista, sta in enumerate(model_camp.stations):
                sta.north.shift = result.statics_syn['displacement.n'][ista]
                sta.north.sigma = 0.

                sta.east.shift = result.statics_syn['displacement.e'][ista]
                sta.east.sigma = 0.

                if sta.up:
                    sta.up.shift = -result.statics_syn['displacement.d'][ista]
                    sta.up.sigma = 0.

            m = automap.Map(
                width=self.size_cm[0],
                height=self.size_cm[1],
                lat=lat,
                lon=lon,
                radius=radius,
                show_topo=self.show_topo,
                show_grid=self.show_grid,
                show_rivers=self.show_rivers,
                color_wet=(216, 242, 254),
                color_dry=(238, 236, 230))

            all_stations = campaign.stations + model_camp.stations
            offset_scale = num.zeros(len(all_stations))

            for ista, sta in enumerate(all_stations):
                for comp in sta.components.values():
                    offset_scale[ista] += comp.shift
            offset_scale = num.sqrt(offset_scale**2).max()

            m.add_gnss_campaign(
                campaign,
                psxy_style={
                    'G': 'black',
                    'W': '0.8p,black',
                },
                offset_scale=offset_scale,
                vertical=vertical)

            m.add_gnss_campaign(
                model_camp,
                psxy_style={
                    'G': 'red',
                    'W': '0.8p,red',
                    't': 30,
                },
                offset_scale=offset_scale,
                vertical=vertical,
                labels=False)

            if isinstance(problem, CMTProblem):
                from pyrocko import moment_tensor
                from pyrocko.plot import gmtpy

                mt = event.moment_tensor.m_up_south_east()
                ev_lat, ev_lon = event.effective_latlon

                xx = num.trace(mt) / 3.
                mc = num.matrix([[xx, 0., 0.], [0., xx, 0.], [0., 0., xx]])
                mc = mt - mc
                mc = mc / event.moment_tensor.scalar_moment() * \
                    moment_tensor.magnitude_to_moment(5.0)
                m6 = tuple(moment_tensor.to6(mc))
                symbol_size = 20.
                m.gmt.psmeca(
                    S='%s%g' % ('d', symbol_size / gmtpy.cm),
                    in_rows=[(ev_lon, ev_lat, 10) + m6 + (1, 0, 0)],
                    M=True,
                    *m.jxyr)

            elif isinstance(problem, RectangularProblem):
                m.gmt.psxy(
                    in_rows=source.outline(cs='lonlat'),
                    L='+p2p,black',
                    W='1p,black',
                    G='black',
                    t=60,
                    *m.jxyr)

            elif isinstance(problem, VLVDProblem):
                ev_lat, ev_lon = event.effective_latlon
                dV = abs(source.volume_change)
                sphere_radius = num.cbrt(dV / (4./3.*num.pi))

                volcanic_circle = [
                    ev_lon,
                    ev_lat,
                    '%fe' % sphere_radius
                ]
                m.gmt.psxy(
                    S='E-',
                    in_rows=[volcanic_circle],
                    W='1p,black',
                    G='orange3',
                    *m.jxyr)

            return (item, m)

        ifig = 0
        for vertical in (False, True):
            for gnss_target, result in zip(problem.gnss_targets, results):
                yield plot_gnss(gnss_target, result, ifig, vertical)
                ifig += 1


class GNSSStationDistribution(StationDistributionPlot):
    ''' Polar plot showing GNSS station distribution and weight '''
    name = 'gnss_station_distribution'

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        mpl_init(fontsize=self.font_size)

        history = environ.get_history(subset='harvest')
        problem = environ.get_problem()
        dataset = environ.get_dataset()

        cm.create_group_mpl(
            self,
            self.draw_figures(problem, dataset, history),
            title=u'GNSS Station Distribution',
            section='checks',
            feather_icon='target',
            description=u'''
Plots showing the GNSS station distribution and their weight.

This polar plot visualises the station distribution in distance and azimuth,
the marker's size is scaled to the stations weight (mean of spatial
components).
''')

    def draw_figures(self, problem, dataset, history):

        event = problem.base_source
        targets = problem.gnss_targets

        for target in targets:
            target.set_dataset(dataset)
            comp_weights = target.component_weights()[0]

            ws_n = comp_weights[:, 0::3] / comp_weights.max()
            ws_e = comp_weights[:, 1::3] / comp_weights.max()
            ws_u = comp_weights[:, 2::3] / comp_weights.max()
            ws_e = num.array(ws_e[0]).flatten()
            ws_n = num.array(ws_n[0]).flatten()
            ws_u = num.array(ws_u[0]).flatten()

            if ws_n.size == 0:
                continue

            distances = target.distance_to(event)
            azimuths = od.azibazi_numpy(
                num.array(event.effective_lat)[num.newaxis],
                num.array(event.effective_lon)[num.newaxis],
                target.get_latlon()[:, 0],
                target.get_latlon()[:, 1])[0]
            labels = target.station_names

            item = PlotItem(name='station_distribution-N-%s' % target.path)
            fig, ax, legend = self.plot_station_distribution(
                azimuths, distances, ws_n, labels)
            legend.set_title('Weight, N components')

            yield (item, fig)

            item = PlotItem(name='station_distribution-E-%s' % target.path)
            fig, ax, legend = self.plot_station_distribution(
                azimuths, distances, ws_e, labels)
            legend.set_title('Weight, E components')

            yield (item, fig)

            item = PlotItem(name='station_distribution-U-%s' % target.path)
            fig, ax, legend = self.plot_station_distribution(
                azimuths, distances, ws_u, labels)
            legend.set_title('Weight, U components')

            yield (item, fig)


def get_plot_classes():
    return [
        GNSSTargetMisfitPlot,
        GNSSStationDistribution
    ]
