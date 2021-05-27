import logging

import numpy as num

from pyrocko import orthodrome as pod
from pyrocko.guts import Float, Bool, Tuple

from pyrocko.plot import automap

from grond.plot.collection import PlotItem
from grond.plot.config import PlotConfig

logger = logging.getLogger('grond.problem.double_sf.plot')

guts_prefix = 'grond'

km = 1e3


class SFForcePlot(PlotConfig):
    ''' Maps showing horizontal and vertical force
        of the best double single force model '''

    name = 'forces_double_singleforce'

    size_cm = Tuple.T(
        2, Float.T(),
        default=(15., 15.),
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
            self.draw_best_sf(ds, history, optimiser),
            title=u'Single Force Source Forces',
            section='solution',
            feather_icon='map',
            description=u'''
Maps show located force vectors of the best double Single Force Source model.

Arrows show the modelled forces (red arrows). The top plot shows the horizontal
forces and the bottom plot the vertical force. The dot indicates the location
of the best double single force source model.
''')

    def draw_best_sf(self, ds, history, optimiser, vertical=False):
        from grond.core import make_stats

        source = history.get_best_source()

        problem = history.problem
        models = history.models

        stats = make_stats(
            problem, models, history.get_primary_chain_misfits())

        def plot_double_sf(source, stats, ifig, vertical=False):
            orient = 'vertical' if vertical else 'horizontal'

            item = PlotItem(
                name='fig_%i' % ifig,
                attributes={},
                title=u'Best %s double single force model force vector' % (
                    orient),
                description=u'''
Double single force source %s force vector for the best model (red). The circle
shows the 95%% confidence ellipse.
''' % orient)

            event = source.pyrocko_event()

            radius = self.radius
            if radius is None or radius < 30.*km:
                logger.warn(
                    'Radius too small, defaulting to 30 km')
                radius = 30*km

            m = automap.Map(
                width=self.size_cm[0],
                height=self.size_cm[1],
                lat=event.lat,
                lon=event.lon,
                radius=radius,
                show_topo=self.show_topo,
                show_grid=self.show_grid,
                show_rivers=self.show_rivers,
                color_wet=(216, 242, 254),
                color_dry=(238, 236, 230))

            offset_scale = source.force
            size = num.linalg.norm(self.size_cm)

            scale = (size / 5.) / offset_scale

            lat, lon = pod.ne_to_latlon(
                event.lat,
                event.lon,
                source.north_shift,
                source.east_shift)

            source.lat, source.lon = lat, lon
            sf1, sf2 = source.split()

            stats_dict = stats.get_values_dict()

            if vertical:
                rows = [[
                    sf1.effective_lon, sf1.effective_lat,
                    0., -sf1.fd * scale,
                    (stats_dict['rfn1.std'] + stats_dict['rfe1.std']),
                    stats_dict['rfd1.std'],
                    0.]]

                rows.append([
                    sf2.effective_lon, sf2.effective_lat,
                    0., -sf2.fd * scale,
                    (stats_dict['rfn2.std'] + stats_dict['rfe2.std']),
                    stats_dict['rfd2.std'],
                    0.])

            else:
                rows = [[
                    sf1.effective_lon, sf1.effective_lat,
                    sf1.fe * scale, sf1.fn * scale,
                    stats_dict['rfe1.std'],
                    stats_dict['rfn1.std'],
                    0.]]

                rows.append([
                    sf2.effective_lon, sf2.effective_lat,
                    sf2.fe * scale, sf2.fn * scale,
                    stats_dict['rfe2.std'],
                    stats_dict['rfn2.std'],
                    0.])

            fontsize = 10.

            default_psxy_style = {
                'h': 0,
                'W': '2.0p,red',
                'A': '+p4p,black+e+a40',
                'G': 'red',
                't': 30,
                'L': True,
                'S': 'e1c/0.95/%d' % fontsize,
            }

            m.gmt.psvelo(
                in_rows=rows,
                *m.jxyr,
                **default_psxy_style)

            m.gmt.psxy(
                S='c10p',
                in_rows=[[lon, lat]],
                W='1p,black',
                G='orange3',
                *m.jxyr)

            return (item, m)

        ifig = 0
        for vertical in (False, True):
            yield plot_double_sf(source, stats, ifig, vertical)
            ifig += 1
