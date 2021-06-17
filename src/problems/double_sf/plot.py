import math
import logging

import numpy as num

from matplotlib import pyplot as plt

from pyrocko.guts import Float, Bool, Tuple
from pyrocko import gf
from pyrocko.plot import automap, mpl_init, beachball, mpl_color

from grond import stats
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
shows the 95%% confidence ellipse. Orange points indicate the location of each
individual single force, while the square shows the combined centroid location.
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
                lat=event.effective_lat,
                lon=event.effective_lon,
                radius=radius,
                show_topo=self.show_topo,
                show_grid=self.show_grid,
                show_rivers=self.show_rivers,
                color_wet=(216, 242, 254),
                color_dry=(238, 236, 230))

            source.lat, source.lon = event.effective_lat, event.effective_lon

            if isinstance(source, gf.DoubleSFSource):
                sf1, sf2 = source.split()
                offset_scale = source.force
                f = 'r'
            elif isinstance(source, gf.CombiSFSource):
                sf1, sf2 = source.subsources

                sf1.lat, sf1.lon = source.lat, source.lon
                sf2.lat, sf2.lon = source.lat, source.lon

                offset_scale = num.max((sf1.force, sf2.force))
                f = ''

            size = num.linalg.norm(self.size_cm)

            scale = (size / 5.) / offset_scale

            stats_dict = stats.get_values_dict()

            if vertical:
                rows = [[
                    sf1.effective_lon, sf1.effective_lat,
                    0., -sf1.fd * scale,
                    (stats_dict[f + 'fn1.std'] + stats_dict[f + 'fe1.std']),
                    stats_dict[f + 'fd1.std'],
                    0.]]

                rows.append([
                    sf2.effective_lon, sf2.effective_lat,
                    0., -sf2.fd * scale,
                    (stats_dict[f + 'fn2.std'] + stats_dict[f + 'fe2.std']),
                    stats_dict[f + 'fd2.std'],
                    0.])

            else:
                rows = [[
                    sf1.effective_lon, sf1.effective_lat,
                    sf1.fe * scale, sf1.fn * scale,
                    stats_dict[f + 'fe1.std'],
                    stats_dict[f + 'fn1.std'],
                    0.]]

                rows.append([
                    sf2.effective_lon, sf2.effective_lat,
                    sf2.fe * scale, sf2.fn * scale,
                    stats_dict[f + 'fe2.std'],
                    stats_dict[f + 'fn2.std'],
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
                S='c8p',
                in_rows=[
                    [sf.effective_lon, sf.effective_lat] for sf in (sf1, sf2)],
                W='1p,black',
                G='orange3',
                *m.jxyr)

            m.gmt.psxy(
                S='s8p',
                in_rows=[[source.effective_lon, source.effective_lat]],
                W='1p,black',
                G='slategray3',
                *m.jxyr)

            return (item, m)

        ifig = 0
        for vertical in (False, True):
            yield plot_double_sf(source, stats, ifig, vertical)
            ifig += 1


class DoubleSFDecompositionPlot(PlotConfig):
    '''
    Double Single Force decomposition plot.
    '''

    name = 'sf_decomposition'
    size_cm = Tuple.T(2, Float.T(), default=(15., 5.))

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        history = environ.get_history(subset='harvest')
        mpl_init(fontsize=self.font_size)
        cm.create_group_mpl(
            self,
            self.draw_figures(history),
            title=u'Single Force Decomposition',
            section='solution',
            feather_icon='sun',
            description=u'''
Double Single Force decomposition of the best-fitting solution into its single
force components.

Shown are the ensemble best and the ensemble mean. The symbol size indicates
the relative strength of the components. The inversion result is consistent
and stable if ensemble mean and ensemble
best have similar symbol size and patterns.
''')

    def draw_figures(self, history):

        fontsize = self.font_size

        fig = plt.figure(figsize=self.size_inch)
        axes = fig.add_subplot(1, 1, 1, aspect=1.0)
        fig.subplots_adjust(left=0., right=1., bottom=0., top=1.)

        problem = history.problem
        models = history.models

        if models.size == 0:
            logger.warn('Empty models vector.')
            return []

        # ref_source = problem.base_source

        mean_source = stats.get_mean_source(
            problem, history.models)

        best_source = history.get_best_source()

        nlines_max = int(round(self.size_cm[1] / 5. * 4. - 1.0))

        def get_deco(source):
            if isinstance(source, gf.DoubleSFSource):
                return [source] + source.split()
            elif isinstance(source, gf.CombiSFSource):
                sf1, sf2 = source.subsources
                f1 = num.array([sf1.fn, sf1.fe, sf1.fd])
                f2 = num.array([sf2.fn, sf2.fe, sf2.fd])

                force = num.linalg.norm([f1 + f2])

                mix = sf2.force / force

                source = gf.DoubleSFSource(
                    force=force,
                    rfn1=sf1.fn / force,
                    rfe1=sf1.fe / force,
                    rfd1=sf1.fd / force,
                    rfn2=sf2.fn / force,
                    rfe2=sf2.fe / force,
                    rfd2=sf2.fd / force,
                    mix=mix)

                return [source, sf1, sf2]

        lines = []
        lines.append(
            ('Ensemble best', get_deco(best_source), mpl_color('aluminium5')))

        lines.append(
            ('Ensemble mean', get_deco(mean_source), mpl_color('aluminium5')))

        force_max = max(sf.force for (_, line, _) in lines for sf in line)

        for xpos, label in [
                (0., 'Double SF'),
                (2., 'SF 1'),
                (4., 'SF 2')]:

            axes.annotate(
                label,
                xy=(1 + xpos, nlines_max),
                xycoords='data',
                xytext=(0., 0.),
                textcoords='offset points',
                ha='center',
                va='center',
                color='black',
                fontsize=fontsize)

        for i, (label, deco, color_t) in enumerate(lines):
            ypos = nlines_max - i - 1.0

            [dsf, sf1, sf2] = deco

            size0 = dsf.force / force_max

            axes.annotate(
                label,
                xy=(-2., ypos),
                xycoords='data',
                xytext=(0., 0.),
                textcoords='offset points',
                ha='left',
                va='center',
                color='black',
                fontsize=fontsize)

            for xpos, sf_part, ratio, ops in [
                    (0., dsf, 1., '='),
                    (2., sf1, sf1.force / force_max, '+'),
                    (4., sf2, sf2.force / force_max, None)]:

                if ratio > 1e-4:
                    try:
                        beachball.plot_singleforce_beachball_mpl(
                            sf_part.fn, sf_part.fe, sf_part.fd, axes,
                            position=(1. + xpos, ypos),
                            size=0.9 * size0 * math.sqrt(ratio),
                            size_units='data',
                            color_t=color_t,
                            linewidth=1.0)

                    except beachball.BeachballError as e:
                        logger.warn(str(e))

                        axes.annotate(
                            'ERROR',
                            xy=(1. + xpos, ypos),
                            ha='center',
                            va='center',
                            color='red',
                            fontsize=fontsize)

                else:
                    axes.annotate(
                        'N/A',
                        xy=(1. + xpos, ypos),
                        ha='center',
                        va='center',
                        color='black',
                        fontsize=fontsize)

                if ops is not None:
                    axes.annotate(
                        ops,
                        xy=(2. + xpos, ypos),
                        ha='center',
                        va='center',
                        color='black',
                        fontsize=fontsize)

        axes.axison = False
        axes.set_xlim(-2.25, 9.75)
        axes.set_ylim(-0.5, nlines_max+0.5)

        item = PlotItem(name='main')
        return [[item, fig]]
