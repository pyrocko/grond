import numpy as num
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import colors, patheffects

from pyrocko.guts import Tuple, Float, Bool
from pyrocko.plot import mpl_init, mpl_graph_color, mpl_color
from pyrocko.plot.dynamic_rupture import RuptureMap

from grond.plot.config import PlotConfig
from grond.plot.collection import PlotItem

km = 1e3


def km_fmt(x, p):
    return '%.1f' % (x / km)


def get_source(history, source_type='mean'):
    if source_type == 'mean':
        return history.get_mean_source()
    elif source_type == 'best':
        return history.get_best_source()
    else:
        raise ValueError('current source_type is not defined.')


class DynamicRuptureSlipMap(PlotConfig):
    '''
    Slip map of best solution.
    '''
    name = 'rupture_slip_map'
    dt_contour = Float.T(
        default=0.5,
        help='Rupture propagation contourline interval in seconds.')
    size_cm = Tuple.T(2, Float.T(), default=(20., 20.))

    def make(self, environ):
        environ.setup_modelling()
        cm = environ.get_plot_collection_manager()
        history = environ.get_history(subset='harvest')
        problem = environ.get_problem()

        mpl_init(fontsize=self.font_size)
        cm.create_group_mpl(
            self,
            self.draw_figures(history, problem),
            title=u'Slip Map',
            section='solution',
            feather_icon='grid',  # alternatively: wind
            description=u'''
Slip distribution and rake of the pseudo dynamic rupture model. The absolute
slip is color coded. The arrows are scaled based on the in plane slip. The
black star marks the nucleation point. Contour lines indicate the rupture
evolution in %.1f s intervals.
''' % self.dt_contour)

    def draw_figures(self, history, problem):
        store_ids = problem.get_gf_store_ids()
        store = problem.get_gf_store(store_ids[0])

        interpolation = 'nearest_neighbor'

        for i, plabel in enumerate(('best', 'mean')):
            source = get_source(history, source_type=plabel)

            fig, ax = plt.subplots(1, 1)

            # ToDo in function with "mean", "best" as arg
            source.discretize_patches(store, interpolation)
            patches = source.patches
            dislocations = source.get_okada_slip(scale_slip=True)

            patches_x = num.array([p.ix for p in patches])\
                .reshape(source.nx, source.ny)
            patches_y = num.array([p.iy for p in patches])\
                .reshape(source.nx, source.ny)
            patches_t = num.array([p.time for p in patches])\
                .reshape(source.nx, source.ny)

            abs_disloc = num.linalg.norm(dislocations, axis=1)
            abs_disloc = abs_disloc.reshape(source.nx, source.ny)

            im = ax.imshow(
                abs_disloc.T,
                cmap='YlOrRd',
                origin='upper',
                aspect='equal',
                extent=(0., source.length, 0., source.width))

            patches_t -= patches_t.min()

            nlevels = patches_t.max() // self.dt_contour
            contours = num.arange(nlevels) * self.dt_contour
            contours += patches_t.min()

            def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
                return colors.LinearSegmentedColormap.from_list(
                    'trunc({n},{a:.2f},{b:.2f})'.format(
                        n=cmap.name, a=minval, b=maxval),
                    cmap(num.linspace(minval, maxval, n)))

            cmap = truncate_colormap(plt.get_cmap('winter'), 0., 0.8)

            contour = ax.contour(
                patches_x + source.length/2, patches_y, patches_t,
                levels=contours, alpha=.8, colors='k')

            labels = ax.clabel(
                contour, contour.levels[::2],
                inline=True, fmt='%.1f s')

            for label in labels:
                label.set_rotation(0.)
                label.set_fontweight('semibold')
                label.set_fontsize('small')
                label.set_path_effects([
                    patheffects.Stroke(linewidth=1.25, foreground='beige'),
                    patheffects.Normal()])

            slip_strike = dislocations[:, 0]
            slip_dip = dislocations[:, 1]

            patch_length = source.length / source.nx
            patch_width = source.width / source.ny

            x_quiver = num.repeat(num.arange(source.nx), source.ny)\
                * patch_length + patch_length/2.
            y_quiver = num.tile(num.arange(source.ny), source.nx)\
                * patch_width + patch_width/2.

            ax.quiver(
                x_quiver, y_quiver, slip_strike, slip_dip,
                facecolor='none', edgecolor='k', linewidth=.7,
                scale=75 / abs_disloc.max(), headwidth=3,
                cmap='YlOrRd', alpha=.6)

            ax.invert_yaxis()

            nucleation_x = ((source.nucleation_x + 1.) / 2.) * source.length
            nucleation_y = ((source.nucleation_y + 1.) / 2.) * source.width

            ax.scatter(
                nucleation_x, nucleation_y,
                s=60, color='w', edgecolor='k', marker='*',
                linewidths=.5, alpha=.7)

            ax.xaxis.set_major_formatter(FuncFormatter(km_fmt))
            ax.yaxis.set_major_formatter(FuncFormatter(km_fmt))

            ax.set_xlabel('Length [km]')
            ax.set_ylabel('Width [km]')

            cmap = fig.colorbar(
                im, orientation='horizontal', pad=0.2, shrink=.8,
                ax=ax, format='%.2f')
            cmap.set_label('Slip [m]')

            item = PlotItem(
                name='_'.join(['ensemble', plabel.lower()]),
                title='ensemble %s slip distribution' % (plabel.lower()),
                description=u'')

            yield item, fig


class DynamicRuptureSTF(PlotConfig):
    '''
    Slip map of best solution.
    '''
    name = 'rupture_source_time_function'
    dt_sampling = Float.T(
        default=2.5,
        help='Moment rate sampling interval in seconds')
    size_cm = Tuple.T(2, Float.T(), default=(15., 15.))

    def make(self, environ):
        environ.setup_modelling()
        cm = environ.get_plot_collection_manager()
        history = environ.get_history(subset='harvest')
        problem = environ.get_problem()

        mpl_init(fontsize=self.font_size)
        cm.create_group_mpl(
            self,
            self.draw_figures(history, problem),
            title=u'Moment Rate Function (STF)',
            section='solution',
            feather_icon='zap',  # alternatively: activity
            description=u'''
Source time function (moment release) of the pseudo dynamic rupture model.
The moment rate function is sampled in %.1f s intervals.
''' % self.dt_sampling)

    def draw_figures(self, history, problem):
        store_ids = problem.get_gf_store_ids()
        store = problem.get_gf_store(store_ids[0])

        interpolation = 'nearest_neighbor'

        for i, plabel in enumerate(('best', 'mean')):
            source = get_source(history, source_type=plabel)

            fig, ax = plt.subplots(1, 1)

            source.discretize_patches(store, interpolation)

            mrate, times = source.get_moment_rate(
                store=store, deltat=self.dt_sampling)

            mrate_max = mrate.max()

            ax.fill_between(
                times,
                mrate / mrate_max,
                color=mpl_color(x='aluminium2'),
                alpha=0.7)

            ax.scatter(
                times,
                mrate / mrate_max,
                c=mpl_graph_color(0))

            ax.set_xlabel('Time [s]')
            ax.set_ylabel(r'$\dot{M}$ / %.2e Nm/s' % mrate_max)

            ax.set_xlim([0, times.max() + self.dt_sampling])
            ax.grid(True)

            item = PlotItem(
                name='_'.join(['ensemble', plabel.lower()]),
                title='ensemble %s source time function' % (plabel.lower()),
                description=u'')

            yield item, fig


class DynamicRuptureMap(PlotConfig):
    '''
    Overview map of best solution.
    '''
    name = 'rupture_overview_map'
    size_cm = Tuple.T(
        2, Float.T(),
        default=(20., 20.),
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
        environ.setup_modelling()
        cm = environ.get_plot_collection_manager()
        history = environ.get_history(subset='harvest')
        problem = environ.get_problem()

        cm.create_group_automap(
            self,
            self.draw_rupture_map(history, problem),
            title=u'Rupture Dislocations',
            section='solution',
            feather_icon='map',
            description=u'''
Maps showing orientation and dislocation of the PseudoDynamicRupture.
''')

    def draw_rupture_map(self, history, problem):
        from pyrocko import orthodrome as pod

        store_ids = problem.get_gf_store_ids()
        store = problem.get_gf_store(store_ids[0])

        def plot_rupture(source, model_name):
            item = PlotItem(
                name='ensemble_%s' % model_name,
                attributes={},
                title=u'Final static rupture dislocation - %s model'
                      % model_name,
                description=u'''
Static rupture dislocation from %s model.''' % model_name)

            lat, lon = pod.ne_to_latlon(
                source.lat,
                source.lon,
                source.north_shift,
                source.east_shift)

            map_kwargs = dict(
                lat=source.lat,
                lon=source.lon,
                radius=self.radius or source.length,
                width=self.size_cm[0],
                height=self.size_cm[1],
                source=source,
                show_topo=self.show_topo,
                show_grid=self.show_grid,
                show_rivers=self.show_rivers,
                color_wet=(216, 242, 254),
                color_dry=(238, 236, 230))

            source.discretize_patches(store, interpolation='nearest_neighbor')

            m = RuptureMap(**map_kwargs)
            m.draw_dislocation(cmap='summer')
            m.draw_time_contour(store)
            m.draw_nucleation_point()
            m.draw_top_edge()

            return (item, m)

        for i, plabel in enumerate(('best', 'mean')):
            source = get_source(history, source_type=plabel)

            yield plot_rupture(source, plabel)
