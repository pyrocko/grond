import numpy as num
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import colors, patheffects

from pyrocko.guts import Tuple, Float
from pyrocko.plot import mpl_init

from grond.plot.config import PlotConfig
from grond.plot.collection import PlotItem

km = 1e3


def km_fmt(x, p):
    return '%.1f' % (x / km)


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
Slip distribution and rake of the finite slip solution. The absolute slip
is color coded. The black star marks the nucleation point.
Contour lines indicate the rupture evolution in %.1f s intervals.
''' % self.dt_contour)

    def draw_figures(self, history, problem):
        source = history.get_best_source()

        store_ids = problem.get_gf_store_ids()
        store = problem.get_gf_store(store_ids[0])

        interpolation = 'nearest_neighbor'

        source.discretize_patches(store, interpolation)
        patches = source.patches
        dislocations = source.get_okada_slip(scale_slip=True)

        patches_x = num.array([p.length_pos for p in patches])\
            .reshape(source.nx, source.ny)
        patches_y = num.array([p.width_pos for p in patches])\
            .reshape(source.nx, source.ny)
        patches_t = num.array([p.time for p in patches])\
            .reshape(source.nx, source.ny)

        fig = plt.figure()
        ax = fig.gca()

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

        for l in labels:
            l.set_rotation(0.)
            l.set_fontweight('semibold')
            l.set_fontsize('small')
            l.set_path_effects([
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
            x_quiver, y_quiver, slip_strike, -slip_dip,
            facecolor='none', edgecolor='k', linewidth=.7,
            scale=10, headwidth=3,
            cmap='YlOrRd', alpha=.6)

        ax.invert_yaxis()

        nucleation_x = ((source.nucleation_x + 1.) / 2.) * source.length
        nucleation_y = ((source.nucleation_y + 1.) / 2.) * source.width

        ax.scatter(
            nucleation_x, nucleation_y,
            s=20, color='k', marker='*',
            alpha=.7)

        ax.xaxis.set_major_formatter(FuncFormatter(km_fmt))
        ax.yaxis.set_major_formatter(FuncFormatter(km_fmt))

        ax.set_xlabel('Length [km]')
        ax.set_ylabel('Width [km]')

        cmap = fig.colorbar(
            im, orientation='horizontal', pad=0.2, shrink=.8,
            format='%.2f m')
        cmap.set_label('Slip')

        yield PlotItem(name='fig_1'), fig
