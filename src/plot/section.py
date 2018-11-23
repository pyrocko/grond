from __future__ import print_function
import numpy as num
from matplotlib.axes import Axes
from matplotlib.ticker import MultipleLocator

from pyrocko.guts import Tuple, Float
from pyrocko import plot

from .config import PlotConfig

guts_prefix = 'grond'


def limits(points):
    lims = num.zeros((3, 2))
    if points.size != 0:
        lims[:, 0] = num.min(points, axis=0)
        lims[:, 1] = num.max(points, axis=0)

    return lims


class NotEnoughSpace(Exception):
    pass


class SectionPlotConfig(PlotConfig):

    size_cm = Tuple.T(
        2, Float.T(), default=(20., 20.))

    margins_em = Tuple.T(
        4, Float.T(), default=(7., 5., 7., 5.))

    separator_em = Float.T(default=1.0)


class SectionPlot(object):

    def __init__(self, config=None):
        if config is None:
            config = SectionPlotConfig()

        self.config = config
        self._disconnect_data = []
        self._width = self._height = self._pixels = None
        self._plt = plot.mpl_init(self.config.font_size)
        self._fig = fig = self._plt.figure(figsize=self.config.size_inch)

        rect = [0., 0., 1., 1.]
        self._axes_xy = Axes(fig, rect)
        self._axes_xz = Axes(fig, rect)
        self._axes_zy = Axes(fig, rect)

        self._view_limits = num.zeros((3, 2))

        self._view_limits[:, :] = num.nan

        self._update_geometry()

        for axes in self.axes_list:
            fig.add_axes(axes)
            self._connect(axes, 'xlim_changed', self.lim_changed_handler)
            self._connect(axes, 'ylim_changed', self.lim_changed_handler)

        self._cid_resize = fig.canvas.mpl_connect(
            'resize_event', self.resize_handler)

        self._connect(fig, 'dpi_changed', self.dpi_changed_handler)

        self._lim_changed_depth = 0

    def _connect(self, obj, sig, handler):
        cid = obj.callbacks.connect(sig, handler)
        self._disconnect_data.append((obj, cid))

    def _disconnect_all(self):
        for obj, cid in self._disconnect_data:
            obj.callbacks.disconnect(cid)

        self._fig.canvas.mpl_disconnect(self._cid_resize)

    def dpi_changed_handler(self, fig):
        self._update_geometry()

    def resize_handler(self, event):
        self._update_geometry()

    def lim_changed_handler(self, axes):
        self._lim_changed_depth += 1
        if self._lim_changed_depth < 2:
            self._update_layout()

        self._lim_changed_depth -= 1

    def _update_geometry(self):
        w, h = self._fig.canvas.get_width_height()
        p = self.get_pixels_factor()

        if (self._width, self._height, self._pixels) != (w, h, p):
            self._width = w
            self._height = h
            self._pixels = p
            self._update_layout()

    @property
    def margins(self):
        return tuple(
            x * self.config.font_size / self._pixels
            for x in self.config.margins_em)

    @property
    def separator(self):
        return self.config.separator_em * self.config.font_size / self._pixels

    def rect_to_figure_coords(self, rect):
        left, bottom, width, height = rect
        return (
            left / self._width,
            bottom / self._height,
            width / self._width,
            height / self._height)

    def point_to_axes_coords(self, axes, point):
        x, y = point
        aleft, abottom, awidth, aheight = axes.get_position().bounds

        x_fig = x / self._width
        y_fig = y / self._height

        x_axes = (x_fig - aleft) / awidth
        y_axes = (y_fig - abottom) / aheight

        return (x_axes, y_axes)

    def get_pixels_factor(self):
        try:
            r = self._fig.canvas.get_renderer()
            return 1.0 / r.points_to_pixels(1.0)
        except AttributeError:
            return 1.0

    def make_limits(self, lims):
        a = plot.AutoScaler(space=0.05)
        return a.make_scale(lims)[:2]

    def get_data_limits(self):
        xs = []
        ys = []
        zs = []
        xs.extend(self._axes_xy.get_xaxis().get_data_interval())
        ys.extend(self._axes_xy.get_yaxis().get_data_interval())
        xs.extend(self._axes_xz.get_xaxis().get_data_interval())
        zs.extend(self._axes_xz.get_yaxis().get_data_interval())
        zs.extend(self._axes_zy.get_xaxis().get_data_interval())
        ys.extend(self._axes_zy.get_yaxis().get_data_interval())
        lims = num.zeros((3, 2))
        lims[0, :] = num.nanmin(xs), num.nanmax(xs)
        lims[1, :] = num.nanmin(ys), num.nanmax(ys)
        lims[2, :] = num.nanmin(zs), num.nanmax(zs)
        lims[num.logical_not(num.isfinite(lims))] = 0.0
        return lims

    def set_xlim(self, xmin, xmax):
        self._view_limits[0, :] = xmin, xmax
        self._update_layout()

    def set_ylim(self, ymin, ymax):
        self._view_limits[1, :] = ymin, ymax
        self._update_layout()

    def set_zlim(self, zmin, zmax):
        self._view_limits[2, :] = zmin, zmax
        self._update_layout()

    def _update_layout(self):
        data_limits = self.get_data_limits()

        limits = num.zeros((3, 2))
        for i in range(3):
            limits[i, :] = self.make_limits(data_limits[i, :])

        mask = num.isfinite(self._view_limits)
        limits[mask] = self._view_limits[mask]

        deltas = limits[:, 1] - limits[:, 0]

        data_w = deltas[0] + deltas[2]
        data_h = deltas[1] + deltas[2]

        ml, mt, mr, mb = self.margins
        ms = self.separator

        data_r = data_h / data_w
        em = self.config.font_size
        w = self._width
        h = self._height
        fig_w_avail = w - mr - ml - ms
        fig_h_avail = h - mt - mb - ms

        if fig_w_avail <= 0.0 or fig_h_avail <= 0.0:
            raise NotEnoughSpace()

        fig_r = fig_h_avail / fig_w_avail

        if data_r < fig_r:
            data_expanded_h = data_w * fig_r
            data_expanded_w = data_w
        else:
            data_expanded_h = data_h
            data_expanded_w = data_h / fig_r

        limits[0, 0] -= 0.5 * (data_expanded_w - data_w)
        limits[0, 1] += 0.5 * (data_expanded_w - data_w)
        limits[1, 0] -= 0.5 * (data_expanded_h - data_h)
        limits[1, 1] += 0.5 * (data_expanded_h - data_h)

        deltas = limits[:, 1] - limits[:, 0]

        w1 = fig_w_avail * deltas[0] / data_expanded_w
        w2 = fig_w_avail * deltas[2] / data_expanded_w

        h1 = fig_h_avail * deltas[1] / data_expanded_h
        h2 = fig_h_avail * deltas[2] / data_expanded_h

        rect_xy = [ml, mb+h2+ms, w1, h1]
        rect_xz = [ml, mb, w1, h2]
        rect_zy = [ml+w1+ms, mb+h2+ms, w2, h1]

        axes_xy, axes_xz, axes_zy = self.axes_list

        axes_xy.set_position(
            self.rect_to_figure_coords(rect_xy), which='both')
        axes_xz.set_position(
            self.rect_to_figure_coords(rect_xz), which='both')
        axes_zy.set_position(
            self.rect_to_figure_coords(rect_zy), which='both')

        def wcenter(rect):
            return rect[0] + rect[2]*0.5

        def hcenter(rect):
            return rect[1] + rect[3]*0.5

        self.set_label_coords(
            axes_xy, 'x', [wcenter(rect_xy), h - 1.0*em])
        self.set_label_coords(
            axes_xy, 'y', [2.0*em, hcenter(rect_xy)])
        self.set_label_coords(
            axes_zy, 'x', [wcenter(rect_zy), h - 1.0*em])
        self.set_label_coords(
            axes_xz, 'y', [2.0*em, hcenter(rect_xz)])

        scaler = plot.AutoScaler()
        inc = scaler.make_scale(
            [0, min(data_expanded_w, data_expanded_h)], override_mode='off')[2]

        axes_xy.set_xlim(*limits[0, :])
        axes_xy.set_ylim(*limits[1, :])
        axes_xy.get_xaxis().set_tick_params(
            bottom=False, top=True, labelbottom=False, labeltop=True)
        axes_xy.get_yaxis().set_tick_params(
            left=True, labelleft=True, right=False, labelright=False)

        axes_xz.set_xlim(*limits[0, :])
        axes_xz.set_ylim(*limits[2, ::-1])
        axes_xz.get_xaxis().set_tick_params(
            bottom=True, top=False, labelbottom=False, labeltop=False)
        axes_xz.get_yaxis().set_tick_params(
            left=True, labelleft=True, right=True, labelright=False)

        axes_zy.set_xlim(*limits[2, :])
        axes_zy.set_ylim(*limits[1, :])
        axes_zy.get_xaxis().set_tick_params(
            bottom=True, top=True, labelbottom=False, labeltop=True)
        axes_zy.get_yaxis().set_tick_params(
            left=False, labelleft=False, right=True, labelright=False)

        for axes in self.axes_list:
            tl = MultipleLocator(inc)
            axes.get_xaxis().set_major_locator(tl)
            tl = MultipleLocator(inc)
            axes.get_yaxis().set_major_locator(tl)

    def set_label_coords(self, axes, which, point):
        axis = axes.get_xaxis() if which == 'x' else axes.get_yaxis()
        axis.set_label_coords(*self.point_to_axes_coords(axes, point))

    @property
    def fig(self):
        return self._fig

    @property
    def axes_xy(self):
        return self._axes_xy

    @property
    def axes_xz(self):
        return self._axes_xz

    @property
    def axes_zy(self):
        return self._axes_zy

    @property
    def axes_list(self):
        return [
            self._axes_xy, self._axes_xz, self._axes_zy]

    def plot(self, points, *args, **kwargs):
        self._axes_xy.plot(points[:, 0], points[:, 1], *args, **kwargs)
        self._axes_xz.plot(points[:, 0], points[:, 2], *args, **kwargs)
        self._axes_zy.plot(points[:, 2], points[:, 1], *args, **kwargs)

    def close(self):
        self._disconnect_all()
        self._plt.close(self._fig)

    def show(self):
        self._plt.show()

    def set_xlabel(self, s):
        self._axes_xy.set_xlabel(s)

    def set_ylabel(self, s):
        self._axes_xy.set_ylabel(s)

    def set_zlabel(self, s):
        self._axes_xz.set_ylabel(s)
        self._axes_zy.set_xlabel(s)
