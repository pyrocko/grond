from pyrocko.guts import Object, Float, Int, List, Tuple, String, load, clone

from grond.meta import GrondError

guts_prefix = 'grond'


inch = 2.54
points = inch / 72.0


class PlotFormat(Object):

    @property
    def extension(self):
        return self.name

    def get_dpi(self, size_cm):
        return None

    def render_mpl(self, fig, path, **kwargs):
        raise NotImplementedError

    def render_automap(self, automap, path, **kwargs):
        raise NotImplementedError


class PNG(PlotFormat):
    name = 'png'

    dpi = Int.T(
        optional=True,
        help='DPI of the figure')
    size_pixels = Int.T(
        optional=True,
        help='Size in pixels')
    width_pixels = Int.T(
        optional=True,
        help='Width in pixels')
    height_pixels = Int.T(
        optional=True,
        help='Height in pixels')

    @property
    def extension(self):
        if self.dpi is not None:
            return 'd%i.png' % self.dpi
        elif self.size_pixels is not None:
            return 's%i.png' % self.size_pixels
        elif self.width_pixels is not None:
            return 'w%i.png' % self.width_pixels
        elif self.height_pixels is not None:
            return 'h%i.png' % self.height_pixels
        else:
            return 'd100.png'

    def get_dpi(self, size_cm):
        w_cm, h_cm = size_cm
        w_inch, h_inch = w_cm/inch, h_cm/inch
        if self.dpi:
            return self.dpi
        elif self.size_pixels is not None:
            return min(self.size_pixels/w_inch, self.size_pixels/h_inch)
        elif self.width_pixels is not None:
            return self.width_pixels/w_inch
        elif self.height_pixels is not None:
            return self.height_pixels/h_inch
        else:
            return 100.0

    def render_mpl(self, fig, path, **kwargs):
        return fig.savefig(path, format=self.name, **kwargs)

    def render_automap(self, automap, path, **kwargs):
        return automap.save(path, **kwargs)


class PDF(PlotFormat):
    name = 'pdf'

    dpi = Int.T(
        default=150,
        help='DPI of the figure')

    def get_dpi(self, size_cm):
        return self.dpi

    def render_mpl(self, fig, path, **kwargs):
        return fig.savefig(path, format=self.name, **kwargs)

    def render_automap(self, automap, path, **kwargs):
        return automap.save(path, **kwargs)


class SVG(PlotFormat):
    name = 'svg'

    dpi = Int.T(
        default=150,
        help='DPI of the figure')

    def get_dpi(self, size_cm):
        return self.dpi

    def render_mpl(self, fig, path, **kwargs):
        return fig.savefig(path, format=self.name, **kwargs)

    def render_automap(self, automap, path, **kwargs):
        return automap.save(path, **kwargs)


class HTML(PlotFormat):
    name = 'html'

    @property
    def extension(self):
        return 'html'

    def render_mpl(self, fig, path, **kwargs):
        import mpld3
        kwargs.pop('dpi')

        mpld3.save_html(
            fig,
            fileobj=path,
            **kwargs)


class PlotConfig(Object):
    name = 'undefined'
    variant = String.T(
        default='default',
        help='Variant of the plot (if applicable)')
    formats = List.T(
        PlotFormat.T(),
        default=[PNG()],
        help='Format of the plot')
    size_cm = Tuple.T(
        2, Float.T(),
        help='size of the plot')
    font_size = Float.T(
        default=10.,
        help='font size')

    @property
    def size_inch(self):
        return self.size_cm[0]/inch, self.size_cm[1]/inch

    @property
    def size_points(self):
        return self.size_cm[0]/points, self.size_cm[1]/points

    def make(self, environ):
        pass


class PlotConfigCollection(Object):
    plot_configs = List.T(PlotConfig.T())

    @classmethod
    def load(cls, path):
        from grond.plot import get_all_plot_classes
        get_all_plot_classes()  # make sure all plot classes are loaded
        collection = load(filename=path)
        if not isinstance(collection, PlotConfigCollection):
            raise GrondError(
                'invalid plot collection configuration in file "%s"' % path)

        return collection

    def get_weeded(self, env):
        '''Get subset of plot configs supported by the current environment.'''

        plot_classes = env.get_plot_classes()
        plot_configs_avail = [
            clone(pc)
            for pc in self.plot_configs
            if pc.__class__ in plot_classes]

        return PlotConfigCollection(plot_configs=plot_configs_avail)


__all__ = [
    'PlotFormat',
    'PNG',
    'PDF',
    'PlotConfig',
    'PlotConfigCollection',
]
