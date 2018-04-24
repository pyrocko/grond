from pyrocko.guts import Object, Float, Int, List, Tuple, String, load

from grond.meta import GrondError

guts_prefix = 'grond'


inch = 2.54


class PlotFormat(Object):

    @property
    def extension(self):
        return self.name

    def get_dpi(self, size_cm):
        return None

    def render_mpl(self, fig, path, **kwargs):
        pass


class PNG(PlotFormat):
    name = 'png'

    dpi = Float.T(
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
        fig.savefig(
            path,
            format=self.name,
            **kwargs)


class PDF(PlotFormat):
    name = 'pdf'

    def render_mpl(self, fig, path, **kwargs):
        fig.savefig(
            path,
            format=self.name,
            **kwargs)


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

    def make(self, environ):
        pass


class PlotConfigCollection(Object):
    plot_configs = List.T(PlotConfig.T())

    @classmethod
    def load(cls, path):
        collection = load(filename=path)
        if not isinstance(collection, PlotConfigCollection):
            raise GrondError(
                'invalid plot collection configuration in file "%s"' % path)

        return collection


__all__ = [
    'PlotFormat',
    'PNG',
    'PDF',
    'PlotConfig',
]
