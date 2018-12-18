import os
import os.path as op
import logging

from pyrocko import guts, util
from pyrocko.guts import Dict, List, Tuple, Float, Unicode, Object, String

from grond.meta import StringID
from grond.plot.config import PlotFormat


guts_prefix = 'grond'

logger = logging.getLogger('grond.plot.collection')


class PlotItem(Object):
    name = StringID.T()
    attributes = Dict.T(
        StringID.T(), List.T(String.T()))
    title = Unicode.T(
        optional=True,
        help='item\'s description')
    description = Unicode.T(
        optional=True,
        help='item\'s description')


class PlotGroup(Object):
    name = StringID.T(
        help='group name')
    section = StringID.T(
        optional=True,
        help='group\'s section path, e.g. results.waveforms')
    title = Unicode.T(
        optional=True,
        help='group\'s title')
    description = Unicode.T(
        optional=True,
        help='group description')
    formats = List.T(
        PlotFormat.T(),
        help='plot format')
    variant = StringID.T(
        help='variant of the group')
    feather_icon = String.T(
        default='bar-chart-2',
        help='Feather icon for the HTML report.')
    size_cm = Tuple.T(2, Float.T())
    items = List.T(PlotItem.T())
    attributes = Dict.T(StringID.T(), List.T(String.T()))

    def filename_image(self, item, format):
        return '%s.%s.%s.%s' % (
            self.name,
            self.variant,
            item.name,
            format.extension)


class PlotCollection(Object):
    group_refs = List.T(Tuple.T(2, StringID.T()))


class PlotCollectionManager(object):

    def __init__(self, path, show=False):
        self._path = path
        self.load_collection()
        self._show = show

    def load_collection(self):
        path = self.path_collection()
        if op.exists(path):
            self._collection = guts.load(filename=self.path_collection())
        else:
            self._collection = PlotCollection()

    def dump_collection(self):
        path = self.path_collection()
        util.ensuredirs(path)
        guts.dump(self._collection, filename=path)

    def path_collection(self):
        return op.join(self._path, 'plot_collection.yaml')

    def path_image(self, group, item, format):
        return op.join(
            self._path, group.name, group.variant,
            group.filename_image(item, format))

    def path_group(self, group_ref=None, group=None):
        if group_ref is not None:
            group_name, group_variant = group_ref
        else:
            group_name = group.name
            group_variant = group.variant

        return op.join(
            self._path, group_name, group_variant, '%s.%s.%s' % (
                group_name, group_variant, 'plot_group.yaml'))

    def paths_group_dirs(self, group_ref=None, group=None):
        if group_ref is not None:
            group_name, group_variant = group_ref
        else:
            group_name = group.name
            group_variant = group.variant

        return [
            op.join(self._path, group_name, group_variant),
            op.join(self._path, group_name)]

    def create_group_mpl(self, config, iter_item_figure, **kwargs):
        from matplotlib import pyplot as plt
        group = PlotGroup(
            formats=guts.clone(config.formats),
            size_cm=config.size_cm,
            name=config.name,
            variant=config.variant,
            **kwargs)

        path_group = self.path_group(group=group)
        if os.path.exists(path_group):
            self.remove_group_files(path_group)

        group_ref = (group.name, group.variant)
        if group_ref in self._collection.group_refs:
            self._collection.group_refs.remove(group_ref)

        self.dump_collection()

        figs_to_close = []
        for item, fig in iter_item_figure:
            group.items.append(item)
            for format in group.formats:
                path = self.path_image(group, item, format)
                util.ensuredirs(path)
                format.render_mpl(
                    fig,
                    path=path,
                    dpi=format.get_dpi(group.size_cm))

                logger.info('Figure saved: %s' % path)

            if not self._show:
                plt.close(fig)
            else:
                figs_to_close.append(fig)

        util.ensuredirs(path_group)
        group.validate()
        group.dump(filename=path_group)
        self._collection.group_refs.append(group_ref)
        self.dump_collection()

        if self._show:
            plt.show()

        for fig in figs_to_close:
            plt.close(fig)

    def create_group_automap(self, config, iter_item_figure, **kwargs):
        group = PlotGroup(
            formats=guts.clone(config.formats),
            size_cm=config.size_cm,
            name=config.name,
            variant=config.variant,
            **kwargs)

        path_group = self.path_group(group=group)
        if os.path.exists(path_group):
            self.remove_group_files(path_group)

        group_ref = (group.name, group.variant)
        if group_ref in self._collection.group_refs:
            self._collection.group_refs.remove(group_ref)

        self.dump_collection()

        for item, automap in iter_item_figure:
            group.items.append(item)
            for format in group.formats:
                path = self.path_image(group, item, format)
                util.ensuredirs(path)
                format.render_automap(
                    automap,
                    path=path,
                    resolution=format.get_dpi(group.size_cm))

                logger.info('Figure saved: %s' % path)

        util.ensuredirs(path_group)
        group.dump(filename=path_group)
        self._collection.group_refs.append(group_ref)
        self.dump_collection()

    def create_group_gmtpy(self, config, iter_item_figure):
        pass

    def remove_group_files(self, path_group):
        group = guts.load(filename=path_group)
        for item in group.items:
            for format in group.formats:
                path = self.path_image(group, item, format)
                try:
                    os.unlink(path)
                except OSError:
                    pass

        os.unlink(path_group)
        for path in self.paths_group_dirs(group=group):
            try:
                os.rmdir(path)
            except OSError:
                pass


__all__ = [
    'PlotItem',
    'PlotGroup',
    'PlotCollection',
    'PlotCollectionManager',
]
