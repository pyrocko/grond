import os
import os.path as op
import logging

from pyrocko import gf, guts, util
from pyrocko.guts import Dict, List, Tuple, Float, Unicode, Object, String

from grond.plot.config import PlotFormat


logger = logging.getLogger('grond.plot.collection')


class StringID(gf.StringID):
    pass


class PlotItem(Object):
    name = StringID.T()
    attributes = Dict.T(StringID.T(), List.T(String.T()))


class PlotGroup(Object):
    name = StringID.T()
    variant = StringID.T()
    description = Unicode.T(optional=True)
    formats = List.T(PlotFormat.T())
    size_cm = Tuple.T(2, Float.T())
    items = List.T(PlotItem.T())
    attributes = Dict.T(StringID.T(), List.T(String.T()))

    def filename_image(self, item, format):
        return '%s.%s.%s' % (
            self.name,
            self.variant,
            item.name,
            format.extension)


class PlotCollection(Object):
    group_refs = List.T(Tuple.T(2, StringID.T()))


class PlotCollectionManager(object):

    def __init__(self, path):
        self._path = path
        self.load_collection()

    def load_collection(self):
        self._collection = guts.load(filename=self.path_index())

    def dump_collection(self):
        guts.dump(filename=self.path_index())

    def path_collection(self):
        return op.join(self._path, 'plot_collection.yaml')

    def path_image(self, group, item, format):
        return op.join(self._path, group.filename_image(item, format))

    def path_group(self, group_ref=None, group=None):
        if group_ref is not None:
            group_name, group_variant = group_ref
        else:
            group_name = group.name
            group_variant = group.variant

        return op.join(
            self._path, group_name, group_variant)

    def create_group(self, config, iter_item_figure, **kwargs):
        group = PlotGroup(
            formats=guts.clone(config.formats),
            size_cm=config.size_cm,
            name=config.plot_name,
            variant=config.plot_variant,
            **kwargs)

        path_group = self.path_group(group=group)
        if os.path.exists(path_group):
            self.remove_group(group.name, group.variant)

        for item, fig in iter_item_figure:
            group.items.append(item)
            for format in group.formats:
                path = self.path_image(group, item, format)
                util.ensuredirs(path)
                fig.savefig(
                    path,
                    format=format.name,
                    dpi=format.get_dpi(group.size_cm))

                logger.info('figure saved: %s' % path)

        group.dump(filename=path_group)
        self.add(group)

    def remove_group(self, group_name, group_variant):
        group = guts.load(filename=self.path_group(group_name, group_variant))
        for item in group.items:
            for format in group.formats:
                path = self.path_image(group, item, format)
                os.unlink(path)

        path_group = self.path_group(group)
        os.unlink(path_group)


__all__ = [
    'StringID',
    'PlotItem',
    'PlotGroup',
    'PlotCollection',
    'PlotCollectionManager',
]
