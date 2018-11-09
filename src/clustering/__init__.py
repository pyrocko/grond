import optparse

from pyrocko.guts import Object, Int, Float
from pyrocko import guts

from grond.version import __version__
from grond.meta import GrondError

guts_prefix = 'grond'


class Clustering(Object):
    '''Base class for clustering method configuration objects.'''

    def perform(self):
        raise NotImplementedError('should be implemented in subclass')

    @classmethod
    def _cli_setup(cls, parser):

        pmap = {
            Float.T: float,
            Int.T: int}

        group = optparse.OptionGroup(
            parser, cls.name,
            'Options specific for the "%s" clustering method' % cls.name)

        for prop in cls.T.properties:
            if isinstance(prop, tuple(pmap.keys())):
                group.add_option(
                    '--%s' % u2d(prop.name),
                    dest=prop.name,
                    type=pmap[prop.__class__],
                    default=prop.default(),
                    help=prop.help + ' (default: %default)')

        parser.add_option_group(group)

    @staticmethod
    def cli_setup(name, setup):
        if name in Clustering.name_to_class:
            def setup_(parser):
                setup(parser)
                Clustering.name_to_class[name]._cli_setup(parser)

            return setup_

        else:
            return setup

    @classmethod
    def _cli_instantiate(cls, options):
        pmap = {
            Float.T: float,
            Int.T: int}

        kwargs = {}
        for prop in cls.T.properties:
            if isinstance(prop, tuple(pmap.keys())):
                kwargs[prop.name] = getattr(options, prop.name)

        return cls(**kwargs)

    @staticmethod
    def cli_instantiate(name, options):
        return Clustering.name_to_class[name]._cli_instantiate(options)


def u2d(u):
    return u.replace('_', '-')


class DBScan(Clustering):
    '''DBSCAN clustering algorithm.'''

    name = 'dbscan'

    nmin = Int.T(
        default=10,
        help='Minimum number of neighbours to define a cluster.')

    eps = Float.T(
        default=0.1,
        help='Maximum distance to search for neighbors.')

    ncluster_limit = Int.T(
        default=None,
        help='Limit maximum number of clusters created to N.')

    def perform(self, similarity_matrix):
        from .dbscan import dbscan
        return dbscan(
            similarity_matrix,
            nmin=self.nmin,
            eps=self.eps,
            ncluster_limit=self.ncluster_limit)


def read_config(path):
    try:
        config = guts.load(filename=path)
    except OSError:
        raise GrondError(
            'cannot read Grond clustering configuration file: %s' % path)

    if not isinstance(config, Clustering):
        raise GrondError(
            'invalid Grond clustering configuration in file "%s"' % path)

    return config


def write_config(config, path):
    try:
        guts.dump(
            config,
            filename=path,
            header='Grond clustering configuration file, version %s'
                   % __version__)

    except OSError:
        raise GrondError(
            'cannot write Grond report configuration file: %s' % path)


Clustering.name_to_class = dict((cls.name, cls) for cls in [DBScan])

methods = sorted(Clustering.name_to_class.keys())

__all__ = [
    'Clustering',
    'DBScan',
]
