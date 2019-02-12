import logging
from pyrocko import guts
from pyrocko.guts import List, Unicode, Object

from grond.version import __version__
from grond.meta import GrondError

guts_prefix = 'grond'
logger = logging.getLogger('grond.report')


class RunInfo(Object):
    tags = List.T(
        Unicode.T(),
        help='List of user defined labels')

    def add_tag(self, tag):
        if tag not in self.tags:
            self.tags.append(tag)
            self.tags.sort()
        else:
            logger.warn('While adding tag: tag already set: %s' % tag)

    def remove_tag(self, tag):
        try:
            self.tags.remove(tag)
        except ValueError:
            logger.warn('While removing tag: tag not set: %s' % tag)


def read_info(path):
    try:
        info = guts.load(filename=path)
    except OSError:
        raise GrondError(
            'Cannot read Grond run info file: %s' % path)

    if not isinstance(info, RunInfo):
        raise GrondError(
            'Invalid Grond run info in file "%s".' % path)

    return info


def write_info(info, path):
    try:
        guts.dump(
            info,
            filename=path,
            header='Grond run info file, version %s' % __version__)

    except OSError:
        raise GrondError(
            'Cannot write Grond run info file: %s' % path)
