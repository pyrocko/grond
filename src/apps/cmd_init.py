import logging

import glob
import os.path as op

logger = logging.getLogger('grond.init')
km = 1e3


class GrondInit(object):

    folder = op.join(op.dirname(__file__), '..', 'data', 'init')

    def __init__(self):
        pass

    def available_inits(self):
        return {
            'examples': tuple((self.filename_to_abbrv(fn),
                               self._get_description(fn))
                              for fn in self.example_files),
            'events': tuple((self.filename_to_abbrv(fn),
                             self._get_description(fn))
                            for fn in self.event_files),
        }

    @staticmethod
    def filename_to_abbrv(filename):
        return op.basename(filename).split('.')[0]

    @property
    def example_files(self):
        return self._get_files('example_*.gronf')

    @property
    def section_files(self):
        return self._get_files('section_*.gronf')

    @property
    def snippet_files(self):
        return self._get_files('snippet_*.gronf')

    @property
    def event_files(self):
        return self._get_files('event_*.txt')

    def _get_files(self, name):
        return glob.glob(op.join(self.folder, name))

    def _get_description(self, filename):
        with open(filename, 'rt') as f:
            for ln in f.readlines():
                if ln.startswith('#'):
                    return ln.split(':')[-1].strip('# \n')
            return 'No description!'
