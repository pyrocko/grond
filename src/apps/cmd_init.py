import logging

import glob
import shutil
import os.path as op
from pyrocko import util

logger = logging.getLogger('grond.init')
km = 1e3


class GrondInit(object):

    folder = op.join(op.dirname(__file__), '..', 'data', 'init')

    def __init__(self):
        pass

    def get_examples(self):
        return {
            self.filename_to_abbrv(fn): self._get_description(fn)
            for fn in self.example_files
        }

    def get_sections(self):
        return {
            self.filename_to_abbrv(fn): self._get_description(fn)
            for fn in self.section_files
        }

    @property
    def example_files(self):
        return self._get_files('example_*.gronf')

    @property
    def section_files(self):
        return self._get_files('section_*.gronf')

    @property
    def snippet_files(self):
        return self._get_files('snippet_*.gronf')

    def _get_files(self, name):
        return glob.glob(op.join(self.folder, name))

    @staticmethod
    def _get_description(filename):
        with open(filename, 'rt') as f:
            for ln in f.readlines():
                if ln.startswith('#'):
                    return ln.split(':')[-1].strip('# \n')
            return 'No description!'

    @staticmethod
    def filename_to_abbrv(filename):
        return op.basename(filename).split('.')[0]

    def init_example(self, abbrv, path, force=False):
        path = op.abspath(path)
        if op.exists(path) and not force:
            raise OSError('Directory already exists: %s' % op.basename(path))

        logger.info('Initialising example configuration for %s...' % abbrv)
        sub_dirs = ['config', 'gf_stores', 'data']

        for d in sub_dirs:
            p = op.join(path, d + '/')
            logger.debug('Creating directory: %s' % p)
            util.ensuredirs(p)

        fn = self.abbrv_to_filename(abbrv)
        logger.debug('Copying config file: %s' % fn)
        shutil.copyfile(fn, op.join(path, 'config', op.basename(fn)))
        logger.debug('Copying README.md: %s' % fn)
        shutil.copyfile(op.join(self.folder, 'README_example.md'),
                        op.join(path, 'README.md'))

    def abbrv_to_filename(self, abbrv):
        ext = '.gronf'
        fn = op.join(self.folder, abbrv + ext)

        if fn not in self._get_files('*.gronf'):
            raise OSError('File not found: %s' % fn)
        return fn

    def get_content(self, abbrv):
        try:
            fn = self.abbrv_to_filename(abbrv)
        except OSError:
            return False

        with open(fn, 'r') as f:
            return f.read()
