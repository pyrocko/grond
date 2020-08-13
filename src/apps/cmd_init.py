from __future__ import print_function, absolute_import

import logging

import glob
import os.path as op

from distutils.dir_util import copy_tree
import distutils.errors

logger = logging.getLogger('grond.init')
km = 1e3


class GrondInit(object):

    snippet_path = op.join(op.dirname(__file__), '..', 'data', 'snippets')
    example_path = op.join(op.dirname(__file__), '..', 'data', 'examples')

    def __init__(self):
        pass

    def get_examples(self):
        return {
            self.filename_to_abbrv(fn): self._get_example_description(fn)
            for fn in self.example_dirs
        }

    def get_sections(self):
        return {
            self.filename_to_abbrv(fn): self._get_description(fn)
            for fn in self.snippet_files
        }

    @property
    def example_dirs(self):
        return [path for path in glob.glob(op.join(self.example_path, '*'))
                if op.isdir(path)]

    @property
    def section_files(self):
        return self._get_snippet_files('section_*.gronf')

    @property
    def snippet_files(self):
        return self._get_snippet_files('*.gronf')

    def _get_snippet_files(self, name):
        files = glob.glob(op.join(self.snippet_path, name))
        files.sort()
        return files

    @staticmethod
    def _get_description(filename):
        with open(filename, 'rt') as f:
            for ln in f.readlines():
                if ln.startswith('#'):
                    return ln.split(':')[-1].strip('# \n')
            return 'No description!'

    def _get_example_description(self, example_dir):
        config_file = self._get_example_config(example_dir)
        with open(config_file, 'rt') as f:
            for ln in f.readlines():
                if ln.startswith('#'):
                    return ln.split(':')[-1].strip('# \n')
            return 'No description!'

    @staticmethod
    def _get_example_config(example_dir):
        fpath_template = op.join(example_dir, 'config', '*.gronf')
        config_file = glob.glob(fpath_template)
        if len(config_file) == 0:
            raise OSError('No example config file found: %s' % fpath_template)
        return config_file[0]

    @staticmethod
    def filename_to_abbrv(filename):
        return op.basename(filename).split('.')[0]

    def init_example(self, abbrv, path, force=False):

        path = op.abspath(path)
        if op.exists(path) and not force:
            raise OSError('Directory already exists: %s' % op.basename(path))
        elif op.exists(path) and force:
            pass
        example_dir = self.abbrv_to_example_dir(abbrv)

        logger.info('Initialising example "%s" in "%s".', abbrv, path)
        try:
            copy_tree(example_dir, path)
        except distutils.errors.DistutilsFileError:
            logger.error('Could not find example %s!', abbrv)

    def abbrv_to_filename(self, abbrv):
        ext = '.gronf'
        fn = op.join(self.snippet_path, abbrv + ext)

        if fn not in self._get_snippet_files('*.gronf'):
            raise OSError('File not found: %s' % fn)
        return fn

    def abbrv_to_example_dir(self, abbrv):
        return op.join(self.example_path, abbrv)

    def get_content_example(self, abbrv):
        try:
            fn = self._get_example_config(
                self.abbrv_to_example_dir(abbrv))
        except OSError:
            return False

        with open(fn, 'r') as f:
            return f.read()

    def get_content_snippet(self, abbrv):
        try:
            fn = self.abbrv_to_filename(abbrv)
        except OSError:
            return False

        with open(fn, 'r') as f:
            return f.read()
