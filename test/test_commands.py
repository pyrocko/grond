from __future__ import absolute_import

import re
import tempfile

from pyrocko import model

from . import common
from .common import grond, run_in_project
from grond import config


def test_usage():
    common.assert_grond_usage()


def test_version():
    assert grond('version').startswith('--- !grond.VersionInfo')
    assert re.match(r'^\d+\.\d+\.\d+$', grond('version', '--short'))
    assert all(len(x.split(': ')) == 2
               for x in grond('version', '--failsafe').split('\n') if x)


def test_init():
    common.assert_grond_usage('init')
    assert grond('init', 'list').startswith('Available')
    with tempfile.NamedTemporaryFile() as f:
        f.write(grond('init', 'example_regional_cmt').encode('utf8'))
        f.flush()
        config.read_config(f.name)


def test_export():

    def main(env, rundir_path):

        def export_test(what, nevents):
            with tempfile.NamedTemporaryFile() as f:
                f.write(grond('export', 'best', rundir_path).encode('utf8'))
                f.flush()
                events = model.load_events(f.name)
                assert len(events) == 1

        export_test('best', 1)
        export_test('mean', 1)
        export_test('ensemble', 1010)

    run_in_project(
        main,
        project_dir_source='example_regional_cmt_full',
        project_dir='example_regional_cmt_full_export',
        event_name='gfz2018pmjk',
        config_path='config/regional_cmt_ampspec.gronf')
