from __future__ import absolute_import
import os
import shutil

from pyrocko import model

from grond import config
from grond import Environment
from grond.meta import expand_template

from . import common
from .common import grond, chdir

_multiprocess_can_split = True


def test_example_regional_cmt():
    playground_dir = common.get_playground_dir()
    with chdir(playground_dir):
        run_example(
            'example_regional_cmt',
            'config/regional_cmt.gronf',
            'config/regional_cmt_quick.gronf',
            'gfz2018pmjk')


def test_example_wphase():
    playground_dir = common.get_playground_dir()
    with chdir(playground_dir):
        run_example(
            'example_wphase',
            'config/wphase_cmt.gronf',
            'config/wphase_cmt_quick.gronf',
            'gfz2015sfdd')


def test_example_insar():
    playground_dir = common.get_playground_dir()
    with chdir(playground_dir):
        run_example(
            'example_insar',
            'config/insar_rectangular.gronf',
            'config/insar_rectangular_quick.gronf',
            '2009laquila')


def run_example(project_name, config_path, quick_config_path, event_name):
    project_dir = project_name
    if os.path.exists(project_dir):
        shutil.rmtree(project_dir)

    grond('init', project_name, project_dir)
    with chdir(project_dir):
        assert os.path.isdir('config')

        common.link_test_data(
            'events/%s/' % event_name, 'data/events/%s/' % event_name)

        env = Environment([config_path, event_name])
        conf = env.get_config()

        store_ids = conf.get_elements('target_groups[:].store_id')
        for store_id in store_ids:
            store_path = 'gf_stores/%s/' % store_id
            if not os.path.exists(store_path):
                common.link_test_data(store_path)

        problem_name = env.get_problem().name
        rundir_path = expand_template(
            conf.rundir_template,
            dict(problem_name=problem_name))

        grond(
            'check', config_path, event_name,
            '--save-stations-used=used_stations.txt')

        sorted(
            s.station for s in model.load_stations('used_stations.txt'))

        mod_conf = conf.clone()
        mod_conf.set_elements(
            'analyser_configs[:].niterations', 100)
        mod_conf.set_elements(
            'optimiser_config.sampler_phases[:].niterations', 100)
        mod_conf.set_elements(
            'optimiser_config.nbootstrap', 10)
        mod_conf.set_basepath(conf.get_basepath())
        config.write_config(mod_conf, quick_config_path)
        grond('go', quick_config_path, event_name)
        grond('harvest', '--force', '--export-fits=best,mean', rundir_path)
        grond('report', rundir_path)
        # assert os.path.isdir(os.path.join('report', event_name))
