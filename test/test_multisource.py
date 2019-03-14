from __future__ import absolute_import, print_function
import os
import shutil

from pyrocko import guts

from grond import config
from grond import Environment, EventGroup, dump_event_groups, \
    MultiCMTProblemConfig, MultiRectangularProblemConfig

from . import common
from .common import grond, chdir

_multiprocess_can_split = True


def test_2cmt_distinct():
    playground_dir = common.get_playground_dir()
    common.get_test_data('gf_stores/crust2_ib/')
    common.get_test_data('gf_stores/crust2_ib_static/')
    gf_stores_path = common.test_data_path('gf_stores')

    with chdir(playground_dir):
        scenario_dir = 'scenario_multisource'
        if os.path.exists(scenario_dir):
            shutil.rmtree(scenario_dir)

        grond('scenario', '--targets=waveforms', '--nevents=2',
              '--nstations=5', '--gf-store-superdirs=%s' % gf_stores_path,
              scenario_dir)

        with chdir(scenario_dir):
            config_path = 'config/scenario.gronf'
            quick_config_path = 'config/scenario_quick_multi.gronf'
            event_names = grond('events', config_path).strip().split('\n')
            event_group_name = 'group'

            eg = EventGroup(
                name=event_group_name,
                event_names=event_names)

            event_groups_path = 'event_group.yaml'
            dump_event_groups([eg], event_groups_path)

            env = Environment([config_path] + event_names)
            conf = env.get_config()

            mod_conf = conf.clone()
            mod_conf.dataset_config.event_groups_path = event_groups_path
            mod_conf.set_elements(
                'analyser_configs[:].niterations', 100)
            mod_conf.set_elements(
                'optimiser_config.sampler_phases[:].niterations', 100)
            mod_conf.set_elements(
                'optimiser_config.nbootstrap', 5)

            kwargs = guts.to_dict(mod_conf.problem_config)
            mod_conf.problem_config = MultiCMTProblemConfig(**kwargs)

            mod_conf.set_basepath(conf.get_basepath())
            config.write_config(mod_conf, quick_config_path)
            grond('diff', config_path, quick_config_path)
            grond('check', quick_config_path, event_group_name)

            grond('go', '--status=quiet', quick_config_path, event_group_name)
            rundir_paths = common.get_rundir_paths(
                quick_config_path, [event_group_name])

            grond('report', *rundir_paths)
            grond('report', '--parallel=2', *rundir_paths)


def test_2rectangular_distinct():
    playground_dir = common.get_playground_dir()
    common.get_test_data('gf_stores/crust2_ib/')
    common.get_test_data('gf_stores/crust2_ib_static/')
    gf_stores_path = common.test_data_path('gf_stores')

    with chdir(playground_dir):
        scenario_dir = 'scenario_multisource'
        if os.path.exists(scenario_dir):
            shutil.rmtree(scenario_dir)

        grond('scenario', '--targets=waveforms', '--nevents=2',
              '--nstations=5', '--problem=rectangular',
              '--gf-store-superdirs=%s' % gf_stores_path,
              scenario_dir)

        with chdir(scenario_dir):
            config_path = 'config/scenario.gronf'
            quick_config_path = 'config/scenario_quick_multi.gronf'
            event_names = grond('events', config_path).strip().split('\n')
            event_group_name = 'group'

            eg = EventGroup(
                name=event_group_name,
                event_names=event_names)

            event_groups_path = 'event_group.yaml'
            dump_event_groups([eg], event_groups_path)

            env = Environment([config_path] + event_names)
            conf = env.get_config()

            mod_conf = conf.clone()
            mod_conf.dataset_config.event_groups_path = event_groups_path
            mod_conf.set_elements(
                'analyser_configs[:].niterations', 100)
            mod_conf.set_elements(
                'optimiser_config.sampler_phases[:].niterations', 100)
            mod_conf.set_elements(
                'optimiser_config.nbootstrap', 5)

            kwargs = guts.to_dict(mod_conf.problem_config)
            mod_conf.problem_config = MultiRectangularProblemConfig(**kwargs)

            mod_conf.set_basepath(conf.get_basepath())
            config.write_config(mod_conf, quick_config_path)
            grond('diff', config_path, quick_config_path)
            grond('check', quick_config_path, event_group_name)

            grond('go', '--status=quiet', quick_config_path, event_group_name)
            rundir_paths = common.get_rundir_paths(
                quick_config_path, [event_group_name])

            grond('report', *rundir_paths)
            grond('report', '--parallel=2', *rundir_paths)
