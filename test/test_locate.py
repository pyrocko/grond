from __future__ import absolute_import
import os
import shutil

from pyrocko import guts

from grond import config
from grond import Environment

from . import common
from .common import grond, chdir

_multiprocess_can_split = True


def test_locate():
    playground_dir = common.get_playground_dir()
    common.get_test_data('gf_stores/crust2_ib/')
    gf_stores_path = common.test_data_path('gf_stores')

    with chdir(playground_dir):
        scenario_dir = 'scenario_locate'
        if os.path.exists(scenario_dir):
            shutil.rmtree(scenario_dir)

        grond('scenario', '--targets=waveforms', '--nevents=1',
              '--nstations=6', '--gf-store-superdirs=%s' % gf_stores_path,
              '--no-map',
              scenario_dir)

        with chdir(scenario_dir):
            config_path = 'config/scenario.gronf'
            quick_config_path = 'config/scenario_quick.gronf'
            event_names = grond('events', config_path).strip().split('\n')

            env = Environment([config_path] + event_names)
            conf = env.get_config()

            mod_conf = conf.clone()
            target_groups = guts.load_all(string='''
--- !grond.PhasePickTargetGroup
normalisation_family: picks
path: pick.p
weight: 1.0
store_id: crust2_ib
distance_min: 10000.0
distance_max: 1000000.0
pick_synthetic_traveltime: '{stored:any_P}'
pick_phasename: 'any_P'
--- !grond.PhasePickTargetGroup
normalisation_family: picks
path: pick.s
weight: 1.0
store_id: crust2_ib
distance_min: 10000.0
distance_max: 1000000.0
pick_synthetic_traveltime: '{stored:any_S}'
pick_phasename: 'any_S'
''')

            mod_conf.dataset_config.picks_paths = [
                'data/scenario/picks/picks.markers']

            mod_conf.target_groups = target_groups

            # mod_conf.set_elements(
            #     'analyser_configs[:].niterations', 100)
            # mod_conf.set_elements(
            #     'optimiser_config.sampler_phases[:].niterations', 100)
            # mod_conf.set_elements(
            #     'optimiser_config.nbootstrap', 5)

            mod_conf.optimiser_config.sampler_phases[-1].niterations = 15000

            mod_conf.set_basepath(conf.get_basepath())
            config.write_config(mod_conf, quick_config_path)

            grond('diff', config_path, quick_config_path)
            grond('check', quick_config_path, *event_names)

            grond('go', quick_config_path, *event_names)
            rundir_paths = common.get_rundir_paths(config_path, event_names)
            grond('report', *rundir_paths)
