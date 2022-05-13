from __future__ import absolute_import

import glob

from .common import grond, run_in_project
from grond import config


def test_starting_point():

    def main(env, rundir_path):

        conf = env.get_config()

        for starting_point in [
                'mean',
                'random',
                'excentricity_compensated']:

            config_path = 'config/starting_point_test.gronf'

            mod_conf = conf.clone()
            mod_conf.set_elements(
                'problem_config.name_template',
                'starting_point_test_%s_${event_name}' % starting_point)
            mod_conf.set_elements(
                'analyser_configs[:].niterations', 20)
            mod_conf.set_elements(
                'optimiser_config.sampler_phases[:].niterations', 20)
            mod_conf.set_elements(
                'optimiser_config.nbootstrap', 10)

            mod_conf.set_elements(
                'optimiser_config.sampler_phases[1].starting_point',
                starting_point)

            mod_conf.set_basepath(conf.get_basepath())
            config.write_config(mod_conf, config_path)

            grond('go', config_path)
            for rundir in glob.glob('runs/starting_point_test_*'):
                grond('plot', 'acceptance', rundir)

    run_in_project(
        main,
        project_dir_source='example_regional_cmt_full',
        project_dir='example_regional_cmt_full_starting_point',
        event_name='gfz2018pmjk',
        config_path='config/regional_cmt_ampspec.gronf')
