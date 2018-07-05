import grond
import logging

import math
import shutil
import os
from os import path as op
from pyrocko.gf import Range

logger = logging.getLogger('grond.init')
km = 1e3

INIT_EVENT = '''name = 2011-myanmar
time = 2011-03-24 13:55:12.010
latitude = 20.687
longitude = 99.822
magnitude = 6.9
moment = 1.9228e+19
depth = 8000
region = Myanmar
--------------------------------------------'''


class GrondProject(object):

    def __init__(self):
        self.dataset_config = grond.DatasetConfig(
            events_path='events.txt',
            stations_path='stations.txt',
            responses_stationxml_paths=['stationxml.xml'])

        self.project_config = None
        self.target_groups = []

        self.sub_dirs = ['gf_stores', 'config']
        self.empty_files = []

    def add_waveforms(self):
        logger.info('Added waveforms')
        dsc = self.dataset_config

        self.sub_dirs.append('data')
        dsc.waveform_paths = ['data']

        self.empty_files.append('stations.xml')
        dsc.stations_path = 'stations.txt'

        self.target_groups.append(
            grond.WaveformTargetGroup(
                normalisation_family='time_domain',
                path='all',
                distance_min=10*km,
                distance_max=1000*km,
                channels=['Z', 'R', 'T'],
                interpolation='multilinear',
                store_id='gf_store_id',
                misfit_config=grond.WaveformMisfitConfig(
                    fmin=0.01,
                    fmax=0.1)))

    def add_insar(self):
        logger.info('Added InSAR')
        dsc = self.dataset_config

        self.sub_dirs.append('scenes')
        dsc.kite_scene_paths = ['scenes']

        self.target_groups.append(
            grond.SatelliteTargetGroup(
                normalisation_family='static',
                path='all',
                interpolation='multilinear',
                store_id='gf_static_store_id',
                kite_scenes=['*all'],
                misfit_config=grond.SatelliteMisfitConfig(
                    optimise_orbital_ramp=True,
                    ranges={
                        'offset': '-0.5 .. 0.5',
                        'ramp_north': '-1e-4 .. 1e-4',
                        'ramp_east': '-1e-4 .. 1e-4'
                        }
                    ))
            )

    def add_gnss(self):
        logger.info('Added InSAR')
        dsc = self.dataset_config

        self.sub_dirs.append('scenes')
        dsc.kite_scene_paths = ['scenes']

        self.target_groups.append(
            grond.GNSSCampaignTargetGroup(
                normalisation_family='static',
                path='all',
                interpolation='multilinear',
                store_id='gf_static_store_id',
                misfit_config=grond.GNSSCampaignMisfitConfig()
                )
            )

    def set_cmt_source(self):
        logger.info('Set problem to rectangular')
        pi2 = math.pi/2
        self.problem_config = grond.CMTProblemConfig(
            name_template='cmt_${event_name}',
            distance_min=2.*km,
            mt_type='deviatoric',
            ranges=dict(
                time=Range(0, 10.0, relative='add'),
                north_shift=Range(-16*km, 16*km),
                east_shift=Range(-16*km, 16*km),
                depth=Range(1*km, 11*km),
                magnitude=Range(4.0, 6.0),
                rmnn=Range(-pi2, pi2),
                rmee=Range(-pi2, pi2),
                rmdd=Range(-pi2, pi2),
                rmne=Range(-1.0, 1.0),
                rmnd=Range(-1.0, 1.0),
                rmed=Range(-1.0, 1.0),
                duration=Range(1.0, 15.0))
            )

    def set_rectangular_source(self):
        logger.info('Set problem to rectangular')
        self.problem_config = grond.RectangularProblemConfig(
            name_template='rect_source',
            ranges=dict(
                north_shift=Range(-20*km, 20*km),
                east_shift=Range(-20*km, 20*km),
                depth=Range(0*km, 10*km),
                length=Range(20*km, 40*km),
                width=Range(5*km, 12*km),
                dip=Range(20, 70),
                strike=Range(0, 180),
                rake=Range(0, 90),
                slip=Range(1, 3))
            )

    def check(self):
        assert self.target_groups, 'No target_groups set!'
        assert self.problem_config, 'No problem set!'

    def get_config(self):
        self.check()

        engine_config = grond.EngineConfig(
            gf_store_superdirs=['gf_stores'])

        optimiser_config = grond.HighScoreOptimiserConfig()

        return grond.Config(
            rundir_template=op.join('rundir', '${problem_name}.grun'),
            dataset_config=self.dataset_config,
            problem_config=self.problem_config,
            target_groups=self.target_groups,
            optimiser_config=optimiser_config,
            engine_config=engine_config)

    def build(self, project_dir, force=False):
        if op.exists(project_dir) and not force:
            raise grond.GrondError(
                'Directory "%s" already exists! Use --force to overwrite'
                % project_dir)
        elif op.exists(project_dir) and force:
            logger.info('Overwriting directory %s.' % project_dir)
            shutil.rmtree(project_dir)

        logger.info('Creating empty project in folder %s' % project_dir)
        config = self.get_config()

        def p(*fn):
            return op.join(project_dir, *fn)

        os.mkdir(op.abspath(project_dir))
        for d in self.sub_dirs:
            os.mkdir(p(d))

        with open(p('config', 'config.gronf'), 'w') as f:
            f.write(str(config))
        with open(p('events.txt'), 'w') as f:
            f.write(INIT_EVENT)

        for fn in self.empty_files:
            open(p(fn), 'w').close()

    def dump(self):
        return self.get_config().dump()
