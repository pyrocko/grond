import logging
import os
import shutil
import os.path as op

import grond

from pyrocko import gf, scenario, util


DEFAULT_STATIC_STORE = 'ak135_static'
DEFAULT_WAVEFORM_STORE = 'global2s'

logger = logging.getLogger('grond.scenario')
km = 1e3


class GrondScenario(object):

    def __init__(self, project_dir,
                 center_lat=23., center_lon=52., radius=230*km,
                 problem=None, observations=[]):
        self.project_dir = project_dir
        self.data_dir = op.join('events', 'scenario')

        self.center_lat = center_lat
        self.center_lon = center_lon
        self.radius = radius

        self.problem = problem
        self.observations = observations

        self.engine = gf.LocalEngine(use_config=True)

    def create_project_dir(self, force=False):
        prj_dir = self.project_dir

        if op.exists(prj_dir) and not force:
            raise EnvironmentError(
                'Directory %s already exists! Use force to overwrite'
                % prj_dir)
        elif op.exists(prj_dir) and force:
            logger.info('Overwriting directory %s.' % prj_dir)
            shutil.rmtree(prj_dir)
        util.ensuredir(prj_dir)

    @property
    def stores_wanted(self):
        return set([obs.store_id for obs in self.observations])

    def symlink_gfstores(self):
        logger.info('Symlinking Green\'s function stores...')
        util.ensuredir(op.join(self.project_dir, 'gf_stores'))

        for store_id in self.stores_wanted:
            store = self.engine.get_store(store_id)
            os.symlink(store.store_dir,
                       op.join(self.project_dir, 'gf_stores', store_id))

    def add_observation(self, observation):
        logger.info('Adding %s' % observation.__class__.__name__)
        self.observations.append(observation)

    def set_problem(self, problem):
        self.problem = problem

    def get_dataset_config(self):
        dataset_config = grond.DatasetConfig(
            events_path=op.join(self.data_dir, 'events.txt'))
        for obs in self.observations:
            obs.update_dataset_config(dataset_config, self.data_dir)
        return dataset_config

    def get_scenario(self):
        if not self.observations:
            raise AttributeError('No observations set,'
                                 ' use .add_observation(Observation)'
                                 ' to add one!')
        if not self.problem:
            raise AttributeError('No Source Problem set'
                                 ' use .set_problem(Problem) to set one')

        return scenario.ScenarioGenerator(
            center_lat=self.center_lat,
            center_lon=self.center_lon,
            radius=self.radius,
            target_generators=[obs.get_scenario_target_generator()
                               for obs in self.observations],
            source_generator=self.problem.get_scenario_source_generator())

    def create_scenario(self, interactive=True):
        logger.info('Creating pyrocko.scenario...')

        scenario = self.get_scenario()
        scenario.init_modelling(engine=self.engine)
        scenario.ensure_gfstores(interactive=interactive)
        scenario.dump(filename=op.join(self.project_dir, 'scenario.yml'))

        data_dir = op.join(self.project_dir, self.data_dir)
        util.ensuredir(data_dir)

        scenario.dump_data(path=data_dir)
        scenario.make_map(op.join(self.project_dir, 'scenario_map.pdf'))

        shutil.move(op.join(data_dir, 'sources.yml'),
                    op.join(data_dir, 'scenario_sources.yml'))

    def get_grond_config(self):
        engine_config = grond.EngineConfig(
            gf_stores_from_pyrocko_config=False,
            gf_store_superdirs=['gf_stores'])

        optimiser_config = grond.HighScoreOptimiserConfig()

        config = grond.Config(
            rundir_template=op.join('rundir', '${problem_name}.grun'),
            dataset_config=self.get_dataset_config(),
            target_groups=[obs.get_grond_target_group()
                           for obs in self.observations],
            problem_config=self.problem.get_grond_problem_config(),
            optimiser_config=optimiser_config,
            engine_config=engine_config)

        return config

    def create_grond_files(self):
        logger.info('Creating grond configuration for %s'
                    % ' and '.join([obs.name for obs in self.observations]))

        with open(op.join(self.project_dir, 'config.yml'), 'w') as cf:
            cf.write(str(self.get_grond_config()))

    def build(self, force=False, interactive=False):
        logger.info('Building Grond scenario using pyrocko.scenario...')

        self.create_project_dir(force)

        self.create_scenario(interactive=interactive)
        self.symlink_gfstores()

        self.create_grond_files()


class Observation(object):

    name = 'Observation Prototype'

    def __init__(self, store_id, *args, **kwargs):
        self.store_id = store_id

    def update_dataset_config(self, dataset, data_dir):
        return dataset

    def get_scenario_target_generator(self):
        pass

    def get_grond_target_group(self):
        pass


class WaveformObservation(Observation):

    name = 'seismic waveforms'

    def __init__(self, store_id=DEFAULT_WAVEFORM_STORE, nstations=25):
        self.nstations = nstations
        self.store_id = store_id

    def update_dataset_config(self, dataset_config, data_dir):
        ds = dataset_config
        ds.waveform_paths = [op.join(data_dir, 'data/waveforms')]
        ds.waveform_path = op.join(data_dir, 'data/waveforms')
        ds.stations_path = op.join(data_dir, 'meta', 'stations.txt')
        ds.responses_stationxml_paths = [
            op.join(data_dir, 'meta', 'stations.xml')]
        return ds

    def get_scenario_target_generator(self):
        return scenario.targets.WaveformGenerator(
            station_generator=scenario.targets.RandomStationGenerator(
                nstations=self.nstations),
            store_id=self.store_id,
            seismogram_quantity='displacement')

    def get_grond_target_group(self):
        return grond.WaveformTargetGroup(
            normalisation_family='time_domain',
            path='all',
            distance_min=10*km,
            distance_max=1000*km,
            channels=['Z', 'R', 'T'],
            interpolation='multilinear',
            store_id=self.store_id,
            misfit_config=grond.WaveformMisfitConfig(
                fmin=0.01,
                fmax=0.1))


class InSARObservation(Observation):

    name = 'InSAR displacement'

    def __init__(self, store_id=DEFAULT_STATIC_STORE):
        self.store_id = store_id

    def update_dataset_config(self, dataset_config, data_dir):
        dataset_config.kite_scene_paths = [op.join(data_dir, 'insar')]
        return dataset_config

    def get_scenario_target_generator(self):
        logger.warning('Inspect the InSAR scenes, covariance and quadtree!')
        return scenario.targets.InSARGenerator(
                store_id=self.store_id)

    def get_grond_target_group(self):
        return grond.SatelliteTargetGroup(
            normalisation_family='insar_target',
            path='all',
            interpolation='multilinear',
            store_id=self.store_id,
            kite_scenes=['*all'],
            misfit_config=grond.SatelliteMisfitConfig(
                use_weight_focal=False,
                optimise_orbital_ramp=True,
                ranges={
                    'offset': '-0.5 .. 0.5',
                    'ramp_north': '-1e-4 .. 1e-4',
                    'ramp_east': '-1e-4 .. 1e-4'
                    })
            )


class GNSSCampaignObservation(Observation):

    name = 'GNSS campaign'

    def __init__(self, store_id=DEFAULT_STATIC_STORE, nstations=25):
        self.store_id = store_id
        self.nstations = nstations

    def update_dataset_config(self, dataset_config, data_dir):
        dataset_config.gnss_campaign_paths = [op.join(data_dir, 'gnss')]
        return dataset_config

    def get_scenario_target_generator(self):
        return scenario.targets.GNSSCampaignGenerator(
            station_generator=scenario.targets.RandomStationGenerator(
                nstations=self.nstations),
            store_id=self.store_id)

    def get_grond_target_group(self):
        return grond.GNSSCampaignTargetGroup(
            gnss_campaigns=['*all'],
            normalisation_family='gnss_target',
            path='all',
            interpolation='multilinear',
            store_id=self.store_id,
            misfit_config=grond.GNSSCampaignMisfitConfig()
        )


class SourceProblem(object):

    def __init__(self, magnitude_min=6., magnitude_max=7., nevents=1):
        self.magnitude_min = magnitude_min
        self.magnitude_max = magnitude_max
        self.nevents = nevents

    def get_scenario_source_generator(self):
        return

    def get_grond_problem_config(self):
        return


class DCSourceProblem(SourceProblem):

    def get_scenario_source_generator(self):
        return scenario.sources.DCSourceGenerator(
            magnitude_min=self.magnitude_min,
            magnitude_max=self.magnitude_max,
            nevents=self.nevents)

    def get_grond_problem_config(self):
        import math

        pi2 = math.pi/2
        return grond.CMTProblemConfig(
            name_template='cmt-${event_name}',
            distance_min=2.*km,
            mt_type='deviatoric',
            ranges=dict(
                time=gf.Range(0, 10.0, relative='add'),
                north_shift=gf.Range(-16*km, 16*km),
                east_shift=gf.Range(-16*km, 16*km),
                depth=gf.Range(1*km, 11*km),
                magnitude=gf.Range(4.0, 6.0),
                rmnn=gf.Range(-pi2, pi2),
                rmee=gf.Range(-pi2, pi2),
                rmdd=gf.Range(-pi2, pi2),
                rmne=gf.Range(-1.0, 1.0),
                rmnd=gf.Range(-1.0, 1.0),
                rmed=gf.Range(-1.0, 1.0),
                duration=gf.Range(1.0, 15.0))
            )


class RectangularSourceProblem(SourceProblem):

    def get_scenario_source_generator(self):
        return scenario.sources.RectangularSourceGenerator(
            magnitude_min=self.magnitude_min,
            magnitude_max=self.magnitude_max,
            nevents=self.nevents)

    def get_grond_problem_config(self):
        return grond.RectangularProblemConfig(
            name_template='rect-source_${event_name}',
            decimation_factor=8,
            ranges=dict(
                north_shift=gf.Range(-20*km, 20*km),
                east_shift=gf.Range(-20*km, 20*km),
                depth=gf.Range(0*km, 10*km),
                length=gf.Range(20*km, 40*km),
                width=gf.Range(5*km, 12*km),
                dip=gf.Range(20, 70),
                strike=gf.Range(0, 180),
                rake=gf.Range(0, 90),
                slip=gf.Range(1, 3))
            )
