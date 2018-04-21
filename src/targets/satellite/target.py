import logging
import numpy as num

from pyrocko import gf
from pyrocko.guts import String, Bool, Dict, List

from grond.meta import Parameter
from ..base import MisfitConfig, MisfitTarget, MisfitResult, TargetGroup

guts_prefix = 'grond'
logger = logging.getLogger('grond.targets.satellite.target')


class SatelliteMisfitConfig(MisfitConfig):
    use_weight_focal = Bool.T(default=False)
    optimize_orbital_ramp = Bool.T(default=True)
    ranges = Dict.T(
        String.T(), gf.Range.T(),
        default={'offset': '-0.5 .. 0.5',
                 'ramp_east': '-1e-4 .. 1e-4',
                 'ramp_north': '-1e-4 .. 1e-4'
                 })


class SatelliteTargetGroup(TargetGroup):
    kite_scenes = List.T(optional=True)
    misfit_config = SatelliteMisfitConfig.T()

    def get_targets(self, ds, event, default_path):
        logger.debug('Selecting satellite targets...')
        targets = []

        for scene in ds.get_kite_scenes():
            if scene.meta.scene_id not in self.kite_scenes and\
               '*all' not in self.kite_scenes:
                continue

            qt = scene.quadtree

            lats = num.empty(qt.nleaves)
            lons = num.empty(qt.nleaves)
            lats.fill(qt.frame.llLat)
            lons.fill(qt.frame.llLon)

            if qt.frame.isDegree():
                logger.debug('Target %s is referenced in degree'
                             % scene.meta.scene_id)
                lons += qt.leaf_focal_points[:, 0]
                lats += qt.leaf_focal_points[:, 1]
                east_shifts = num.zeros_like(lats)
                north_shifts = num.zeros_like(lats)
            elif qt.frame.isDegree():
                logger.debug('Target %s is referenced in meter'
                             % scene.meta.scene_id)
                east_shifts = qt.leaf_focal_points[:, 0]
                north_shifts = qt.leaf_focal_points[:, 0]

            sat_target = SatelliteMisfitTarget(
                quantity='displacement',
                scene_id=scene.meta.scene_id,
                lats=lats,
                lons=lons,
                east_shifts=east_shifts,
                north_shifts=north_shifts,
                theta=qt.leaf_thetas,
                phi=qt.leaf_phis,
                tsnapshot=None,
                interpolation=self.interpolation,
                store_id=self.store_id,
                normalisation_family=self.normalisation_family,
                path=self.path or default_path,
                misfit_config=self.misfit_config)

            sat_target.set_dataset(ds)
            targets.append(sat_target)

        return targets


class SatelliteMisfitResult(gf.Result, MisfitResult):
    statics_syn = Dict.T(optional=True)
    statics_obs = Dict.T(optional=True)


class SatelliteMisfitTarget(gf.SatelliteTarget, MisfitTarget):
    scene_id = String.T()
    available_parameters = [
        Parameter('offset', 'm'),
        Parameter('ramp_north', 'm/m'),
        Parameter('ramp_east', 'm/m')
        ]
    misfit_config = SatelliteMisfitConfig.T()

    def __init__(self, *args, **kwargs):
        gf.SatelliteTarget.__init__(self, *args, **kwargs)
        MisfitTarget.__init__(self, **kwargs)
        if not self.misfit_config.optimize_orbital_ramp:
            self.parameters = []
        else:
            self.parameters = self.available_parameters

        self.parameter_values = {}

    @property
    def target_ranges(self):
        return self.misfit_config.ranges

    def string_id(self):
        return '.'.join([self.path, self.scene_id])

    def set_dataset(self, ds):
        MisfitTarget.set_dataset(self, ds)

    def post_process(self, engine, source, statics):
        scene = self._ds.get_kite_scene(self.scene_id)
        quadtree = scene.quadtree

        stat_obs = quadtree.leaf_medians

        if self.misfit_config.optimize_orbital_ramp:
            stat_level = num.zeros_like(stat_obs)
            stat_level.fill(self.parameter_values['offset'])
            stat_level += (quadtree.leaf_center_distance[:, 0]
                           * self.parameter_values['ramp_east'])
            stat_level += (quadtree.leaf_center_distance[:, 1]
                           * self.parameter_values['ramp_north'])
            statics['displacement.los'] += stat_level

        stat_syn = statics['displacement.los']

        res = stat_obs - stat_syn

        misfit_value = num.sqrt(
            num.sum((res * scene.covariance.weight_vector)**2))
        misfit_norm = num.sqrt(
            num.sum((stat_obs * scene.covariance.weight_vector)**2))

        result = SatelliteMisfitResult(
            misfits=num.array([[misfit_value, misfit_norm]], dtype=num.float))

        if self._result_mode == 'full':
            result.statics_syn = statics
            result.statics_obs = quadtree.leaf_medians

        return result

    def get_combined_weight(self, apply_balancing_weights=False):
        return num.array([self.manual_weight], dtype=num.float)

    def prepare_modelling(self, engine, source):
        return [self]

    def finalize_modelling(
            self, engine, source, modelling_targets, modelling_results):
        return modelling_results[0]

    @classmethod
    def get_plots(cls):
        from . import plot
        return [plot.SatelliteTargetPlotter]


__all__ = '''
    SatelliteTargetGroup
    SatelliteMisfitConfig
    SatelliteMisfitTarget
    SatelliteMisfitResult
'''.split()
