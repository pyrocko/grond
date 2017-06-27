import logging
import numpy as num

from pyrocko import gf
from pyrocko.guts import String, Bool, Dict, List

from .base import MisfitTarget, MisfitConfig, MisfitResult, TargetGroup
from ..meta import Parameter

guts_prefix = 'grond'
logger = logging.getLogger('grond.target').getChild('satellite')


class SatelliteTargetGroup(TargetGroup):
    kite_scenes = List.T(optional=True)

    def get_targets(self, ds, event, default_path):
        logger.debug('Selecting satellite targets...')
        targets = []

        for scene in ds.get_kite_scenes():
            if scene.meta.scene_id not in self.kite_scenes and\
               '*all' not in self.kite_scenes:
                continue

            if not isinstance(self.misfit_config,
                              SatelliteMisfitConfig):
                raise AttributeError('misfit_config must be of type'
                                     ' SatelliteMisfitConfig')

            qt = scene.quadtree

            lats = num.empty(qt.nleaves)
            lons = num.empty(qt.nleaves)
            lats.fill(qt.frame.llLat)
            lons.fill(qt.frame.llLon)

            north_shifts = qt.leaf_focal_points[:, 1]
            east_shifts = qt.leaf_focal_points[:, 0]

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


class SatelliteMisfitResult(MisfitResult):
    statics_syn = Dict.T(optional=True)
    statics_obs = Dict.T(optional=True)


class SatelliteMisfitConfig(MisfitConfig):
    use_weight_focal = Bool.T(default=False)
    optimize_orbital_ramp = Bool.T(default=True)
    ranges = Dict.T(String.T(), gf.Range.T(),
                    default={'offset': '-0.5 .. 0.5',
                             'ramp_north': '-1e-4 .. 1e-4',
                             'ramp_east': '-1e-4 .. 1e-4'})


class SatelliteMisfitTarget(gf.SatelliteTarget, MisfitTarget):
    scene_id = String.T()

    available_parameters = [
        Parameter('offset', 'm'),
        Parameter('ramp_north', 'm/m'),
        Parameter('ramp_east', 'm/m'),
        ]

    def __init__(self, *args, **kwargs):
        gf.SatelliteTarget.__init__(self, *args, **kwargs)
        MisfitTarget.__init__(self)
        if not self.misfit_config.optimize_orbital_ramp:
            self.parameters = []
        else:
            self.parameters = self.available_parameters

        self.parameter_values = {}

    @property
    def id(self):
        return self.scene_id

    def set_dataset(self, ds):
        MisfitTarget.set_dataset(self, ds)
        scene = self._ds.get_kite_scene(self.scene_id)
        self.nmisfits = scene.quadtree.nleaves

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
            misfit_value=misfit_value,
            misfit_norm=misfit_norm)

        if self._result_mode == 'full':
            result.statics_syn = statics
            result.statics_obs = quadtree.leaf_medians

        return result

    def get_combined_weight(self, apply_balancing_weights=False):
        return self.manual_weight
