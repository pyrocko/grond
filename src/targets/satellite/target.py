import logging
import numpy as num
from scipy import linalg as splinalg

from pyrocko import gf
from pyrocko.guts import String, Bool, Dict, List
from pyrocko.guts_array import Array

import os

from grond.meta import Parameter, has_get_plot_classes
from ..base import MisfitConfig, MisfitTarget, MisfitResult, TargetGroup

guts_prefix = 'grond'
logger = logging.getLogger('grond.targets.satellite.target')


class SatelliteMisfitConfig(MisfitConfig):
    """Carries the misfit configuration."""
    optimise_orbital_ramp = Bool.T(
        default=True,
        help='Switch to account for a linear orbital ramp or not')
    ranges = Dict.T(
        String.T(), gf.Range.T(),
        default={'offset': '-0.5 .. 0.5',
                 'ramp_east': '-1e-7 .. 1e-7',
                 'ramp_north': '-1e-7 .. 1e-7'
                 },
        help='These parameters give bounds for an offset [m], a linear'
             ' gradient in east direction [m/m] and a linear gradient in north'
             ' direction [m/m]. Note, while the optimisation of these ramps'
             ' is individual for each target, the ranges set here are common'
             ' for all satellite targets.')


class SatelliteTargetGroup(TargetGroup):
    r"""Handles maps of static ground motion from satellite observations (InSAR)

    The InSAR displacement maps post-processed by the `pyrocko` module `kite`
    are usually `Quadtree` downsampled (Jonsson, 2002). Each data point has a
    latitude, longitude, Line-of-sight displacement value [m] as well as an
    orientation and elevation angle, which define the Line-of-Sight. The data
    are associated with a weight matrix, which is the inverse of a full
    variance-covariance matrix (Sudhaus \& Jonsson, 2009). In principle, these
    data sets could stem from pixel offset maps. See also the documentation of
    the `kite` module.
    """
    kite_scenes = List.T(
        optional=True,
        help='List of InSAR data files prepared \
              by the ``pyrocko`` module ``kite``')
    misfit_config = SatelliteMisfitConfig.T(
        help='Settings for the objective function of these targets')

    def get_targets(self, ds, event, default_path='none'):
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
                logger.debug('Target "%s" is referenced in degree.'
                             % scene.meta.scene_id)
                lons += qt.leaf_focal_points[:, 0]
                lats += qt.leaf_focal_points[:, 1]
                east_shifts = num.zeros_like(lats)
                north_shifts = num.zeros_like(lats)
            elif qt.frame.isMeter():
                logger.debug('Target "%s" is referenced in meter.'
                             % scene.meta.scene_id)
                east_shifts = qt.leaf_focal_points[:, 0]
                north_shifts = qt.leaf_focal_points[:, 1]
            else:
                assert False

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
    """Carries the observations for a target and corresponding synthetics."""
    statics_syn = Dict.T(
        String.T(),
        Array.T(dtype=num.float, shape=(None,), serialize_as='base64'),
        optional=True,
        help='Predicted static displacements for a target (synthetics).')
    statics_obs = Dict.T(
        String.T(),
        Array.T(dtype=num.float, shape=(None,), serialize_as='base64'),
        optional=True,
        help='Observed static displacement for a target.')


@has_get_plot_classes
class SatelliteMisfitTarget(gf.SatelliteTarget, MisfitTarget):
    """Handles and carries out operations related to the objective functions.

    Standard operations are the calculation of the weighted misfit between
    observed and predicted (synthetic) data. If enabled in the misfit
    configuration, orbital ramps are optimized for.
    """
    scene_id = String.T(
        help='UID string each Kite displacemente scene.'
             ' Corresponds to Kite scene_id.')
    misfit_config = SatelliteMisfitConfig.T(
        help='Configuration of the ``SatelliteTarget``')

    can_bootstrap_residuals = True

    available_parameters = [
        Parameter('offset', 'm'),
        Parameter('ramp_north', 'm/m'),
        Parameter('ramp_east', 'm/m')]

    def __init__(self, *args, **kwargs):
        gf.SatelliteTarget.__init__(self, *args, **kwargs)
        MisfitTarget.__init__(self, **kwargs)
        if not self.misfit_config.optimise_orbital_ramp:
            self.parameters = []
        else:
            self.parameters = self.available_parameters

        self.parameter_values = {}

        self._noise_weight_matrix = None

    @property
    def target_ranges(self):
        return self.misfit_config.ranges

    def string_id(self):
        return '.'.join([self.path, self.scene_id])

    @property
    def id(self):
        return self.scene_id

    def set_dataset(self, ds):
        MisfitTarget.set_dataset(self, ds)

    @property
    def nmisfits(self):
        return self.lats.size

    def get_correlated_weights(self, nthreads=0):
        ''' is for L2-norm weighting, the square-rooted, inverse covar '''
        if self._noise_weight_matrix is None:
            logger.info(
                'Inverting scene covariance matrix (nthreads=%i)...'
                % nthreads)
            cov = self.scene.covariance
            cov.nthreads = nthreads

            self._noise_weight_matrix = splinalg.sqrtm(
                num.linalg.inv(cov.covariance_matrix))

            logger.info('Inverting scene covariance matrix done.')

        return self._noise_weight_matrix

    @property
    def scene(self):
        return self._ds.get_kite_scene(self.scene_id)

    def post_process(self, engine, source, statics):
        """Applies the objective function.

        As a result the weighted misfits are given and the observed and
        synthetic data. For the satellite target the orbital ramp is
        calculated and applied here."""
        scene = self.scene
        quadtree = scene.quadtree

        obs = quadtree.leaf_medians

        if self.misfit_config.optimise_orbital_ramp:
            stat_level = num.full_like(obs, self.parameter_values['offset'])

            stat_level += (quadtree.leaf_center_distance[:, 0]
                           * self.parameter_values['ramp_east'])
            stat_level += (quadtree.leaf_center_distance[:, 1]
                           * self.parameter_values['ramp_north'])
            statics['displacement.los'] += stat_level

        stat_syn = statics['displacement.los']

        res = obs - stat_syn

        misfit_value = res
        misfit_norm = obs

        mf = num.vstack([misfit_value, misfit_norm]).T
        result = SatelliteMisfitResult(
            misfits=mf)

        if self._result_mode == 'full':
            result.statics_syn = statics
            result.statics_obs = {'displacement.los': obs}

        return result

    def get_combined_weight(self):
        if self._combined_weight is None:
            # invcov = self.scene.covariance.weight_matrix
            # self._combined_weight = invcov * self.manual_weight
            self._combined_weight = num.full(self.nmisfits, self.manual_weight)

        return self._combined_weight

    def prepare_modelling(self, engine, source, targets):
        return [self]

    def finalize_modelling(
            self, engine, source, modelling_targets, modelling_results):
        return modelling_results[0]

    def init_bootstrap_residuals(self, nbootstraps, rstate=None, nthreads=0):
        logger.info(
            'Scene "%s", initializing bootstrapping residuals from noise '
            'pertubations...' % self.scene_id)

        if rstate is None:
            rstate = num.random.RandomState()

        scene = self.scene
        qt = scene.quadtree
        cov = scene.covariance
        bootstraps = num.zeros((nbootstraps, qt.nleaves))

        try:
            # TODO:mi Signal handler is not given back to the main task!
            # This is a python3.7 bug
            logger.warning('Using multi-threading for SatelliteTargets. '
                           'Python 3.7 is buggy and needs to be killed hard:'
                           ' `killall grond`')
            from concurrent.futures import ThreadPoolExecutor
            nthreads = os.cpu_count() if not nthreads else nthreads

            with ThreadPoolExecutor(max_workers=nthreads) as executor:
                res = executor.map(
                    cov.getQuadtreeNoise,
                    [rstate for _ in range(nbootstraps)])

                for ibs, bs in enumerate(res):
                    bootstraps[ibs, :] = bs
        except ImportError:
            for ibs in range(nbootstraps):
                if not (ibs+1) % 5:
                    logger.info('Calculating noise realisation %d/%d.'
                                % (ibs+1, nbootstraps))
                bootstraps[ibs, :] = cov.getQuadtreeNoise(rstate=rstate)

        self.set_bootstrap_residuals(bootstraps)

    @classmethod
    def get_plot_classes(cls):
        from . import plot
        plots = super(SatelliteMisfitTarget, cls).get_plot_classes()
        plots.extend(plot.get_plot_classes())
        return plots


__all__ = '''
    SatelliteTargetGroup
    SatelliteMisfitConfig
    SatelliteMisfitTarget
    SatelliteMisfitResult
'''.split()
