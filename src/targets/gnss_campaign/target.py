import logging
import numpy as num

from pyrocko import gf
from pyrocko.guts import String, Dict, List

from ..base import MisfitConfig, MisfitTarget, MisfitResult, TargetGroup

guts_prefix = 'grond'
logger = logging.getLogger('grond.target').getChild('gnss_campaign')


class GNSSCampaignMisfitResult(MisfitResult):
    """Carries the observations for a target and corresponding synthetics. """
    statics_syn = Dict.T(optional=True,
                         help='synthetic gnss surface displacements')
    statics_obs = Dict.T(optional=True,
                         help='observed gnss surface displacements')


class GNSSCampaignMisfitConfig(MisfitConfig):

    pass


class GNSSCampaignTargetGroup(TargetGroup):
    """Handles static displacements from campaign GNSS observations, e.g GPS.

    Station information, displacements and measurement errors are provided in
    a `yaml`-file (please find an example in the documentation). The
    measurement errors may consider correlations between components of a
    station, but correlations between single stations is not considered.
    """
    gnss_campaigns = List.T(
        optional=True,
        help='List of individual campaign names'
             ' (`name` in `gnss.yaml` files).')
    misfit_config = GNSSCampaignMisfitConfig.T()

    def get_targets(self, ds, event, default_path):
        logger.debug('Selecting GNSS targets...')
        targets = []

        for camp in ds.get_gnss_campaigns():
            if camp.name not in self.gnss_campaigns and\
               '*all' not in self.gnss_campaigns:
                continue

            if not isinstance(self.misfit_config,
                              GNSSCampaignMisfitConfig):
                raise AttributeError('misfit_config must be of type'
                                     ' GNSSCampaignMisfitConfig')

            lats = num.array([s.lat for s in camp.stations])
            lons = num.array([s.lon for s in camp.stations])

            north_shifts = num.array([s.north_shift for s in camp.stations])
            east_shifts = num.array([s.east_shift for s in camp.stations])

            gnss_target = GNSSCampaignMisfitTarget(
                quantity='displacement',
                campaign_name=camp.name,
                lats=lats,
                lons=lons,
                east_shifts=east_shifts,
                north_shifts=north_shifts,
                tsnapshot=None,
                interpolation=self.interpolation,
                store_id=self.store_id,
                normalisation_family=self.normalisation_family,
                path=self.path or default_path,
                misfit_config=self.misfit_config)

            gnss_target.set_dataset(ds)
            targets.append(gnss_target)

        return targets


class GNSSCampaignMisfitTarget(gf.GNSSCampaignTarget, MisfitTarget):
    """Handles and carries out operations related to the objective functions.

    The objective function is here the weighted misfit between observed
    and predicted surface displacements.
    """
    campaign_name = String.T()
    misfit_config = GNSSCampaignMisfitConfig.T()

    can_bootstrap_weights = True

    def __init__(self, **kwargs):
        gf.GNSSCampaignTarget.__init__(self, **kwargs)
        MisfitTarget.__init__(self, **kwargs)
        self._obs_data = None
        self._sigma = None
        self._weights = None

    @property
    def id(self):
        return self.campaign_name

    def string_id(self):
        return self.campaign_name

    @property
    def nmisfits(self):
        return self.lats.size

    def set_dataset(self, ds):
        MisfitTarget.set_dataset(self, ds)

    @property
    def campaign(self):
        return self._ds.get_gnss_campaign(self.campaign_name)

    @property
    def obs_data(self):
        if self._obs_data is None:
            self._obs_data = num.array([
                [s.north.shift for s in self.campaign.stations],
                [s.east.shift for s in self.campaign.stations],
                [s.up.shift for s in self.campaign.stations]])\
              .ravel(order='F')
        return self._obs_data

    @property
    def obs_sigma(self):
        if self._sigma is None:
            self._sigma = num.array([
                [s.north.sigma for s in self.campaign.stations],
                [s.east.sigma for s in self.campaign.stations],
                [s.up.sigma for s in self.campaign.stations]])\
              .ravel(order='F')
        return self._sigma

    @property
    def weights(self):
        """Weights are the inverse of the data error variance-covariance.

        The single component variances, and if provided the component
        covariances, are used to build a data variance matrix or
        variance-covariance matrix. Correlations between stations are
        not implemented.
        """
        if self._weights is None:
            covar = self.campaign.get_covariance_matrix()

            if not num.all(covar.diagonal()):
                logger.warning('GNSS Stations have an empty covariance matrix.'
                               ' Weights will be all equal.')
                num.fill_diagonal(covar, 1.)
            self._weights = num.asmatrix(covar).I
        return self._weights

    def post_process(self, engine, source, statics):
        """Applies the objective function.

        As a result the weighted misfits are given and the observed and
        synthetic data.
        """
        obs = self.obs_data
        weights = self.weights
        nstations = self.campaign.nstations

        # All data is ordered in vectors as
        # S1_n, S1_e, S1_u, ..., Sn_n, Sn_e, Sn_u. Hence (.ravel(order='F'))
        syn = num.array([statics['displacement.n'],
                         statics['displacement.e'],
                         -statics['displacement.d']])\
            .ravel(order='F')

        res = obs - syn

        misfit_value = res * weights
        misfit_norm = obs * weights

        misfit_value = num.sum(
            misfit_value.reshape((nstations, 3)), axis=1)
        misfit_norm = num.sum(
            misfit_norm.reshape((nstations, 3)), axis=1)

        mf = num.hstack((misfit_value, misfit_norm))
        result = GNSSCampaignMisfitResult(
            misfits=mf)

        if self._result_mode == 'full':
            result.statics_syn = statics
            result.statics_obs = obs

        return result

    def get_combined_weight(self):
        """A given manual weight in the configuration is applied."""
        if self._combined_weight is None:
            self._combined_weight = num.full(self.nmisfits, self.manual_weight)

        return self._combined_weight

    def prepare_modelling(self, engine, source, targets):
        return [self]

    def finalize_modelling(
            self, engine, source, modelling_targets, modelling_results):

        return modelling_results[0]

    @classmethod
    def get_plot_classes(cls):
        from . import plot
        plots = super(GNSSCampaignMisfitTarget, cls).get_plot_classes()
        plots.extend(plot.get_plot_classes())
        return plots
