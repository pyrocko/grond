import logging
import numpy as num
from scipy import linalg as splinalg

from pyrocko import gf
from pyrocko.guts import String, Dict, List, Int

from ..base import MisfitConfig, MisfitTarget, MisfitResult, TargetGroup
from grond.meta import has_get_plot_classes

guts_prefix = 'grond'
logger = logging.getLogger('grond.target').getChild('gnss_campaign')


class GNSSCampaignMisfitResult(MisfitResult):
    """Carries the observations for a target and corresponding synthetics. """
    statics_syn = Dict.T(optional=True,
                         help='Synthetic gnss surface displacements')
    statics_obs = Dict.T(optional=True,
                         help='Observed gnss surface displacements')


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

    def get_targets(self, ds, event, default_path='none'):
        logger.debug('Selecting GNSS targets...')
        targets = []

        for camp in ds.get_gnss_campaigns():
            if camp.name not in self.gnss_campaigns and\
               '*all' not in self.gnss_campaigns:
                continue

            if not isinstance(self.misfit_config,
                              GNSSCampaignMisfitConfig):
                raise AttributeError('misfit_config must be of type'
                                     ' GNSSCampaignMisfitConfig.')

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
                ncomponents=camp.ncomponents,
                tsnapshot=None,
                interpolation=self.interpolation,
                store_id=self.store_id,
                normalisation_family=self.normalisation_family,
                path=self.path or default_path,
                misfit_config=self.misfit_config)

            gnss_target.set_dataset(ds)
            targets.append(gnss_target)

        return targets


@has_get_plot_classes
class GNSSCampaignMisfitTarget(gf.GNSSCampaignTarget, MisfitTarget):
    """Handles and carries out operations related to the objective functions.

    The objective function is here the weighted misfit between observed
    and predicted surface displacements.
    """
    campaign_name = String.T()
    ncomponents = Int.T(optional=True)
    misfit_config = GNSSCampaignMisfitConfig.T()

    can_bootstrap_weights = True
    can_bootstrap_residuals = True

    plot_misfits_cumulative = False

    def __init__(self, **kwargs):
        gf.GNSSCampaignTarget.__init__(self, **kwargs)
        MisfitTarget.__init__(self, **kwargs)
        self._obs_data = None
        self._sigma = None
        self._weights = None
        self._correlated_weights = None
        self._station_component_mask = None

    @property
    def id(self):
        return self.campaign_name

    def string_id(self):
        return self.campaign_name

    def misfits_string_ids(self):
        id_list = []
        for station in self.campaign.stations:
            for name in ('north', 'east', 'up'):
                component = station.__getattribute__(name)
                if not component:
                    continue
                id_list.append('%s.%s.%s' %
                               (self.path, station.code, name[0].upper()))
        return id_list

    @property
    def station_names(self):
        return ['%s' % (station.code)
                for station in self.campaign.stations]

    @property
    def nmisfits(self):
        if self.ncomponents is None:
            try:
                self.campaign.ncomponents
            except ValueError:
                raise ValueError('Set the dataset!')
        return self.ncomponents

    @property
    def nstations(self):
        return self.lats.size

    def set_dataset(self, ds):
        MisfitTarget.set_dataset(self, ds)

    def get_correlated_weights(self, nthreads=0):
        if self._correlated_weights is None:
            self._correlated_weights = splinalg.sqrtm(self.weights)
        return self._correlated_weights

    @property
    def campaign(self):
        return self._ds.get_gnss_campaign(self.campaign_name)

    @property
    def obs_data(self):
        if self._obs_data is None:
            self._obs_data = num.concatenate(
                [s.get_displacement_data() for s in self.campaign.stations])
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
        """Weights are the square-rooted, inverted the data error variance-covariance.

        The single component variances, and if provided the component
        covariances, are used to build a data variance matrix or
        variance-covariance matrix.

        This matrix has the size for all possible NEU components,
        but throws zeros for not given components, also recorded in
        the _station_component_mask.
        """
        if self._weights is None:
            covar = self.campaign.get_covariance_matrix()

            if not num.any(covar.diagonal()):
                logger.warning('GNSS Stations have an empty covariance matrix.'
                               ' Weights will be all equal.')
                num.fill_diagonal(covar, 1.)

            self._weights = num.asmatrix(covar).I

        return self._weights

    @property
    def station_component_mask(self):
        if self._station_component_mask is None:
            self._station_component_mask = self.campaign.get_component_mask()
        return self._station_component_mask

    @property
    def station_weights(self):
        weights = num.diag(self.weights)
        return num.mean([weights[0::3], weights[1::3], weights[2::3]], axis=0)

    def component_weights(self):
        ws = num.sum(self.weights, axis=0)
        return ws

    def post_process(self, engine, source, statics):
        """Applies the objective function.

        As a result the weighted misfits are given and the observed and
        synthetic data.
        """
        obs = self.obs_data

        # All data is ordered in vectors as
        # S1_n, S1_e, S1_u, ..., Sn_n, Sn_e, Sn_u. Hence (.ravel(order='F'))
        syn = num.array([
              statics['displacement.n'],
              statics['displacement.e'],
              -statics['displacement.d']])\
            .ravel(order='F')

        syn = syn[self.station_component_mask]
        res = obs - syn

        misfit_value = res
        misfit_norm = obs

        mf = num.vstack((misfit_value, misfit_norm)).T

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

    def init_bootstrap_residuals(self, nbootstraps, rstate=None, nthreads=0):
        logger.info('GNSS campaign %s, bootstrapping residuals'
                    ' from measurement uncertainties ...'
                    % self.campaign.name)
        if rstate is None:
            rstate = num.random.RandomState()

        campaign = self.campaign
        bootstraps = num.empty((nbootstraps, self.ncomponents))

        sigmas = num.diag(num.sqrt(campaign.get_covariance_matrix()))
        if not num.all(sigmas):
            logger.warning('Bootstrapping GNSS stations is meaningless,'
                           ' all station\'s sigma are 0.0!')

        for ibs in range(nbootstraps):
            syn_noise = rstate.normal(scale=sigmas.ravel())
            bootstraps[ibs, :] = syn_noise

        self.set_bootstrap_residuals(bootstraps)

    @classmethod
    def get_plot_classes(cls):
        from . import plot
        plots = super(GNSSCampaignMisfitTarget, cls).get_plot_classes()
        plots.extend(plot.get_plot_classes())
        return plots
