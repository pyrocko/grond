import logging
import numpy as num

from pyrocko import gf
from pyrocko.guts import String, Dict, List

from ..base import MisfitConfig, MisfitTarget, MisfitResult, TargetGroup

guts_prefix = 'grond'
logger = logging.getLogger('grond.target').getChild('gnss_campaign')


class GNSSCampaignMisfitResult(MisfitResult):
    statics_syn = Dict.T(optional=True)
    statics_obs = Dict.T(optional=True)


class GNSSCampaignMisfitConfig(MisfitConfig):
    pass


class GNSSCampaignTargetGroup(TargetGroup):
    gnss_campaigns = List.T(optional=True)
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

    @classmethod
    def get_plot_classes(cls):
        from . import plot
        plots = super(GNSSCampaignMisfitTarget, cls).get_plot_classes()
        plots.extend(plot.get_plot_classes())
        return plots


class GNSSCampaignMisfitTarget(gf.GNSSCampaignTarget, MisfitTarget):
    campaign_name = String.T()
    misfit_config = GNSSCampaignMisfitConfig.T()

    def __init__(self, **kwargs):
        gf.GNSSCampaignTarget.__init__(self, **kwargs)
        MisfitTarget.__init__(self, **kwargs)
        self._obs_data = None
        self._sigma = None
        self._weights = None

    @property
    def id(self):
        return self.campaign_name

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
                [s.up.shift for s in self.campaign.stations]])
        return self._obs_data

    @property
    def obs_sigma(self):
        if self._sigma is None:
            self._sigma = num.array([
                [s.north.sigma for s in self.campaign.stations],
                [s.east.sigma for s in self.campaign.stations],
                [s.up.sigma for s in self.campaign.stations]])
        return self._obs_data

    @property
    def weights(self):
        if self._weights is None:
            self._weights = 1./self.obs_sigma
            self._weights[self._weights == num.inf] = 1.
        return self._weights

    def post_process(self, engine, source, statics):
        obs = self.obs_data
        weights = self.weights
        syn = num.array([statics['displacement.n'],
                         statics['displacement.e'],
                         -statics['displacement.d']])

        res = obs - syn

        misfit_value = num.sqrt(
            num.sum((res * weights)**2))
        misfit_norm = num.sqrt(
            num.sum((obs * weights)**2))
        result = GNSSCampaignMisfitResult(
            misfits=num.array([misfit_value, misfit_norm],
                              dtype=num.float))

        if self._result_mode == 'full':
            result.statics_syn = statics
            result.statics_obs = obs

        return result

    def get_combined_weight(self, apply_balancing_weights=False,
                            apply_station_noise_weights=False):
        return num.array([self.manual_weight])

    def prepare_modelling(self, engine, source):
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
