import logging

import numpy as num

from pyrocko.guts import Float, Tuple, String, Bool
from pyrocko import gf

from ..base import (
    MisfitTarget, TargetGroup, MisfitResult)
from . import measure as fm
from grond import dataset
from grond.meta import has_get_plot_classes

from ..waveform.target import StoreIDSelector

guts_prefix = 'grond'
logger = logging.getLogger('grond.targets.waveform_phase_ratio.target')


def log_exclude(target, reason):
    logger.debug('Excluding potential target %s: %s' % (
        target.string_id(), reason))


class PhaseRatioTargetGroup(TargetGroup):

    '''
    Generate targets to compare ratios or log ratios of two seismic phases.

      misfit = | a_obs / (a_obs + b_obs)  - a_syn / (a_syn + b_syn) |

    or

      misfit = | log(a_obs / (a_obs + b_obs) + waterlevel)  -
                 log(a_syn / (a_syn + b_syn) + waterlevel) |
    '''

    distance_min = Float.T(optional=True)
    distance_max = Float.T(optional=True)
    distance_3d_min = Float.T(optional=True)
    distance_3d_max = Float.T(optional=True)
    depth_min = Float.T(optional=True)
    depth_max = Float.T(optional=True)
    measure_a = fm.FeatureMeasure.T()
    measure_b = fm.FeatureMeasure.T()
    interpolation = gf.InterpolationMethod.T()
    store_id = gf.StringID.T(optional=True)
    store_id_selector = StoreIDSelector.T(
            optional=True,
            help='select GF store based on event-station geometry.')

    fit_log_ratio = Bool.T(
        default=True,
        help='If true, compare synthetic and observed log ratios')

    fit_log_ratio_waterlevel = Float.T(
        default=0.01,
        help='Waterlevel added to both ratios when comparing on logarithmic '
             'scale, to avoid log(0)')

    def get_targets(self, ds, event, default_path='none'):
        logger.debug('Selecting phase ratio targets...')
        origin = event
        targets = []

        for st in ds.get_stations():
            blacklisted = False
            for measure in [self.measure_a, self.measure_b]:
                for cha in measure.channels:
                    if ds.is_blacklisted((st.nsl() + (cha,))):
                        blacklisted = True

            if self.store_id_selector:
                store_id = self.store_id_selector.get_store_id(
                    event, st, cha)
            else:
                store_id = self.store_id

            target = PhaseRatioTarget(
                codes=st.nsl(),
                lat=st.lat,
                lon=st.lon,
                north_shift=st.north_shift,
                east_shift=st.east_shift,
                depth=st.depth,
                interpolation=self.interpolation,
                store_id=store_id,
                measure_a=self.measure_a,
                measure_b=self.measure_b,
                manual_weight=self.weight,
                normalisation_family=self.normalisation_family,
                path=self.path or default_path,
                backazimuth=0.0,
                fit_log_ratio=self.fit_log_ratio,
                fit_log_ratio_waterlevel=self.fit_log_ratio_waterlevel)

            if blacklisted:
                log_exclude(target, 'blacklisted')
                continue

            if self.distance_min is not None and \
               target.distance_to(origin) < self.distance_min:
                log_exclude(target, 'distance < distance_min')
                continue

            if self.distance_max is not None and \
               target.distance_to(origin) > self.distance_max:
                log_exclude(target, 'distance > distance_max')
                continue

            if self.distance_3d_min is not None and \
               target.distance_3d_to(origin) < self.distance_3d_min:
                log_exclude(target, 'distance_3d < distance_3d_min')
                continue

            if self.distance_3d_max is not None and \
               target.distance_3d_to(origin) > self.distance_3d_max:
                log_exclude(target, 'distance_3d > distance_3d_max')
                continue

            if self.depth_min is not None and \
               target.depth < self.depth_min:
                log_exclude(target, 'depth < depth_min')
                continue

            if self.depth_max is not None and \
               target.depth > self.depth_max:
                log_exclude(target, 'depth > depth_max')
                continue

            bazi, _ = target.azibazi_to(origin)
            target.backazimuth = bazi
            target.set_dataset(ds)
            targets.append(target)

        return targets


class PhaseRatioResult(MisfitResult):
    a_obs = Float.T(optional=True)
    b_obs = Float.T(optional=True)
    a_syn = Float.T(optional=True)
    b_syn = Float.T(optional=True)


@has_get_plot_classes
class PhaseRatioTarget(gf.Location, MisfitTarget):

    '''
    Target to compare ratios or log ratios of two seismic phases.

      misfit = | a_obs / (a_obs + b_obs)  - a_syn / (a_syn + b_syn) |

    '''

    codes = Tuple.T(
        3, String.T(),
        help='network, station, location codes.')

    store_id = gf.StringID.T(
        help='ID of Green\'s function store to use for the computation.')

    backazimuth = Float.T(optional=True)

    interpolation = gf.InterpolationMethod.T()

    measure_a = fm.FeatureMeasure.T()
    measure_b = fm.FeatureMeasure.T()

    fit_log_ratio = Bool.T(
        default=True,
        help='if true, compare synthetic and observed log ratios')

    fit_log_ratio_waterlevel = Float.T(
        default=0.01,
        help='Waterlevel added to both ratios when comparing on logarithmic '
             'scale, to avoid log(0)')

    can_bootstrap_weights = True

    def __init__(self, **kwargs):
        gf.Location.__init__(self, **kwargs)
        MisfitTarget.__init__(self, **kwargs)

    @classmethod
    def get_plot_classes(cls):
        from . import plot
        plots = super(PhaseRatioTarget, cls).get_plot_classes()
        plots.extend(plot.get_plot_classes())
        return plots

    def string_id(self):
        return '.'.join(x for x in (self.path,) + self.codes)

    def get_plain_targets(self, engine, source):
        return self.prepare_modelling(engine, source, None)

    def prepare_modelling(self, engine, source, targets):
        modelling_targets = []
        for measure in [self.measure_a, self.measure_b]:
            modelling_targets.extend(measure.get_modelling_targets(
                self.codes, self.lat, self.lon, self.depth, self.store_id,
                self.backazimuth))

        return modelling_targets

    def finalize_modelling(
            self, engine, source, modelling_targets, modelling_results):

        ds = self.get_dataset()

        try:
            imt = 0
            amps = []
            for measure in [self.measure_a, self.measure_b]:
                nmt_this = measure.get_nmodelling_targets()
                amp_obs, _ = measure.evaluate(
                    engine, source,
                    modelling_targets[imt:imt+nmt_this],
                    dataset=ds)

                amp_syn, _ = measure.evaluate(
                    engine, source,
                    modelling_targets[imt:imt+nmt_this],
                    trs=[r.trace.pyrocko_trace()
                         for r
                         in modelling_results[imt:imt+nmt_this]])

                amps.append((amp_obs, amp_syn))

                imt += nmt_this

            (a_obs, a_syn), (b_obs, b_syn) = amps

            eps = self.fit_log_ratio_waterlevel
            if self.fit_log_ratio:
                res_a = num.log(a_obs / (a_obs + b_obs) + eps) \
                    - num.log(a_syn / (a_syn + b_syn) + eps)
            else:
                res_a = a_obs / (a_obs + b_obs) - a_syn / (a_syn + b_syn)

            misfit = num.abs(res_a)
            norm = 1.0

            result = PhaseRatioResult(
                misfits=num.array([[misfit, norm]], dtype=num.float),
                a_obs=a_obs,
                b_obs=b_obs,
                a_syn=a_syn,
                b_syn=b_syn)

            return result

        except dataset.NotFound as e:
            logger.debug(str(e))
            return gf.SeismosizerError('No waveform data: %s' % str(e))


__all__ = '''
    PhaseRatioTargetGroup
    PhaseRatioTarget
    PhaseRatioResult
'''.split()
