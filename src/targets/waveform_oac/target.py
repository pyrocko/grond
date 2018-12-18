import logging

import numpy as num

from pyrocko.guts import Float, Int
from pyrocko import gf

from ..base import (
    MisfitTarget, TargetGroup, MisfitResult)

from ..waveform.target import WaveformPiggybackSubtarget, \
    WaveformPiggybackSubresult

guts_prefix = 'grond'
logger = logging.getLogger('grond.targets.waveform_piggyback.target')


class WOACSubtarget(WaveformPiggybackSubtarget):

    def evaluate(
            self, tr_proc_obs, trspec_proc_obs, tr_proc_syn, trspec_proc_syn):

        a_obs = num.sqrt(num.mean(num.abs(tr_proc_obs.ydata**2)))
        a_syn = num.sqrt(num.mean(num.abs(tr_proc_syn.ydata**2)))

        res = WOACSubresult(
            piggy_id=self.piggy_id,
            amplitude_obs=float(a_obs),
            amplitude_syn=float(a_syn))

        return res


class WOACSubresult(WaveformPiggybackSubresult):
    amplitude_obs = Float.T()
    amplitude_syn = Float.T()


class WaveformOverallAmplitudeConstraint(TargetGroup):
    associated_path = gf.StringID.T()
    norm_exponent = Int.T(default=2)

    def get_targets(self, ds, event, default_path='none'):
        logger.debug('Selecting waveform piggyback targets...')

        target = WOACTarget(
            path=self.path,
            norm_exponent=self.norm_exponent,
            associated_path=self.associated_path)

        return [target]


class WOACTarget(MisfitTarget):
    associated_path = gf.StringID.T()
    norm_exponent = Int.T(default=2)

    can_bootstrap_weights = True

    def __init__(self, **kwargs):
        MisfitTarget.__init__(self, **kwargs)
        self.piggy_ids = set()

    def string_id(self):
        return self.path

    def get_plain_targets(self, engine, source):
        return []

    def distance_to(self, other):
        return num.array([], dtype=float)

    def prepare_modelling(self, engine, source, targets):
        from ..waveform.target import WaveformMisfitTarget
        relevant = []
        for target in targets:
            if isinstance(target, WaveformMisfitTarget) \
                    and target.path == self.associated_path:

                st = WOACSubtarget()
                self.piggy_ids.add(st.piggy_id)
                target.add_piggyback_subtarget(st)
                relevant.append(target)

        return relevant

    def finalize_modelling(
            self, engine, source, modelling_targets, modelling_results):

        from ..waveform.target import WaveformMisfitResult
        amps = []
        for mtarget, mresult in zip(modelling_targets, modelling_results):
            if isinstance(mresult, WaveformMisfitResult):
                for sr in list(mresult.piggyback_subresults):
                    try:
                        self.piggy_ids.remove(sr.piggy_id)
                        amps.append((sr.amplitude_obs, sr.amplitude_syn))
                        mresult.piggyback_subresults.remove(sr)
                    except KeyError:
                        logger.error(
                            'Found inconsistency while gathering piggyback '
                            'results.')
                        pass

        amps = num.array(amps, dtype=num.float)
        mask = num.all(num.isfinite(amps), axis=1)

        amp_obs = num.median(amps[mask, 0])
        amp_syn = num.median(amps[mask, 1])
        m = num.abs(num.log(amp_obs / amp_syn))**self.norm_exponent

        result = WOACMisfitResult(
            misfits=num.array([[m, 1.]], dtype=num.float))

        return result


class WOACMisfitResult(MisfitResult):
    pass


__all__ = '''
    WOACSubtarget
    WOACSubresult
    WaveformOverallAmplitudeConstraint
    WOACTarget
    WOACMisfitResult
'''.split()
