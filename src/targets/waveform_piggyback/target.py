import logging

import numpy as num

from pyrocko.guts import Float, Object, Int
from pyrocko import gf

from ..base import (
    MisfitTarget, TargetGroup, MisfitResult)

guts_prefix = 'grond'
logger = logging.getLogger('grond.targets.waveform_piggyback.target')


class WaveformPiggybackTargetGroup(TargetGroup):
    associated_path = gf.StringID.T()
    norm_exponent = Int.T(default=2)

    def get_targets(self, ds, event, default_path):
        logger.debug('Selecting waveform piggyback targets...')

        target = WaveformPiggybackTarget(
            path=self.path,
            norm_exponent=self.norm_exponent,
            associated_path=self.associated_path)

        return [target]


class WaveformPiggybackSubresult(Object):
    piggy_id = Int.T()
    amplitude_obs = Float.T()
    amplitude_syn = Float.T()


class WaveformPiggybackSubtarget(Object):

    piggy_id = Int.T()

    def evaluate(
            self, tr_proc_obs, trspec_proc_obs, tr_proc_syn, trspec_proc_syn):

        a_obs = num.sqrt(num.mean(num.abs(tr_proc_obs.ydata**2)))
        a_syn = num.sqrt(num.mean(num.abs(tr_proc_syn.ydata**2)))

        res = WaveformPiggybackSubresult(
            piggy_id=self.piggy_id,
            amplitude_obs=float(a_obs),
            amplitude_syn=float(a_syn))

        return res


class WaveformPiggybackMisfitResult(MisfitResult):
    pass


class WaveformPiggybackTarget(MisfitTarget):
    associated_path = gf.StringID.T()
    norm_exponent = Int.T(default=2)

    _next_piggy_id = 0

    def __init__(self, **kwargs):
        MisfitTarget.__init__(self, **kwargs)
        self.piggy_ids = set()

    @classmethod
    def new_piggy_id(cls):
        piggy_id = cls._next_piggy_id
        cls._next_piggy_id += 1
        return piggy_id

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

                piggy_id = self.new_piggy_id()

                self.piggy_ids.add(piggy_id)

                target.add_piggyback_subtarget(
                    WaveformPiggybackSubtarget(
                        piggy_id=piggy_id))

                relevant.append(target)

        return relevant

    def finalize_modelling(
            self, engine, source, modelling_targets, modelling_results):

        from ..waveform.target import WaveformMisfitResult
        amps_obs = []
        amps_syn = []
        for mtarget, mresult in zip(modelling_targets, modelling_results):
            if isinstance(mresult, WaveformMisfitResult):
                for sr in mresult.piggyback_subresults:
                    try:
                        self.piggy_ids.remove(sr.piggy_id)
                        amps_obs.append(sr.amplitude_obs)
                        amps_syn.append(sr.amplitude_syn)
                    except KeyError:
                        pass

        amps_obs = num.array(amps_obs)
        amps_syn = num.array(amps_syn)

        amp_obs = num.median(amps_obs[num.isfinite(amps_obs)])
        amp_syn = num.median(amps_syn[num.isfinite(amps_syn)])
        m = num.abs(num.log(amp_obs / amp_syn))**self.norm_exponent

        result = WaveformPiggybackMisfitResult(
            misfits=num.array([[m, 1.]], dtype=num.float))

        return result


__all__ = '''
    WaveformPiggybackTargetGroup
    WaveformPiggybackTarget
    WaveformPiggybackMisfitResult
    WaveformPiggybackSubtarget
    WaveformPiggybackSubresult
'''.split()
