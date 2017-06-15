
import copy
import numpy as num
from pyrocko.guts import Object, Int
from ..targets import TargetAnalysisResult
from ..meta import Forbidden


guts_prefix = 'grond'


class Analyser(object):

    def __init__(self, niter):
        self.niter = niter

    def analyse(self, problem, notifier):
        if self.niter == 0:
            return

        wtargets = []
        if not problem.has_waveforms:
            return

        for target in problem.waveform_targets:
            wtarget = copy.copy(target)
            wtarget.flip_norm = True
            wtarget.weight = 1.0
            wtargets.append(wtarget)

        wproblem = problem.copy()
        wproblem.targets = wtargets

        xbounds = num.array(wproblem.get_parameter_bounds(), dtype=num.float)
        npar = xbounds.shape[0]

        mss = num.zeros((self.niter, wproblem.ntargets))
        rstate = num.random.RandomState(123)

        notifier.emit('progress_start', 'analysing problem', self.niter)

        isbad_mask = None
        for iiter in xrange(self.niter):
            while True:
                x = []
                for ipar in xrange(npar):
                    v = rstate.uniform(xbounds[ipar, 0], xbounds[ipar, 1])
                    x.append(v)

                try:
                    x = wproblem.preconstrain(x)
                    break

                except Forbidden:
                    pass

            if isbad_mask is not None and num.any(isbad_mask):
                isok_mask = num.logical_not(isbad_mask)
            else:
                isok_mask = None
            _, ms = wproblem.evaluate(x, mask=isok_mask)
            mss[iiter, :] = ms

            isbad_mask = num.isnan(ms)
            notifier.emit('progress_update', 'analysing problem', iiter)

        notifier.emit('progress_finish', 'analysing problem')

        mean_ms = num.mean(mss, axis=0)
        weights = 1. / mean_ms
        groups, ngroups = wproblem.get_group_mask()

        for igroup in xrange(ngroups):
            weights[groups == igroup] /= (
                num.nansum(weights[groups == igroup]) /
                num.nansum(num.isfinite(weights[groups == igroup])))

        for weight, target in zip(weights, problem.waveform_targets):
            target.analysis_result = TargetAnalysisResult(
                balancing_weight=float(weight))


class AnalyserConfig(Object):
    niter = Int.T(default=1000)

    def get_analyser(self):
        return Analyser(niter=self.niter)


__all__ = '''
    Analyser
    AnalyserConfig
'''.split()
