import copy
import time
import logging
import numpy as num
from pyrocko.guts import Int, Bool

from grond.targets import TargetAnalysisResult
from grond.meta import Forbidden

from ..base import Analyser, AnalyserConfig

logger = logging.getLogger('grond.analysers.target_balancer')


guts_prefix = 'grond'


class TargetBalancingAnalyser(Analyser):

    def __init__(self, niter):
        Analyser.__init__(self)
        self.niter = niter

    def log_progress(self, problem, iiter, niter):
        t = time.time()
        if self._tlog_last < t - 10. \
                or iiter == 0 \
                or iiter == niter - 1:

            logger.info(
                '%s at %i/%i' % (
                    problem.name,
                    iiter, niter))

            self._tlog_last = t

    def analyse(self, problem, ds):
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

        xbounds = wproblem.get_parameter_bounds()
        npar = xbounds.shape[0]

        mss = num.zeros((self.niter, wproblem.ntargets))
        rstate = num.random.RandomState(123)

        isbad_mask = None

        self._tlog_last = 0
        for iiter in range(self.niter):
            self.log_progress(problem, iiter, self.niter)
            while True:
                x = []
                for ipar in range(npar):
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
            ms = wproblem.misfits(x, mask=isok_mask)[:, 1]
            mss[iiter, :] = ms

            isbad_mask = num.isnan(ms)

        mean_ms = num.mean(mss, axis=0)
        weights = 1. / mean_ms
        families, nfamilies = wproblem.get_family_mask()

        for ifamily in range(nfamilies):
            weights[families == ifamily] /= (
                num.nansum(weights[families == ifamily]) /
                num.nansum(num.isfinite(weights[families == ifamily])))

        for weight, target in zip(weights, problem.waveform_targets):
            target.analysis_result = TargetAnalysisResult(
                balancing_weight=float(weight))


class TargetBalancingAnalyserConfig(AnalyserConfig):
    niterations = Int.T(default=1000)

    def get_analyser(self):
        return TargetBalancingAnalyser(niter=self.niterations)


__all__ = '''
    TargetBalancingAnalyser
    TargetBalancingAnalyserConfig
'''.split()
