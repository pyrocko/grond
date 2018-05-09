import copy
import time
import logging
import numpy as num
from pyrocko.guts import Int, Float

from grond.meta import Forbidden

from ..base import Analyser, AnalyserConfig, AnalyserResult

logger = logging.getLogger('grond.analysers.target_balancer')


guts_prefix = 'grond'


class TargetBalancingAnalyser(Analyser):
    """ Estimating target weights that balance the signal amplitudes.

     Signal amplitudes depend on the source-receiver distance, on the
     phase type and the taper used. Large signals have in general
     a higher contribution to the misfit than smaller signals,
     without carrying more information. With this function, weights are
     estimated that shall balance the phase contributions.

     The weight estimation is based on synthetic waveforms that stem from
     a given number of random forward models. The inverse of the mean
     synthetic signal amplitudes gives the balancing weight. This is
     described as adaptive station weighting in Heimann (2011).
     """

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
            target.analyser_results['target_balancing'] = \
                TargetBalancingAnalyserResult(weight=float(weight))


class TargetBalancingAnalyserResult(AnalyserResult):
    weight = Float.T()


class TargetBalancingAnalyserConfig(AnalyserConfig):
    """Configuration parameters of the target balancing."""
    niterations = Int.T(default=1000,
                        help='Number of random forward models for mean \
                             phase amplitude estimation')

    def get_analyser(self):
        return TargetBalancingAnalyser(niter=self.niterations)


__all__ = '''
    TargetBalancingAnalyser
    TargetBalancingAnalyserConfig
'''.split()
