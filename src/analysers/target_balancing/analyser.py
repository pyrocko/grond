import copy
import time
import logging
import numpy as num
from pyrocko.guts import Int, Float, Bool
from pyrocko import gf

from grond.meta import Forbidden, GrondError

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

    def __init__(self, niter, use_reference_magnitude, cutoff):
        Analyser.__init__(self)
        self.niter = niter
        self.use_reference_magnitude = use_reference_magnitude
        self.cutoff = cutoff

    def log_progress(self, problem, iiter, niter):
        t = time.time()
        if self._tlog_last < t - 10. \
                or iiter == 0 \
                or iiter == niter - 1:

            logger.info(
                'Target balancing for "%s" at %i/%i.' % (
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

        misfits = num.zeros((self.niter, wproblem.ntargets, 2))
        rstate = num.random.RandomState(123)

        isbad_mask = None

        self._tlog_last = 0
        for iiter in range(self.niter):
            self.log_progress(problem, iiter, self.niter)
            while True:
                if self.use_reference_magnitude:
                    try:
                        fixed_magnitude = wproblem.base_source.get_magnitude()
                    except gf.DerivedMagnitudeError:
                        raise GrondError(
                            'Cannot use use_reference_magnitude for this type '
                            'of source model.')
                else:
                    fixed_magnitude = None

                x = wproblem.random_uniform(
                    xbounds, rstate, fixed_magnitude=fixed_magnitude)

                try:
                    x = wproblem.preconstrain(x)
                    break

                except Forbidden:
                    pass

            if isbad_mask is not None and num.any(isbad_mask):
                isok_mask = num.logical_not(isbad_mask)
            else:
                isok_mask = None
            misfits[iiter, :, :] = wproblem.misfits(x, mask=isok_mask)

            isbad_mask = num.isnan(misfits[iiter, :, 1])

        mean_ms = num.mean(misfits[:, :, 0], axis=0)

        mean_ps = num.mean(misfits[:, :, 1], axis=0)

        weights = 1. / mean_ps
        families, nfamilies = wproblem.get_family_mask()

        for ifamily in range(nfamilies):
            weights[families == ifamily] /= (
                num.nansum(weights[families == ifamily]) /
                num.nansum(num.isfinite(weights[families == ifamily])))

        if self.cutoff is not None:
            weights[mean_ms / mean_ps > self.cutoff] = 0.0

        for weight, target in zip(weights, problem.waveform_targets):
            target.analyser_results['target_balancing'] = \
                TargetBalancingAnalyserResult(weight=float(weight))

        for itarget, target in enumerate(problem.waveform_targets):
            logger.info((
                'Balancing analysis for target "%s":\n'
                '  m/p: %g\n'
                '  weight: %g\n'
                ) % (
                    target.string_id(),
                    mean_ms[itarget] / mean_ps[itarget],
                    weights[itarget]))


class TargetBalancingAnalyserResult(AnalyserResult):
    weight = Float.T()


class TargetBalancingAnalyserConfig(AnalyserConfig):
    """Configuration parameters of the target balancing."""
    niterations = Int.T(
        default=1000,
        help='Number of random forward models for mean phase amplitude '
             'estimation')

    use_reference_magnitude = Bool.T(
        default=False,
        help='Fix magnitude of random sources to the magnitude of the '
             'reference event.')

    cutoff = Float.T(
        optional=True,
        help='Remove targets where ratio m/p > cutoff, where m is the misfit '
             'between synthetics and observations and p is the misfit between '
             'synthetics and zero-traces. Magnitude should be fixed to use '
             'this.')

    def get_analyser(self):
        return TargetBalancingAnalyser(
            niter=self.niterations,
            use_reference_magnitude=self.use_reference_magnitude,
            cutoff=self.cutoff)


__all__ = '''
    TargetBalancingAnalyser
    TargetBalancingAnalyserConfig
'''.split()
