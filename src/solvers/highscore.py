import math
import logging

import numpy as num

from collections import OrderedDict

from pyrocko.guts import StringChoice, Int, Float, Bool, Object
from pyrocko.guts_array import Array

from ..meta import GrondError, Forbidden
from .base import Optimizer, OptimizerConfig

guts_prefix = 'grond'

logger = logging.getLogger('grond.optimizers.highscore')


def excentricity_compensated_probabilities(xs, sbx, factor):
    inonflat = num.where(sbx != 0.0)[0]
    scale = num.zeros_like(sbx)
    scale[inonflat] = 1.0 / (sbx[inonflat] * (factor if factor != 0. else 1.0))
    distances_sqr_all = num.sum(
        ((xs[num.newaxis, :, :] - xs[:, num.newaxis, :]) *
         scale[num.newaxis, num.newaxis, :])**2, axis=2)
    probabilities = 1.0 / num.sum(distances_sqr_all < 1.0, axis=1)
    print num.sort(num.sum(distances_sqr_all < 1.0, axis=1))
    probabilities /= num.sum(probabilities)
    return probabilities


def excentricity_compensated_choice(xs, sbx, factor):
    probabilities = excentricity_compensated_probabilities(
        xs, sbx, factor)
    r = num.random.random()
    ichoice = num.searchsorted(num.cumsum(probabilities), r)
    ichoice = min(ichoice, xs.shape[0]-1)
    return ichoice


def local_std(xs):
    ssbx = num.sort(xs, axis=0)
    dssbx = num.diff(ssbx, axis=0)
    mdssbx = num.median(dssbx, axis=0)
    return mdssbx * dssbx.shape[0] / 2.6


class Phase(Object):
    niterations = Int.T()
    ntries_preconstrain_limit = Int.T(default=1000)

    def set_problem(self, problem):
        self._problem = problem

    def get_sample(
            self, iiter, accept_sum, xhist, chains_i, local_sxs, mxs, covs,
            nlinks):

        assert 0 <= iiter < self.niterations

        ntries_preconstrain = 0
        for ntries_preconstrain in xrange(self.ntries_preconstrain_limit):
            try:
                return self._problem.preconstrain(self.get_raw_sample(
                    iiter, accept_sum, xhist, chains_i, local_sxs, mxs, covs, 
                    nlinks))

            except Forbidden:
                pass

        raise GrondError(
            'could not find any suitable candidate sample within %i tries' % (
                self.ntries_preconstrain_limit))


class InjectionPhase(Phase):
    xs_inject = Array.T(dtype=num.float, shape=(None, None))

    def get_raw_sample(
            self, iiter, accept_sum, xhist, chains_i, local_sxs, mxs, covs,
            nlinks):
        return self.xs_inject[iiter, :]


class UniformPhase(Phase):

    def get_raw_sample(
            self, iiter, accept_sum, xhist, chains_i, local_sxs, mxs, covs,
            nlinks):
        xbounds = self._problem.get_parameter_bounds()
        return self._problem.random_uniform(xbounds)


class DirectedPhase(Phase):
    scatter_scale = Float.T(optional=True)
    scatter_scale_begin = Float.T(optional=True)
    scatter_scale_end = Float.T(optional=True)
    starting_point = StringChoice.T(
        choices=['excentricity_compensated', 'random', 'mean'])

    sampler_distribution = StringChoice.T(
        choices=['normal', 'multivariate_normal'])

    def get_scatter_scale_factor(self, iiter):
        s = self.scatter_scale
        sa = self.scatter_scale_begin
        sb = self.scatter_scale_end

        assert s is None or (sa is None and sb is None)

        if sa != sb:
            tb = float(self.niterations-1)
            tau = tb/(math.log(sa) - math.log(sb))
            t0 = math.log(sa) * tau
            t = float(iiter)
            return num.exp(-(t-t0) / tau)

        else:
            return s or 1.0

    def get_raw_sample(
            self, iiter, accept_sum, xhist, chains_i, local_sxs, mxs, covs,
            nlinks):

        factor = self.get_scatter_scale_factor()
        npar = self._problem.nparameters
        pnames = self._problem.parameter_names
        xbounds = self._problem.get_parameter_bounds()

        ichain_choice = num.argmin(accept_sum)

        if self.starting_point == 'excentricity_compensated':
            ichoice = excentricity_compensated_choice(
                xhist[chains_i[ichain_choice, :], :],
                local_sxs[ichain_choice], 2.)

            xchoice = xhist[
                chains_i[ichain_choice, ichoice], :]

        elif self.starting_point == 'random':
            ichoice = num.random.randint(0, nlinks)
            xchoice = xhist[
                chains_i[ichain_choice, ichoice], :]

        elif self.starting_point == 'mean':
            xchoice = mxs[ichain_choice]

        else:
            assert False, 'invalid starting_point choice: %s' % (
                self.starting_point)

        ntries_sample = 0
        if self.sampler_distribution == 'normal':
            x = num.zeros(npar, dtype=num.float)
            for ipar in xrange(npar):
                ntries = 0
                while True:
                    if local_sxs[ichain_choice][ipar] > 0.:
                        v = num.random.normal(
                            xchoice[ipar],
                            factor*local_sxs[
                                ichain_choice][ipar])
                    else:
                        v = xchoice[ipar]

                    if xbounds[ipar, 0] <= v and \
                            v <= xbounds[ipar, 1]:

                        break

                    if ntries > self.ntries_sample_limit:
                        raise GrondError(
                            'failed to produce a suitable '
                            'candidate sample from normal '
                            'distribution')

                    ntries += 1

                x[ipar] = v

        elif self.sampler_distribution == 'multivariate_normal':
            ok_mask_sum = num.zeros(npar, dtype=num.int)
            while True:
                ntries_sample += 1
                xcandi = num.random.multivariate_normal(
                    xchoice, factor**2 * covs[ichain_choice])

                ok_mask = num.logical_and(
                    xbounds[:, 0] <= xcandi, xcandi <= xbounds[:, 1])

                if num.all(ok_mask):
                    break

                ok_mask_sum += ok_mask

                if ntries_sample > self.ntries_sample_limit:
                    raise GrondError(
                        'failed to produce a suitable candidate '
                        'sample from multivariate normal '
                        'distribution, (%s)' %
                        ', '.join('%s:%i' % xx for xx in
                                  zip(pnames, ok_mask_sum)))

            x = xcandi

        return x, ichain_choice, ichoice, ntries_sample

class OptimizerRunData(object):

    nmodels_capacity_min = 1024

    def __init__(self, problem):
        self.problem = problem
        self.nmodels_capacity = self.nmodels_capacity_min

    @property
    def nmodels(self):
        if self.models is None:
            return 0
        else:
            return self.models.shape[0]

    @nmodels.setter
    def nmodels(self, nmodels_new):
        assert 0 <= nmodels_new <= self.nmodels
        self.models = self._models_buffer[:nmodels_new, :]
        self.misfits = self._misfits_buffer[:nmodels_new, :, :]

    @property
    def nmodels_capacity(self):
        if self._models_buffer is None:
            return 0
        else:
            return self._models_buffer.shape[0]

    @nmodels_capacity.setter
    def nmodels_capacity(self, nmodels_capacity_new):
        if self.nmodels_capacity != self.nmodels_capacity_new:
            models_buffer = num.zeros(
                (nmodels_capacity_new, self.problem.nparameters),
                dtype=num.float)

            misfits_buffer = num.zeros(
                (nmodels_capacity_new, self.problem.ntargets, 2),
                dtype=num.float)

            ncopy = min(self.nmodels, nmodels_capacity_new)

            models_buffer[:ncopy, :] = self._models_buffer[:ncopy, :]
            misfits_buffer[:ncopy, :, :] = self._misfits_buffer[:ncopy, :, :]
            self._models_buffer = models_buffer
            self._misfits_buffer = misfits_buffer

    def clear(self):
        self.nmodels = 0
        self.nmodels_capacity = self.nmodels_capacity_min

    def append(self, models, misfits):
        assert models.shape[0] == misfits.shape[0]

        nmodels_add = models.shape[0]

        nmodels = self.nmodels
        nmodels_new = nmodels + nmodels_add

        nmodels_capacity_want = max(
            nmodels_capacity_min, nextpow2(nmodels_new))

        if nmodels_capacity_want != self.nmodels_capacity:
            self.nmodels_capacity = nmodels_capacity_want

        self._models_buffer[nmodels:nmodels+nmodels_add, :] = models
        self._misfits_buffer[nmodels:nmodels+nmodels_add, :, :] = misfits

        nmodels = nmodels_new

        self.models = self._models_buffer[:nmodels, :]
        self.misfits = self._misfits_buffer[:nmodels, :, :]


class Chains(object):
    def __init__(self, problem, nchains, nlinks_cap):
        problem


class BadProblem(GrondError):
    pass


def solve(problem,
          phases,
          rundir=None,

          chain_length_factor=8.0,
          standard_deviation_estimator='median_density_single_chain',

          status=(),
          plot=None,
          notifier=None,
          state=None):

    xbounds = problem.get_parameter_bounds()
    npar = problem.nparameters

    nlinks_cap = int(round(chain_length_factor * npar + 1))
    chains_m = num.zeros((1 + problem.nchains, nlinks_cap), num.float)
    chains_i = num.zeros((1 + problem.nchains, nlinks_cap), num.int)
    nlinks = 0
    mbx = None

    niter = sum(phase.nitererations for phase in phases)

    iiter = 0
    sbx = None
    mxs = None
    covs = None
    local_sxs = None
    xhist = num.zeros((niter, npar))
    isbad_mask = None
    accept_sum = num.zeros(1 + problem.nchains, dtype=num.int)
    accept_hist = num.zeros(niter, dtype=num.int)
    pnames = problem.parameter_names


    for par in state.parameter_sets.values():
        par.fill(num.nan)

    state.niter = niter

    if plot:
        plot.start(problem)

    phase = phases.pop(0)
    iiter_phase_start = 0
    while iiter < niter:

        if iiter - iiter_phase_starte == phase.niterations:
            phase = phases.pop(0)
            iiter_phase_start = iiter

        iiter_phase = iiter - iiter_phase_start

        x, ichain_choice, ichoice, ntries_preconstrain, ntries_sample = \
            phase.get_sample(
                iiter_phase, accept_sum, xhist, chains_i, local_sxs, mxs, covs,
                nlinks):

        ibase = None
        if ichoice is not None and ichain_choice is not None:
            ibase = chains_i[ichain_choice, ichoice]

        if isbad_mask is not None and num.any(isbad_mask):
            isok_mask = num.logical_not(isbad_mask)
        else:
            isok_mask = None

        ms, ns = problem.evaluate(x, mask=isok_mask)

        isbad_mask_new = num.isnan(ms)
        if isbad_mask is not None and num.any(isbad_mask != isbad_mask_new):
            errmess = [
                'problem %s: inconsistency in data availability' %
                problem.name]

            for target, isbad_new, isbad in zip(
                    problem.targets, isbad_mask_new, isbad_mask):

                if isbad_new != isbad:
                    errmess.append('  %s, %s -> %s' % (
                        target.string_id(), isbad, isbad_new))

            raise BadProblem(errmess)

        isbad_mask = isbad_mask_new

        if num.all(isbad_mask):
            raise BadProblem(
                'problem %s: all target misfit values are NaN' % problem.name)

        gm = problem.global_misfit(ms, ns)
        bms = problem.bootstrap_misfit(ms, ns)

        chains_m[0, nlinks] = gm
        chains_m[1:, nlinks] = bms
        chains_i[:, nlinks] = iiter

        nlinks += 1

        for ichain in xrange(chains_m.shape[0]):
            isort = num.argsort(chains_m[ichain, :nlinks])
            chains_m[ichain, :nlinks] = chains_m[ichain, isort]
            chains_i[ichain, :nlinks] = chains_i[ichain, isort]

        if nlinks == nlinks_cap:
            accept = (chains_i[:, nlinks_cap-1] != iiter).astype(num.int)
            nlinks -= 1
        else:
            accept = num.ones(1 + problem.nchains, dtype=num.int)

        if rundir:
            problem.dump_problem_data(
                rundir, x, ms, ns, accept,
                ichain_choice if ichain_choice is not None else -1,
                ibase if ibase is not None else -1)

        accept_sum += accept
        accept_hist[iiter] = num.sum(accept)

        xhist[iiter, :] = x

        bxs = xhist[chains_i[:, :nlinks].ravel(), :]
        gxs = xhist[chains_i[0, :nlinks], :]
        gms = chains_m[0, :nlinks]

        if nlinks > (nlinks_cap-1)/2:
            # mean and std of all bootstrap ensembles together
            mbx = num.mean(bxs, axis=0)
            sbx = num.std(bxs, axis=0)

            # mean and std of global configuration
            mgx = num.mean(gxs, axis=0)
            sgx = num.std(gxs, axis=0)

            # best in global configuration
            bgx = xhist[chains_i[0, 0], :]

            covs = []
            mxs = []
            local_sxs = []

            for i in xrange(1 + problem.nchains):
                xs = xhist[chains_i[i, :nlinks], :]
                mx = num.mean(xs, axis=0)
                cov = num.cov(xs.T)

                mxs.append(mx)
                covs.append(cov)
                if standard_deviation_estimator == \
                        'median_density_single_chain':
                    local_sx = local_std(xs)
                    local_sxs.append(local_sx)
                elif standard_deviation_estimator == \
                        'standard_deviation_all_chains':
                    local_sxs.append(sbx)
                elif standard_deviation_estimator == \
                        'standard_deviation_single_chain':
                    sx = num.std(xs, axis=0)
                    local_sxs.append(sx)
                else:
                    assert False, 'invalid standard_deviation_estimator choice'

            state.parameter_sets['BS mean'][:-1-problem.ndependants] = mbx
            state.parameter_sets['BS std'][:-1-problem.ndependants] = sbx
            state.parameter_sets['Global mean'][:-1-problem.ndependants] = mgx
            state.parameter_sets['Global std'][:-1-problem.ndependants] = sgx
            state.parameter_sets['Global best'][:-1-problem.ndependants] = bgx

            state.parameter_sets['Global mean'][-1] = num.mean(gms)
            state.parameter_sets['Global std'][-1] = num.std(gms)
            state.parameter_sets['Global best'][-1] = num.min(gms)

            state.iiter = iiter + 1
            state.extra_text =\
                'Phase: %s (factor %d); ntries %d, ntries_preconstrain %d'\
                % (phase, factor, ntries_sample, ntries_preconstrain)

        if 'state' in status:
            notifier.emit('state', state)

        if 'matrix' in status:
            lines = []
            matrix = (chains_i[:, :30] % 94 + 32).T
            for row in matrix[::-1]:
                lines.append(''.join(chr(xxx) for xxx in row))
            print '\n'.join(lines)

        if plot and plot.want_to_update(iiter):
            plot.update(
                xhist[:iiter+1, :],
                chains_i[:, :nlinks],
                ibase,
                ichain_choice,
                local_sxs,
                factor)

        iiter += 1

    if plot:
        plot.finish()


class HighScoreOptimizer(Optimizer):

    def __init__(self, kwargs):
        Optimizer.__init__(self)
        self._kwargs = kwargs

    def solve(
            self, problem, rundir=None, status=(), plot=None, xs_inject=None,
            notifier=None):

        solve(
            problem,
            rundir=rundir,
            status=status,
            plot=plot,
            xs_inject=xs_inject,
            notifier=notifier,
            state=self.state,
            **self._kwargs)


class SamplerDistributionChoice(StringChoice):
    choices = ['multivariate_normal', 'normal']


class StandardDeviationEstimatorChoice(StringChoice):
    choices = [
        'median_density_single_chain',
        'standard_deviation_all_chains',
        'standard_deviation_single_chain']


class HighScoreOptimizerConfig(OptimizerConfig):
    niter_uniform = Int.T(default=1000)
    niter_transition = Int.T(default=0)
    niter_explorative = Int.T(default=10000)
    niter_non_explorative = Int.T(default=0)
    sampler_distribution = SamplerDistributionChoice.T(
        default='multivariate_normal')
    standard_deviation_estimator = StandardDeviationEstimatorChoice.T(
        default='median_density_single_chain')
    scatter_scale_transition = Float.T(default=2.0)
    scatter_scale = Float.T(default=1.0)
    chain_length_factor = Float.T(default=8.0)
    compensate_excentricity = Bool.T(default=True)

    def get_optimizer_kwargs(self):
        return dict(
            niter_uniform=self.niter_uniform,
            niter_transition=self.niter_transition,
            niter_explorative=self.niter_explorative,
            niter_non_explorative=self.niter_non_explorative,
            sampler_distribution=self.sampler_distribution,
            standard_deviation_estimator=self.standard_deviation_estimator,
            scatter_scale_transition=self.scatter_scale_transition,
            scatter_scale=self.scatter_scale,
            chain_length_factor=self.chain_length_factor,
            compensate_excentricity=self.compensate_excentricity)

    def get_optimizer(self):
        return HighScoreOptimizer(self.get_optimizer_kwargs())


__all__ = '''
    HighScoreOptimizer
    SamplerDistributionChoice
    StandardDeviationEstimatorChoice
    HighScoreOptimizerConfig
'''.split()
