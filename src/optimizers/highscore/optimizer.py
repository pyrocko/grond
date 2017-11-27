from __future__ import print_function
import math
import os.path as op
import os
import logging

import numpy as num

from pyrocko.guts import StringChoice, Int, Float, Object, List
from pyrocko.guts_array import Array

from grond.meta import GrondError, Forbidden

from ..base import Optimizer, OptimizerConfig, BadProblem

from grond.problems.base import ModelHistory

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
    print(num.sort(num.sum(distances_sqr_all < 1.0, axis=1)))
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


class SamplerDistributionChoice(StringChoice):
    choices = ['multivariate_normal', 'normal']


class StandardDeviationEstimatorChoice(StringChoice):
    choices = [
        'median_density_single_chain',
        'standard_deviation_all_chains',
        'standard_deviation_single_chain']


class SamplerStartingPointChoice(StringChoice):
    choices = ['excentricity_compensated', 'random', 'mean']


class SamplerPhase(Object):
    niterations = Int.T()
    ntries_preconstrain_limit = Int.T(default=1000)

    def get_sample(self, problem, iiter, chains):

        assert 0 <= iiter < self.niterations

        ntries_preconstrain = 0
        for ntries_preconstrain in range(self.ntries_preconstrain_limit):
            try:
                return problem.preconstrain(
                    self.get_raw_sample(problem, iiter, chains))

            except Forbidden:
                pass

        raise GrondError(
            'could not find any suitable candidate sample within %i tries' % (
                self.ntries_preconstrain_limit))


class InjectionSamplerPhase(SamplerPhase):
    xs_inject = Array.T(dtype=num.float, shape=(None, None))

    def get_raw_sample(self, problem, iiter, chains):
        return self.xs_inject[iiter, :]


class UniformSamplerPhase(SamplerPhase):

    def get_raw_sample(self, problem, iiter, chains):
        xbounds = problem.get_parameter_bounds()
        return problem.random_uniform(xbounds)


class DirectedSamplerPhase(SamplerPhase):
    scatter_scale = Float.T(optional=True)
    scatter_scale_begin = Float.T(optional=True)
    scatter_scale_end = Float.T(optional=True)
    starting_point = SamplerStartingPointChoice.T(
        default='excentricity_compensated')

    sampler_distribution = SamplerDistributionChoice.T(
        default='normal')

    standard_deviation_estimator = StandardDeviationEstimatorChoice.T(
        default='median_density_single_chain')

    ntries_sample_limit = Int.T(default=1000)

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

    def get_raw_sample(self, problem, iiter, chains):

        factor = self.get_scatter_scale_factor(iiter)
        npar = problem.nparameters
        pnames = problem.parameter_names
        xbounds = problem.get_parameter_bounds()

        ichain_choice = num.argmin(chains.accept_sum)

        if self.starting_point == 'excentricity_compensated':
            models = chains.models(ichain_choice)
            ilink_choice = excentricity_compensated_choice(
                models,
                chains.standard_deviation_models(
                    ichain_choice, self.standard_deviation_estimator),
                2.)

            xchoice = chains.model(ichain_choice, ilink_choice)

        elif self.starting_point == 'random':
            ilink_choice = num.random.randint(0, chains.nlinks)
            xchoice = chains.model(ichain_choice, ilink_choice)

        elif self.starting_point == 'mean':
            xchoice = chains.mean_model(ichain_choice)

        else:
            assert False, 'invalid starting_point choice: %s' % (
                self.starting_point)

        ntries_sample = 0
        if self.sampler_distribution == 'normal':
            x = num.zeros(npar, dtype=num.float)
            sx = chains.standard_deviation_models(
                ichain_choice, self.standard_deviation_estimator)

            for ipar in range(npar):
                ntries = 0
                while True:
                    if sx[ipar] > 0.:
                        v = num.random.normal(
                            xchoice[ipar],
                            factor*sx[ipar])
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
                    xchoice, factor**2 * chains.cov(ichain_choice))

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

        return x


class Chains(object):
    def __init__(self, problem, history, nchains, nlinks_cap):

        self.problem = problem
        self.history = history
        self.nchains = nchains
        self.nlinks_cap = nlinks_cap
        self.chains_m = num.zeros(
            (nchains, nlinks_cap), num.float)
        self.chains_i = num.zeros(
            (nchains, nlinks_cap), num.int)
        self.nlinks = 0
        self.accept_sum = num.zeros(nchains, dtype=num.int)
        history.add_listener(self)

    def append(self, iiter, model, misfits):

        gm = self.problem.global_misfit(misfits)
        bms = self.problem.bootstrap_misfit(misfits, self.nchains - 1)
        gbms = num.concatenate(([gm], bms))

        self.chains_m[:, self.nlinks] = gbms
        self.chains_i[:, self.nlinks] = iiter

        self.nlinks += 1
        chains_m = self.chains_m
        chains_i = self.chains_i

        for ichain in range(chains_m.shape[0]):
            isort = num.argsort(chains_m[ichain, :self.nlinks])
            chains_m[ichain, :self.nlinks] = chains_m[ichain, isort]
            chains_i[ichain, :self.nlinks] = chains_i[ichain, isort]

        if self.nlinks == self.nlinks_cap:
            accept = (chains_i[:, self.nlinks_cap-1] != iiter).astype(num.int)
            self.nlinks -= 1
        else:
            accept = num.ones(self.nchains, dtype=num.int)

        self.accept_sum += accept

        return accept

    def extend(self, ioffset, n, models, misfits):
        for i in range(ioffset, ioffset+n):
            self.append(i, models[i-ioffset, :], misfits[i-ioffset, :, :])

    def models(self, ichain=None):
        if ichain is not None:
            return self.history.models[
                self.chains_i[ichain, :self.nlinks], :]
        else:
            return self.history.models[
                self.chains_i[:, :self.nlinks].ravel(), :]

    def model(self, ichain, ilink):
        return self.history.models[self.chains_i[ichain, ilink], :]

    def mean_model(self, ichain):
        xs = self.models(ichain)
        return num.mean(xs, axis=0)

    def standard_deviation_models(self, ichain, estimator):

        if estimator == 'median_density_single_chain':
            xs = self.models(ichain)
            return local_std(xs)
        elif estimator == 'standard_deviation_all_chains':
            bxs = self.models()
            return num.std(bxs, axis=0)
        elif estimator == 'standard_deviation_single_chain':
            xs = self.models(ichain)
            return num.std(xs, axis=0)
        else:
            assert False, 'invalid standard_deviation_estimator choice'

    def covariance_models(self, ichain):
        xs = self.models(ichain)
        return num.cov(xs.T)


class HighScoreOptimizer(Optimizer):

    sampler_phases = List.T(SamplerPhase.T())
    chain_length_factor = Float.T(default=8.)
    nbootstrap = Int.T(default=10)

    def optimize(self, problem, rundir=None):

        if rundir is not None:
            self.dump(filename=op.join(rundir, 'optimizer.yaml'))

        niter = sum(phase.niterations for phase in self.sampler_phases)

        iiter = 0

        isbad_mask = None

        nlinks_cap = int(round(
            self.chain_length_factor * problem.nparameters + 1))

        history = ModelHistory(problem, path=rundir, mode='w')
        chains = Chains(problem, history, self.nbootstrap+1, nlinks_cap)

        phases = list(self.sampler_phases)
        phase = phases.pop(0)
        iiter_phase_start = 0
        while iiter < niter:

            if iiter - iiter_phase_start == phase.niterations:
                phase = phases.pop(0)
                iiter_phase_start = iiter

            iiter_phase = iiter - iiter_phase_start

            x = phase.get_sample(problem, iiter_phase, chains)

            if isbad_mask is not None and num.any(isbad_mask):
                isok_mask = num.logical_not(isbad_mask)
            else:
                isok_mask = None

            misfits = problem.evaluate(x, mask=isok_mask)

            isbad_mask_new = num.isnan(misfits[:, 0])
            if isbad_mask is not None and num.any(
                    isbad_mask != isbad_mask_new):

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
                    'problem %s: all target misfit values are NaN'
                    % problem.name)

            history.append(x, misfits)

            iiter += 1


class HighScoreOptimizerConfig(OptimizerConfig):

    sampler_phases = List.T(
        SamplerPhase.T(),
        default=[UniformSamplerPhase(niterations=1000),
                 DirectedSamplerPhase(niterations=5000)])
    chain_length_factor = Float.T(default=8.)
    nbootstrap = Int.T(default=10)

    def get_optimizer(self):
        return HighScoreOptimizer(
            sampler_phases=list(self.sampler_phases),
            chain_length_factor=self.chain_length_factor,
            nbootstrap=self.nbootstrap)


def load_optimizer_history(dirname, problem):
    fn = op.join(dirname, 'accepted')
    with open(fn, 'r') as f:
        nmodels = os.fstat(f.fileno()).st_size // (problem.nbootstrap+1)
        data1 = num.fromfile(
            f,
            dtype='<i1',
            count=nmodels*(problem.nbootstrap+1)).astype(num.bool)

    accepted = data1.reshape((nmodels, problem.nbootstrap+1))

    fn = op.join(dirname, 'choices')
    with open(fn, 'r') as f:
        data2 = num.fromfile(
            f,
            dtype='<i8',
            count=nmodels*2).astype(num.int64)

    ibootstrap_choices, imodel_choices = data2.reshape((nmodels, 2)).T
    return ibootstrap_choices, imodel_choices, accepted


__all__ = '''
    SamplerDistributionChoice
    StandardDeviationEstimatorChoice
    SamplerPhase
    InjectionSamplerPhase
    UniformSamplerPhase
    DirectedSamplerPhase
    Chains
    HighScoreOptimizerConfig
    HighScoreOptimizer
'''.split()
