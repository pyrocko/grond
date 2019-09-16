from __future__ import print_function
import math
import os.path as op
import os
import logging
import time
import numpy as num
from collections import OrderedDict

from pyrocko.guts import StringChoice, Int, Float, Object, List
from pyrocko.guts_array import Array

from grond.meta import GrondError, Forbidden, has_get_plot_classes
from grond.problems.base import ModelHistory
from grond.optimisers.base import Optimiser, OptimiserConfig, BadProblem, \
    OptimiserStatus

guts_prefix = 'grond'

logger = logging.getLogger('grond.optimisers.highscore.optimiser')


def nextpow2(i):
    return 2**int(math.ceil(math.log(i)/math.log(2.)))


def excentricity_compensated_probabilities(xs, sbx, factor):
    inonflat = num.where(sbx != 0.0)[0]
    scale = num.zeros_like(sbx)
    scale[inonflat] = 1.0 / (sbx[inonflat] * (factor if factor != 0. else 1.0))
    distances_sqr_all = num.sum(
        ((xs[num.newaxis, :, :] - xs[:, num.newaxis, :]) *
         scale[num.newaxis, num.newaxis, :])**2, axis=2)
    probabilities = 1.0 / num.sum(distances_sqr_all < 1.0, axis=1)
    # print(num.sort(num.sum(distances_sqr_all < 1.0, axis=1)))
    probabilities /= num.sum(probabilities)
    return probabilities


def excentricity_compensated_choice(xs, sbx, factor, rstate):
    probabilities = excentricity_compensated_probabilities(
        xs, sbx, factor)
    r = rstate.random_sample()
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


class BootstrapTypeChoice(StringChoice):
    choices = ['bayesian', 'classic']


def fnone(i):
    return i if i is not None else -1


class Sample(Object):

    '''Sample model with context about how it was generated.'''

    model = Array.T(shape=(None,), dtype=num.float, serialize_as='list')
    iphase = Int.T(optional=True)
    ichain_base = Int.T(optional=True)
    ilink_base = Int.T(optional=True)
    imodel_base = Int.T(optional=True)

    def preconstrain(self, problem):
        self.model = problem.preconstrain(self.model)

    def pack_context(self):
        i = num.zeros(4, dtype=num.int)
        i[:] = (
            fnone(self.iphase),
            fnone(self.ichain_base),
            fnone(self.ilink_base),
            fnone(self.imodel_base))

        return i


class SamplerPhase(Object):
    niterations = Int.T(
        help='Number of iteration for this phase.')
    ntries_preconstrain_limit = Int.T(
        default=1000,
        help='Tries to find a valid preconstrained sample.')
    seed = Int.T(
        optional=True,
        help='Random state seed.')

    def __init__(self, *args, **kwargs):
        Object.__init__(self, *args, **kwargs)
        self._rstate = None

    def get_rstate(self):
        if self._rstate is None:
            self._rstate = num.random.RandomState(self.seed)

        return self._rstate

    def get_raw_sample(self, problem, iiter, chains):
        raise NotImplementedError

    def get_sample(self, problem, iiter, chains):
        assert 0 <= iiter < self.niterations

        ntries_preconstrain = 0
        for ntries_preconstrain in range(self.ntries_preconstrain_limit):
            try:
                sample = self.get_raw_sample(problem, iiter, chains)
                sample.preconstrain(problem)
                return sample

            except Forbidden:
                pass

        raise GrondError(
            'could not find any suitable candidate sample within %i tries' % (
                self.ntries_preconstrain_limit))


class InjectionSamplerPhase(SamplerPhase):
    xs_inject = Array.T(
        dtype=num.float, shape=(None, None),
        help='Array with the reference model.')

    def get_raw_sample(self, problem, iiter, chains):
        return Sample(model=self.xs_inject[iiter, :])


class UniformSamplerPhase(SamplerPhase):

    def get_raw_sample(self, problem, iiter, chains):
        xbounds = problem.get_parameter_bounds()
        return Sample(model=problem.random_uniform(xbounds, self.get_rstate()))


class DirectedSamplerPhase(SamplerPhase):
    scatter_scale = Float.T(
        optional=True,
        help='Scales search radius around the current `highscore` models')
    scatter_scale_begin = Float.T(
        optional=True,
        help='Scaling factor at beginning of the phase.')
    scatter_scale_end = Float.T(
        optional=True,
        help='Scaling factor at the end of the directed phase.')
    starting_point = SamplerStartingPointChoice.T(
        default='excentricity_compensated',
        help='Tunes to the center value of the sampler distribution.'
             'May increase the likelihood to draw a highscore member model'
             ' off-center to the mean value')

    sampler_distribution = SamplerDistributionChoice.T(
        default='normal',
        help='Distribution new models are drawn from.')

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
        rstate = self.get_rstate()
        factor = self.get_scatter_scale_factor(iiter)
        npar = problem.nparameters
        pnames = problem.parameter_names
        xbounds = problem.get_parameter_bounds()

        ilink_choice = None
        ichain_choice = num.argmin(chains.accept_sum)

        if self.starting_point == 'excentricity_compensated':
            models = chains.models(ichain_choice)
            ilink_choice = excentricity_compensated_choice(
                models,
                chains.standard_deviation_models(
                    ichain_choice, self.standard_deviation_estimator),
                2., rstate)

            xchoice = chains.model(ichain_choice, ilink_choice)

        elif self.starting_point == 'random':
            ilink_choice = rstate.randint(0, chains.nlinks)
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
                        v = rstate.normal(
                            xchoice[ipar],
                            factor*sx[ipar])
                    else:
                        v = xchoice[ipar]

                    if xbounds[ipar, 0] <= v and \
                            v <= xbounds[ipar, 1]:

                        break

                    if ntries > self.ntries_sample_limit:
                        logger.warning(
                            'failed to produce a suitable '
                            'candidate sample from normal '
                            'distribution for parameter \'%s\''
                            '- drawing from uniform instead.' %
                            pnames[ipar])
                        v = rstate.uniform(xbounds[ipar, 0],
                                           xbounds[ipar, 1])
                        break

                    ntries += 1

                x[ipar] = v

        elif self.sampler_distribution == 'multivariate_normal':
            ok_mask_sum = num.zeros(npar, dtype=num.int)
            while True:
                ntries_sample += 1
                xcandi = rstate.multivariate_normal(
                    xchoice, factor**2 * chains.cov(ichain_choice))

                ok_mask = num.logical_and(
                    xbounds[:, 0] <= xcandi, xcandi <= xbounds[:, 1])

                if num.all(ok_mask):
                    break

                ok_mask_sum += ok_mask

                if ntries_sample > self.ntries_sample_limit:
                    logger.warning(
                        'failed to produce a suitable candidate '
                        'sample from multivariate normal '
                        'distribution, (%s) - drawing from uniform instead' %
                        ', '.join('%s:%i' % xx for xx in
                                  zip(pnames, ok_mask_sum)))
                    xbounds = problem.get_parameter_bounds()
                    xcandi = problem.random_uniform(xbounds, rstate)
                    break

            x = xcandi

        imodel_base = None
        if ilink_choice is not None:
            imodel_base = chains.imodel(ichain_choice, ilink_choice)

        return Sample(
            model=x,
            ichain_base=ichain_choice,
            ilink_base=ilink_choice,
            imodel_base=imodel_base)


def make_bayesian_weights(nbootstrap, nmisfits,
                          type='bayesian', rstate=None):
    ws = num.zeros((nbootstrap, nmisfits))
    if rstate is None:
        rstate = num.random.RandomState()

    for ibootstrap in range(nbootstrap):
        if type == 'classic':
            ii = rstate.randint(0, nmisfits, size=nmisfits)
            ws[ibootstrap, :] = num.histogram(
                ii, nmisfits, (-0.5, nmisfits - 0.5))[0]
        elif type == 'bayesian':
            f = rstate.uniform(0., 1., size=nmisfits+1)
            f[0] = 0.
            f[-1] = 1.
            f = num.sort(f)
            g = f[1:] - f[:-1]
            ws[ibootstrap, :] = g * nmisfits
        else:
            assert False
    return ws


class Chains(object):
    def __init__(
            self, problem, history, nchains, nlinks_cap):

        self.problem = problem
        self.history = history
        self.nchains = nchains
        self.nlinks_cap = nlinks_cap
        self.chains_m = num.zeros(
            (self.nchains, nlinks_cap), dtype=num.float)
        self.chains_i = num.zeros(
            (self.nchains, nlinks_cap), dtype=num.int)
        self.nlinks = 0
        self.nread = 0

        self.accept_sum = num.zeros(self.nchains, dtype=num.int)
        self._acceptance_history = num.zeros(
            (self.nchains, 1024), dtype=num.bool)

        history.add_listener(self)

    def goto(self, n=None):
        if n is None:
            n = self.history.nmodels

        n = min(self.history.nmodels, n)

        assert self.nread <= n

        while self.nread < n:
            nread = self.nread
            gbms = self.history.bootstrap_misfits[nread, :]

            self.chains_m[:, self.nlinks] = gbms
            self.chains_i[:, self.nlinks] = nread
            nbootstrap = self.chains_m.shape[0]

            self.nlinks += 1
            chains_m = self.chains_m
            chains_i = self.chains_i

            for ichain in range(nbootstrap):
                isort = num.argsort(chains_m[ichain, :self.nlinks])
                chains_m[ichain, :self.nlinks] = chains_m[ichain, isort]
                chains_i[ichain, :self.nlinks] = chains_i[ichain, isort]

            if self.nlinks == self.nlinks_cap:
                accept = (chains_i[:, self.nlinks_cap-1] != nread) \
                    .astype(num.bool)
                self.nlinks -= 1
            else:
                accept = num.ones(self.nchains, dtype=num.bool)

            self._append_acceptance(accept)
            self.accept_sum += accept
            self.nread += 1

    def load(self):
        return self.goto()

    def extend(self, ioffset, n, models, misfits, sampler_contexts):
        self.goto(ioffset + n)

    def indices(self, ichain):
        if ichain is not None:
            return self.chains_i[ichain, :self.nlinks]
        else:
            return self.chains_i[:, :self.nlinks].ravel()

    def models(self, ichain=None):
        return self.history.models[self.indices(ichain), :]

    def model(self, ichain, ilink):
        return self.history.models[self.chains_i[ichain, ilink], :]

    def imodel(self, ichain, ilink):
        return self.chains_i[ichain, ilink]

    def misfits(self, ichain=0):
        return self.chains_m[ichain, :self.nlinks]

    def misfit(self, ichain, ilink):
        assert ilink < self.nlinks
        return self.chains_m[ichain, ilink]

    def mean_model(self, ichain=None):
        xs = self.models(ichain)
        return num.mean(xs, axis=0)

    def best_model(self, ichain=0):
        xs = self.models(ichain)
        return xs[0]

    def best_model_misfit(self, ichain=0):
        return self.chains_m[ichain, 0]

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

    @property
    def acceptance_history(self):
        return self._acceptance_history[:, :self.nread]

    def _append_acceptance(self, acceptance):
        if self.nread >= self._acceptance_history.shape[1]:
            new_buf = num.zeros(
                (self.nchains, nextpow2(self.nread+1)), dtype=num.bool)
            new_buf[:, :self._acceptance_history.shape[1]] = \
                self._acceptance_history
            self._acceptance_history = new_buf
        self._acceptance_history[:, self.nread] = acceptance


@has_get_plot_classes
class HighScoreOptimiser(Optimiser):
    '''Monte-Carlo-based directed search optimisation with bootstrap.'''

    sampler_phases = List.T(SamplerPhase.T())
    chain_length_factor = Float.T(default=8.)
    nbootstrap = Int.T(default=100)
    bootstrap_type = BootstrapTypeChoice.T(default='bayesian')
    bootstrap_seed = Int.T(default=23)

    SPARKS = u'\u2581\u2582\u2583\u2584\u2585\u2586\u2587\u2588'
    ACCEPTANCE_AVG_LEN = 100

    def __init__(self, **kwargs):
        Optimiser.__init__(self, **kwargs)
        self._bootstrap_weights = None
        self._bootstrap_residuals = None
        self._correlated_weights = None
        self._status_chains = None
        self._rstate_bootstrap = None

    def get_rstate_bootstrap(self):
        if self._rstate_bootstrap is None:
            self._rstate_bootstrap = num.random.RandomState(
                self.bootstrap_seed)

        return self._rstate_bootstrap

    def init_bootstraps(self, problem):
        self.init_bootstrap_weights(problem)
        self.init_bootstrap_residuals(problem)

    def init_bootstrap_weights(self, problem):
        logger.info('Initializing Bayesian bootstrap weights.')

        nmisfits_w = sum(
            t.nmisfits for t in problem.targets if t.can_bootstrap_weights)

        ws = make_bayesian_weights(
            self.nbootstrap,
            nmisfits=nmisfits_w,
            rstate=self.get_rstate_bootstrap())

        imf = 0
        for t in problem.targets:
            if t.can_bootstrap_weights:
                t.set_bootstrap_weights(ws[:, imf:imf+t.nmisfits])
                imf += t.nmisfits
            else:
                t.set_bootstrap_weights(
                    num.ones((self.nbootstrap, t.nmisfits)))

    def init_bootstrap_residuals(self, problem):
        logger.info('Initializing Bayesian bootstrap residuals.')

        for t in problem.targets:
            if t.can_bootstrap_residuals:
                t.init_bootstrap_residuals(
                    self.nbootstrap, rstate=self.get_rstate_bootstrap(),
                    nthreads=self._nthreads)
            else:
                t.set_bootstrap_residuals(
                    num.zeros((self.nbootstrap, t.nmisfits)))

    def get_bootstrap_weights(self, problem):
        if self._bootstrap_weights is None:
            try:
                problem.targets[0].get_bootstrap_weights()
            except Exception:
                self.init_bootstraps(problem)

            bootstrap_weights = num.hstack(
                [t.get_bootstrap_weights()
                 for t in problem.targets])

            self._bootstrap_weights = num.vstack((
                num.ones((1, problem.nmisfits)),
                bootstrap_weights))

        return self._bootstrap_weights

    def get_bootstrap_residuals(self, problem):
        if self._bootstrap_residuals is None:
            try:
                problem.targets[0].get_bootstrap_residuals()
            except Exception:
                self.init_bootstraps(problem)

            bootstrap_residuals = num.hstack(
                [t.get_bootstrap_residuals()
                 for t in problem.targets])

            self._bootstrap_residuals = num.vstack((
                num.zeros((1, problem.nmisfits)),
                bootstrap_residuals))

        return self._bootstrap_residuals

    def get_correlated_weights(self, problem):
        if self._correlated_weights is None:
            corr = dict()
            misfit_idx = num.cumsum(
                [0.] + [t.nmisfits for t in problem.targets], dtype=num.int)

            for it, target in enumerate(problem.targets):
                weights = target.get_correlated_weights(
                    nthreads=self._nthreads)
                if weights is None:
                    continue
                corr[misfit_idx[it]] = weights

            self._correlated_weights = corr

        return self._correlated_weights

    @property
    def nchains(self):
        return self.nbootstrap + 1

    def chains(self, problem, history):
        nlinks_cap = int(round(
            self.chain_length_factor * problem.nparameters + 1))

        return Chains(
            problem, history,
            nchains=self.nchains, nlinks_cap=nlinks_cap)

    def get_sampler_phase(self, iiter):
        niter = 0
        for iphase, phase in enumerate(self.sampler_phases):
            if iiter < niter + phase.niterations:
                return iphase, phase, iiter - niter

            niter += phase.niterations

        assert False, 'sample out of bounds'

    def log_progress(self, problem, iiter, niter, phase, iiter_phase):
        t = time.time()
        if self._tlog_last < t - 10. \
                or iiter_phase == 0 \
                or iiter_phase == phase.niterations - 1:

            logger.info(
                '%s at %i/%i (%s, %i/%i)' % (
                    problem.name,
                    iiter+1, niter,
                    phase.__class__.__name__, iiter_phase, phase.niterations))

            self._tlog_last = t

    def optimise(self, problem, rundir=None):
        if rundir is not None:
            self.dump(filename=op.join(rundir, 'optimiser.yaml'))

        history = ModelHistory(problem,
                               nchains=self.nchains,
                               path=rundir, mode='w')
        chains = self.chains(problem, history)

        niter = self.niterations
        isbad_mask = None
        self._tlog_last = 0
        for iiter in range(niter):
            iphase, phase, iiter_phase = self.get_sampler_phase(iiter)
            self.log_progress(problem, iiter, niter, phase, iiter_phase)

            sample = phase.get_sample(problem, iiter_phase, chains)
            sample.iphase = iphase

            if isbad_mask is not None and num.any(isbad_mask):
                isok_mask = num.logical_not(isbad_mask)
            else:
                isok_mask = None

            misfits = problem.misfits(
                sample.model, mask=isok_mask, nthreads=self._nthreads)

            bootstrap_misfits = problem.combine_misfits(
                misfits,
                extra_weights=self.get_bootstrap_weights(problem),
                extra_residuals=self.get_bootstrap_residuals(problem),
                extra_correlated_weights=self.get_correlated_weights(problem))
            isbad_mask_new = num.isnan(misfits[:, 0])
            if isbad_mask is not None and num.any(
                    isbad_mask != isbad_mask_new):

                errmess = [
                    'problem %s: inconsistency in data availability'
                    ' at iteration %i' %
                    (problem.name, iiter)]

                for target, isbad_new, isbad in zip(
                        problem.targets, isbad_mask_new, isbad_mask):

                    if isbad_new != isbad:
                        errmess.append('  %s, %s -> %s' % (
                            target.string_id(), isbad, isbad_new))

                raise BadProblem('\n'.join(errmess))

            isbad_mask = isbad_mask_new

            if num.all(isbad_mask):
                raise BadProblem(
                    'Problem %s: all target misfit values are NaN.'
                    % problem.name)

            history.append(
                sample.model, misfits,
                bootstrap_misfits,
                sample.pack_context())

    @property
    def niterations(self):
        return sum([ph.niterations for ph in self.sampler_phases])

    def get_status(self, history):
        if self._status_chains is None:
            self._status_chains = self.chains(history.problem, history)

        self._status_chains.goto(history.nmodels)

        chains = self._status_chains
        problem = history.problem

        row_names = [p.name_nogroups for p in problem.parameters]
        row_names.append('Misfit')

        def colum_array(data):
            arr = num.full(len(row_names), fill_value=num.nan)
            arr[:data.size] = data
            return arr

        phase = self.get_sampler_phase(history.nmodels-1)[1]

        bs_mean = colum_array(chains.mean_model(ichain=None))
        bs_std = colum_array(chains.standard_deviation_models(
            ichain=None, estimator='standard_deviation_all_chains'))

        glob_mean = colum_array(chains.mean_model(ichain=0))
        glob_mean[-1] = num.mean(chains.misfits(ichain=0))

        glob_std = colum_array(chains.standard_deviation_models(
            ichain=0, estimator='standard_deviation_single_chain'))
        glob_std[-1] = num.std(chains.misfits(ichain=0))

        glob_best = colum_array(chains.best_model(ichain=0))
        glob_best[-1] = chains.best_model_misfit()

        glob_misfits = chains.misfits(ichain=0)

        acceptance_latest = chains.acceptance_history[
            :, -min(chains.acceptance_history.shape[1], self.ACCEPTANCE_AVG_LEN):]  # noqa
        acceptance_avg = acceptance_latest.mean(axis=1)

        def spark_plot(data, bins):
            hist, _ = num.histogram(data, bins)
            hist_max = num.max(hist)
            if hist_max == 0.0:
                hist_max = 1.0
            hist = hist / hist_max
            vec = num.digitize(hist, num.linspace(0., 1., len(self.SPARKS)))
            return ''.join([self.SPARKS[b-1] for b in vec])

        return OptimiserStatus(
            row_names=row_names,
            column_data=OrderedDict(
                zip(['BS mean', 'BS std',
                     'Glob mean', 'Glob std', 'Glob best'],
                    [bs_mean, bs_std, glob_mean, glob_std, glob_best])),
            extra_header=  # noqa
                u'Optimiser phase: {phase}, exploring {nchains} BS chains\n'  # noqa
                u'Global chain misfit distribution: \u2080{mf_dist}\xb9\n'
                u'Acceptance rate distribution:     \u2080{acceptance}'
                u'\u2081\u2080\u2080\ufe6a (Median {acceptance_med:.1f}%)'
                .format(
                    phase=phase.__class__.__name__,
                    nchains=chains.nchains,
                    mf_dist=spark_plot(
                        glob_misfits, num.linspace(0., 1., 25)),
                    acceptance=spark_plot(
                        acceptance_avg,
                        num.linspace(0., 1., 25)),
                    acceptance_med=num.median(acceptance_avg) * 100.
                    ))

    def get_movie_maker(
            self, problem, history, xpar_name, ypar_name, movie_filename):

        from . import plot
        return plot.HighScoreOptimiserPlot(
            self, problem, history, xpar_name, ypar_name, movie_filename)

    @classmethod
    def get_plot_classes(cls):
        from .plot import HighScoreAcceptancePlot
        plots = Optimiser.get_plot_classes()
        plots.append(HighScoreAcceptancePlot)
        return plots


class HighScoreOptimiserConfig(OptimiserConfig):

    sampler_phases = List.T(
        SamplerPhase.T(),
        default=[UniformSamplerPhase(niterations=1000),
                 DirectedSamplerPhase(niterations=5000)],
        help='Stages of the sampler: Start with uniform sampling of the model'
             ' model space and narrow down through directed sampling.')
    chain_length_factor = Float.T(
        default=8.,
        help='Controls the length of each chain: '
             'chain_length_factor * nparameters + 1')
    nbootstrap = Int.T(
        default=100,
        help='Number of bootstrap realisations to be tracked simultaneously in'
             ' the optimisation.')

    def get_optimiser(self):
        return HighScoreOptimiser(
            sampler_phases=list(self.sampler_phases),
            chain_length_factor=self.chain_length_factor,
            nbootstrap=self.nbootstrap)


def load_optimiser_history(dirname, problem):
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
    HighScoreOptimiserConfig
    HighScoreOptimiser
'''.split()
