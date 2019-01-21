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
from shapely.geometry import Polygon

from grond.meta import GrondError, Forbidden
from grond.problems.base import ModelHistory
from grond.optimisers.base import Optimiser, OptimiserConfig, BadProblem, \
    OptimiserStatus

guts_prefix = 'grond'

logger = logging.getLogger('grond.optimisers.highscore.optimiser')
d2r = math.pi / 180.

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


class BootstrapTypeChoice(StringChoice):
    choices = ['bayesian', 'classic']


class SamplerPhase(Object):
    niterations = Int.T()
    ntries_preconstrain_limit = Int.T(default=1000)

    def get_raw_sample(self, problem, iiter, chains):
        raise NotImplementedError

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
        intersect = True
        while intersect is True:
            pars = problem.random_uniform(xbounds)
            nsources = 2
            sources = []
            if nsources == 2:
                for i in range(nsources):
                    source = problem.get_source(pars, i)
                    sources.append(source)
            src1 = sources[0].outline('xy')
            src2 = sources[1].outline('xy')
            p1 = Polygon(src1)
            p2 = Polygon(src2)
            if not p1.intersects(p2):
                intersect = False
            else:
                intersect = True
                print('intersection, uniform phase redraw')

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
        intersect = True
        ichain_choice = num.argmin(chains.accept_sum)

        while intersect is True:
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
            pars = x
            nsources = 2
            sources = []
            if nsources == 2:
                for i in range(nsources):
                    source = problem.get_source(pars, i)
                    sources.append(source)
            src1 = sources[0].outline('xy')
            src2 = sources[1].outline('xy')
            p1 = Polygon(src1)
            p2 = Polygon(src2)
            if not p1.intersects(p2):
                intersect = False
            else:
                intersect = True
                print('intersection, directed phase redraw')

            return x


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
            (self.nchains, nlinks_cap), num.float)
        self.chains_i = num.zeros(
            (self.nchains, nlinks_cap), num.int)
        self.nlinks = 0
        self.accept_sum = num.zeros(self.nchains, dtype=num.int)
        self.nread = 0
        history.add_listener(self)

    def goto(self, n=None):
        if n is None:
            n = self.history.nmodels

        n = min(self.history.nmodels, n)

        assert self.nread <= n

        while self.nread < n:
            gbms = self.history.bootstrap_misfits[self.nread, :]

            self.chains_m[:, self.nlinks] = gbms
            self.chains_i[:, self.nlinks] = n-1
            nbootstrap = self.chains_m.shape[0]

            self.nlinks += 1
            chains_m = self.chains_m
            chains_i = self.chains_i

            for ichain in range(nbootstrap):
                isort = num.argsort(chains_m[ichain, :self.nlinks])
                chains_m[ichain, :self.nlinks] = chains_m[ichain, isort]
                chains_i[ichain, :self.nlinks] = chains_i[ichain, isort]

            if self.nlinks == self.nlinks_cap:
                accept = (
                    chains_i[:, self.nlinks_cap-1] != n-1).astype(num.int)
                self.nlinks -= 1
            else:
                accept = num.ones(self.nchains, dtype=num.int)

            self.accept_sum += accept
            self.nread += 1

    def append(self, iiter, model, misfits):
        self.goto(iiter)

    def extend(self, ioffset, n, models, misfits):
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


class HighScoreOptimiser(Optimiser):
    '''Monte-Carlo-based directed search optimisation with bootstrap.'''

    sampler_phases = List.T(SamplerPhase.T())
    chain_length_factor = Float.T(default=8.)
    nbootstrap = Int.T(default=100)
    bootstrap_type = BootstrapTypeChoice.T(default='bayesian')
    bootstrap_seed = Int.T(default=23)

    def __init__(self, **kwargs):
        Optimiser.__init__(self, **kwargs)
        self._bootstrap_weights = None
        self._bootstrap_residuals = None
        self._status_chains = None
        self.rstate = num.random.RandomState(self.bootstrap_seed)

    def init_bootstraps(self, problem):
        self.init_bootstrap_weights(problem)
        self.init_bootstrap_residuals(problem)

    def init_bootstrap_weights(self, problem):
        logger.info('Initializing Bayesian bootstrap weights')
        bootstrap_targets = set([t for t in problem.targets
                                 if t.can_bootstrap_weights])

        ws = make_bayesian_weights(
            self.nbootstrap,
            nmisfits=problem.nmisfits,
            rstate=self.rstate)

        imf = 0
        for it, t in enumerate(bootstrap_targets):
            t.set_bootstrap_weights(ws[:, imf:imf+t.nmisfits])
            imf += t.nmisfits

        for t in set(problem.targets) - bootstrap_targets:
            t.set_bootstrap_weights(
                num.ones((self.nbootstrap, t.nmisfits)))

    def init_bootstrap_residuals(self, problem):
        logger.info('Initializing Bayesian bootstrap residuals')
        residual_targets = set([t for t in problem.targets
                                if t.can_bootstrap_residuals])

        for t in residual_targets:
            t.init_bootstrap_residuals(self.nbootstrap, rstate=self.rstate)

        for t in set(problem.targets) - residual_targets:
            t.set_bootstrap_residuals(num.zeros((self.nbootstrap, t.nmisfits)))

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
        for phase in self.sampler_phases:
            if iiter < niter + phase.niterations:
                return phase, iiter - niter

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
                    iiter, niter,
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
            phase, iiter_phase = self.get_sampler_phase(iiter)
            self.log_progress(problem, iiter, niter, phase, iiter_phase)

            x = phase.get_sample(problem, iiter_phase, chains)

            if isbad_mask is not None and num.any(isbad_mask):
                isok_mask = num.logical_not(isbad_mask)
            else:
                isok_mask = None

            misfits = problem.misfits(x, mask=isok_mask)
            bootstrap_misfits = problem.combine_misfits(
                misfits,
                extra_weights=self.get_bootstrap_weights(problem),
                extra_residuals=self.get_bootstrap_residuals(problem))

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
                    'problem %s: all target misfit values are NaN'
                    % problem.name)

            history.append(x, misfits, bootstrap_misfits)

    @property
    def niterations(self):
        return sum([ph.niterations for ph in self.sampler_phases])

    def get_status(self, history):
        sparks = u'\u2581\u2582\u2583\u2584\u2585\u2586\u2587\u2588'

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

        phase = self.get_sampler_phase(history.nmodels-1)[0]

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

        def spark_plot(data, bins):
            hist, _ = num.histogram(data, bins)
            hist_max = num.max(hist)
            if hist_max == 0.0:
                hist_max = 1.0
            hist = hist / hist_max
            vec = num.digitize(hist, num.linspace(0., 1., len(sparks)))
            return ''.join([sparks[b-1] for b in vec])

        return OptimiserStatus(
            row_names=row_names,
            column_data=OrderedDict(
                zip(['BS mean', 'BS std',
                     'Glob mean', 'Glob std', 'Glob best'],
                    [bs_mean, bs_std, glob_mean, glob_std, glob_best])),
            extra_header=u'Optimiser phase: %s, exploring %d BS chains\n'
                         u'Global chain misfit distribution: \u2080%s\xb9'
                         % (phase.__class__.__name__,
                            chains.nchains,
                            spark_plot(
                                glob_misfits,
                                num.linspace(0., 1., 25))
                            ))

    def get_movie_maker(
            self, problem, history, xpar_name, ypar_name, movie_filename):

        from . import plot
        return plot.HighScoreOptimiserPlot(
            self, problem, history, xpar_name, ypar_name, movie_filename)


class HighScoreOptimiserConfig(OptimiserConfig):

    sampler_phases = List.T(
        SamplerPhase.T(),
        default=[UniformSamplerPhase(niterations=1000),
                 DirectedSamplerPhase(niterations=5000)])
    chain_length_factor = Float.T(default=8.)
    nbootstrap = Int.T(default=100)

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
