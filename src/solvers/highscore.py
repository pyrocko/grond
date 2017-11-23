import math
import logging

import numpy as num

from collections import OrderedDict

from pyrocko.guts import StringChoice, Int, Float, Bool

from ..meta import GrondError, Forbidden
from .base import Solver, SolverConfig

guts_prefix = 'grond'

logger = logging.getLogger('grond.solver.highscore')


def excentricity_compensated_probabilities(xs, sbx, factor):
    inonflat = num.where(sbx != 0.0)[0]
    scale = num.zeros_like(sbx)
    scale[inonflat] = 1.0 / (sbx[inonflat] * (factor if factor != 0. else 1.0))
    distances_sqr_all = num.sum(
        ((xs[num.newaxis, :, :] - xs[:, num.newaxis, :]) *
         scale[num.newaxis, num.newaxis, :])**2, axis=2)
    probabilities = 1.0 / num.sum(distances_sqr_all < 1.0, axis=1)
    # print num.sort(num.sum(distances_sqr_all < 1.0, axis=1))
    probabilities /= num.sum(probabilities)
    return probabilities


def excentricity_compensated_choice(xs, sbx, factor):
    probabilities = excentricity_compensated_probabilities(
        xs, sbx, factor)
    r = num.random.random()
    ichoice = num.searchsorted(num.cumsum(probabilities), r)
    ichoice = min(ichoice, xs.shape[0]-1)
    return ichoice


def select_most_excentric(xcandidates, xs, sbx, factor):
    inonflat = num.where(sbx != 0.0)[0]
    scale = num.zeros_like(sbx)
    scale[inonflat] = 1.0 / (sbx[inonflat] * (factor if factor != 0. else 1.0))
    distances_sqr_all = num.sum(
        ((xcandidates[num.newaxis, :, :] - xs[:, num.newaxis, :]) *
         scale[num.newaxis, num.newaxis, :])**2, axis=2)
    # print num.sort(num.sum(distances_sqr_all < 1.0, axis=0))
    ichoice = num.argmin(num.sum(distances_sqr_all < 1.0, axis=0))
    return xcandidates[ichoice]


def local_std(xs):
    ssbx = num.sort(xs, axis=0)
    dssbx = num.diff(ssbx, axis=0)
    mdssbx = num.median(dssbx, axis=0)
    return mdssbx * dssbx.shape[0] / 2.6


def solve(problem,
          rundir=None,
          nbootstrap=100,
          niter_uniform=1000,
          niter_transition=1000,
          niter_explorative=10000,
          niter_non_explorative=0,
          scatter_scale_transition=2.0,
          scatter_scale=1.0,
          chain_length_factor=8.0,
          sampler_distribution='multivariate_normal',
          standard_deviation_estimator='median_density_single_chain',
          compensate_excentricity=True,
          xs_inject=None,
          status=(),
          plot=None,
          notifier=None,
          state=None):

    xbounds = num.array(problem.get_parameter_bounds(), dtype=num.float)
    npar = problem.nparameters

    nlinks_cap = int(round(chain_length_factor * npar + 1))
    chains_m = num.zeros((1 + nbootstrap, nlinks_cap), num.float)
    chains_i = num.zeros((1 + nbootstrap, nlinks_cap), num.int)
    nlinks = 0
    mbx = None

    if xs_inject is not None and xs_inject.size != 0:
        niter_inject = xs_inject.shape[0]
    else:
        niter_inject = 0

    niter = niter_inject + niter_uniform + niter_transition + \
        niter_explorative + niter_non_explorative

    iiter = 0
    sbx = None
    mxs = None
    covs = None
    local_sxs = None
    xhist = num.zeros((niter, npar))
    isbad_mask = None
    accept_sum = num.zeros(1 + nbootstrap, dtype=num.int)
    accept_hist = num.zeros(niter, dtype=num.int)
    pnames = problem.parameter_names

    state.problem_name = problem.name
    state.parameter_names = problem.parameter_names + ['Misfit']

    state.parameter_sets = OrderedDict()

    state.parameter_sets['BS mean'] = num.zeros(state.nparameters)
    state.parameter_sets['BS std'] = num.zeros(state.nparameters)
    state.parameter_sets['Global mean'] = num.zeros(state.nparameters)
    state.parameter_sets['Global std'] = num.zeros(state.nparameters)
    state.parameter_sets['Global best'] = num.zeros(state.nparameters)

    for par in state.parameter_sets.values():
        par.fill(num.nan)

    state.niter = niter

    if plot:
        plot.start(problem)

    while iiter < niter:
        ibootstrap_choice = None
        ichoice = None

        if iiter < niter_inject:
            phase = 'inject'
        elif iiter < niter_inject + niter_uniform:
            phase = 'uniform'
        elif iiter < niter_inject + niter_uniform + niter_transition:
            phase = 'transition'
        elif iiter < niter_inject + niter_uniform + niter_transition + \
                niter_explorative:
            phase = 'explorative'
        else:
            phase = 'non_explorative'

        factor = 0.0
        if phase == 'transition':
            T = float(niter_transition)
            A = scatter_scale_transition
            B = scatter_scale
            tau = T/(math.log(A) - math.log(B))
            t0 = math.log(A) * T / (math.log(A) - math.log(B))
            t = float(iiter - niter_uniform - niter_inject)
            factor = num.exp(-(t-t0) / tau)

        elif phase in ('explorative', 'non_explorative'):
            factor = scatter_scale

        ntries_preconstrain = 0
        ntries_sample = 0

        if phase == 'inject':
            x = xs_inject[iiter, :]
        else:
            while True:
                ntries_preconstrain += 1

                if mbx is None or phase == 'uniform':
                    x = problem.random_uniform(xbounds)
                else:
                    # ibootstrap_choice = num.random.randint(
                    #     0, 1 + nbootstrap)
                    ibootstrap_choice = num.argmin(accept_sum)

                    if phase in ('transition', 'explorative'):

                        if compensate_excentricity:
                            ichoice = excentricity_compensated_choice(
                                xhist[chains_i[ibootstrap_choice, :], :],
                                local_sxs[ibootstrap_choice], 2.)

                            xchoice = xhist[
                                chains_i[ibootstrap_choice, ichoice], :]

                        else:
                            ichoice = num.random.randint(0, nlinks)
                            xchoice = xhist[
                                chains_i[ibootstrap_choice, ichoice], :]
                    else:
                        xchoice = mxs[ibootstrap_choice]

                    if sampler_distribution == 'multivariate_normal':
                        ntries_sample = 0

                        ntry = 0
                        ok_mask_sum = num.zeros(npar, dtype=num.int)
                        while True:
                            ntries_sample += 1
                            vs = num.random.multivariate_normal(
                                xchoice, factor**2 * covs[ibootstrap_choice])

                            ok_mask = num.logical_and(
                                xbounds[:, 0] <= vs, vs <= xbounds[:, 1])

                            if num.all(ok_mask):
                                break

                            ok_mask_sum += ok_mask

                            if ntry > 1000:
                                raise GrondError(
                                    'failed to produce a suitable candidate '
                                    'sample from multivariate normal '
                                    'distribution, (%s)' %
                                    ', '.join('%s:%i' % xx for xx in
                                              zip(pnames, ok_mask_sum)))

                            ntry += 1

                        x = vs.tolist()

                    if sampler_distribution == 'normal':
                        ncandidates = 1
                        xcandidates = num.zeros((ncandidates, npar))
                        for icandidate in range(ncandidates):
                            for ipar in range(npar):
                                ntry = 0
                                while True:
                                    if local_sxs[ibootstrap_choice][ipar] > 0.:
                                        v = num.random.normal(
                                            xchoice[ipar],
                                            factor*local_sxs[
                                                ibootstrap_choice][ipar])
                                    else:
                                        v = xchoice[ipar]

                                    if xbounds[ipar, 0] <= v and \
                                            v <= xbounds[ipar, 1]:

                                        break

                                    if ntry > 1000:
                                        raise GrondError(
                                            'failed to produce a suitable '
                                            'candidate sample from normal '
                                            'distribution')

                                    ntry += 1

                                xcandidates[icandidate, ipar] = v

                        x = select_most_excentric(
                            xcandidates,
                            xhist[chains_i[ibootstrap_choice, :], :],
                            local_sxs[ibootstrap_choice],
                            factor)

                try:
                    x = problem.preconstrain(x)
                    break

                except Forbidden:
                    pass

        ibase = None
        if ichoice is not None and ibootstrap_choice is not None:
            ibase = chains_i[ibootstrap_choice, ichoice]

        if isbad_mask is not None and num.any(isbad_mask):
            isok_mask = num.logical_not(isbad_mask)
        else:
            isok_mask = None

        ms, ns = problem.evaluate(x, mask=isok_mask)

        isbad_mask_new = num.isnan(ms)
        if isbad_mask is not None and num.any(isbad_mask != isbad_mask_new):
            logger.error(
                'skipping problem %s: inconsistency in data availability' %
                problem.name)

            for target, isbad_new, isbad in zip(
                    problem.targets, isbad_mask_new, isbad_mask):

                if isbad_new != isbad:
                    logger.error('%s, %s -> %s' % (
                        target.string_id(), isbad, isbad_new))

            return

        isbad_mask = isbad_mask_new

        if num.all(isbad_mask):
            logger.error(
                'skipping problem %s: all target misfit values are NaN' %
                problem.name)
            return

        gm = problem.global_misfit(ms, ns)
        bms = problem.bootstrap_misfit(ms, ns, nbootstrap)

        chains_m[0, nlinks] = gm
        chains_m[1:, nlinks] = bms
        chains_i[:, nlinks] = iiter

        nlinks += 1

        for ichain in range(chains_m.shape[0]):
            isort = num.argsort(chains_m[ichain, :nlinks])
            chains_m[ichain, :nlinks] = chains_m[ichain, isort]
            chains_i[ichain, :nlinks] = chains_i[ichain, isort]

        if nlinks == nlinks_cap:
            accept = (chains_i[:, nlinks_cap-1] != iiter).astype(num.int)
            nlinks -= 1
        else:
            accept = num.ones(1 + nbootstrap, dtype=num.int)

        if rundir:
            problem.dump_problem_data(
                rundir, x, ms, ns, accept,
                ibootstrap_choice if ibootstrap_choice is not None else -1,
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

            for i in range(1 + nbootstrap):
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
            print('\n'.join(lines))

        if plot and plot.want_to_update(iiter):
            plot.update(
                xhist[:iiter+1, :],
                chains_i[:, :nlinks],
                ibase,
                ibootstrap_choice,
                local_sxs,
                factor)

        iiter += 1

    if plot:
        plot.finish()


class HighScoreSolver(Solver):

    def __init__(self, kwargs):
        Solver.__init__(self)
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


class HighScoreSolverConfig(SolverConfig):
    nbootstrap = Int.T(default=100)
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

    def get_solver_kwargs(self):
        return dict(
            nbootstrap=self.nbootstrap,
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

    def get_solver(self):
        return HighScoreSolver(self.get_solver_kwargs())


__all__ = '''
    HighScoreSolver
    SamplerDistributionChoice
    StandardDeviationEstimatorChoice
    HighScoreSolverConfig
'''.split()
