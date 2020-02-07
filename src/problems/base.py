'''
Base classes for Grond's problem definition and the model history container.

Common behaviour of all source models offered by Grond is implemented here.
Source model specific details are implemented in the respective submodules.
'''

import numpy as num
import math
import copy
import logging
import os.path as op
import os
import time

from pyrocko import gf, util, guts
from pyrocko.guts import Object, String, List, Dict, Int

from grond.meta import ADict, Parameter, GrondError, xjoin, Forbidden, \
    StringID, has_get_plot_classes
from ..targets import MisfitResult, MisfitTarget, TargetGroup, \
    WaveformMisfitTarget, SatelliteMisfitTarget, GNSSCampaignMisfitTarget

from grond import stats

from grond.version import __version__

guts_prefix = 'grond'
logger = logging.getLogger('grond.problems.base')
km = 1e3
as_km = dict(scale_factor=km, scale_unit='km')

g_rstate = num.random.RandomState()


def nextpow2(i):
    return 2**int(math.ceil(math.log(i)/math.log(2.)))


def correlated_weights(values, weight_matrix):
    '''
    Applies correlated weights to values

    The resulting weighed values have to be squared! Check out
    :meth:`Problem.combine_misfits` for more information.

    :param values: Misfits or norms as :class:`numpy.Array`
    :param weight: Weight matrix, commonly the inverse of covariance matrix

    :returns: :class:`numpy.Array` weighted values
    '''
    return num.matmul(values, weight_matrix)


class ProblemConfig(Object):
    '''
    Base class for config section defining the objective function setup.

    Factory for :py:class:`Problem` objects.
    '''
    name_template = String.T()
    norm_exponent = Int.T(default=2)
    nthreads = Int.T(default=1)

    def get_problem(self, event, target_groups, targets):
        '''
        Instantiate the problem with a given event and targets.

        :returns: :py:class:`Problem` object
        '''
        raise NotImplementedError


@has_get_plot_classes
class Problem(Object):
    '''
    Base class for objective function setup.

    Defines the *problem* to be solved by the optimiser.
    '''
    name = String.T()
    ranges = Dict.T(String.T(), gf.Range.T())
    dependants = List.T(Parameter.T())
    norm_exponent = Int.T(default=2)
    base_source = gf.Source.T(optional=True)
    targets = List.T(MisfitTarget.T())
    target_groups = List.T(TargetGroup.T())
    grond_version = String.T(optional=True)
    nthreads = Int.T(default=1)

    def __init__(self, **kwargs):
        Object.__init__(self, **kwargs)

        if self.grond_version is None:
            self.grond_version = __version__

        self._target_weights = None
        self._engine = None
        self._family_mask = None

        if hasattr(self, 'problem_waveform_parameters') and self.has_waveforms:
            self.problem_parameters =\
                self.problem_parameters + self.problem_waveform_parameters

        unused_parameters = []
        for p in self.problem_parameters:
            if p.optional and p._name not in self.ranges.keys():
                unused_parameters.append(p)

        for p in unused_parameters:
            self.problem_parameters.remove(p)

        self.check()

    @classmethod
    def get_plot_classes(cls):
        from . import plot
        return plot.get_plot_classes()

    def check(self):
        paths = set()
        for grp in self.target_groups:
            if grp.path == 'all':
                continue
            if grp.path in paths:
                raise ValueError('Path %s defined more than once! In %s'
                                 % (grp.path, grp.__class__.__name__))
            paths.add(grp.path)
        logger.debug('TargetGroup check OK.')

    def copy(self):
        o = copy.copy(self)
        o._target_weights = None
        return o

    def set_target_parameter_values(self, x):
        nprob = len(self.problem_parameters)
        for target in self.targets:
            target.set_parameter_values(x[nprob:nprob+target.nparameters])
            nprob += target.nparameters

    def get_parameter_dict(self, model, group=None):
        params = []
        for ip, p in enumerate(self.parameters):
            if group in p.groups or group is None:
                params.append((p.name, model[ip]))
        return ADict(params)

    def get_parameter_array(self, d):
        arr = num.zeros(self.nparameters, dtype=num.float)
        for ip, p in enumerate(self.parameters):
            if p.name in d.keys():
                arr[ip] = d[p.name]
        return arr

    def dump_problem_info(self, dirname):
        fn = op.join(dirname, 'problem.yaml')
        util.ensuredirs(fn)
        guts.dump(self, filename=fn)

    def dump_problem_data(
            self, dirname, x, misfits, chains=None,
            sampler_context=None):

        fn = op.join(dirname, 'models')
        if not isinstance(x, num.ndarray):
            x = num.array(x)
        with open(fn, 'ab') as f:
            x.astype('<f8').tofile(f)

        fn = op.join(dirname, 'misfits')
        with open(fn, 'ab') as f:
            misfits.astype('<f8').tofile(f)

        if chains is not None:
            fn = op.join(dirname, 'chains')
            with open(fn, 'ab') as f:
                chains.astype('<f8').tofile(f)

        if sampler_context is not None:
            fn = op.join(dirname, 'choices')
            with open(fn, 'ab') as f:
                num.array(sampler_context, dtype='<i8').tofile(f)

    def name_to_index(self, name):
        pnames = [p.name for p in self.combined]
        return pnames.index(name)

    @property
    def parameters(self):
        target_parameters = []
        for target in self.targets:
            target_parameters.extend(target.target_parameters)
        return self.problem_parameters + target_parameters

    @property
    def parameter_names(self):
        return [p.name for p in self.combined]

    @property
    def dependant_names(self):
        return [p.name for p in self.dependants]

    @property
    def nparameters(self):
        return len(self.parameters)

    @property
    def ntargets(self):
        return len(self.targets)

    @property
    def nwaveform_targets(self):
        return len(self.waveform_targets)

    @property
    def nsatellite_targets(self):
        return len(self.satellite_targets)

    @property
    def ngnss_targets(self):
        return len(self.gnss_targets)

    @property
    def nmisfits(self):
        nmisfits = 0
        for target in self.targets:
            nmisfits += target.nmisfits
        return nmisfits

    @property
    def ndependants(self):
        return len(self.dependants)

    @property
    def ncombined(self):
        return len(self.parameters) + len(self.dependants)

    @property
    def combined(self):
        return self.parameters + self.dependants

    @property
    def satellite_targets(self):
        return [t for t in self.targets
                if isinstance(t, SatelliteMisfitTarget)]

    @property
    def gnss_targets(self):
        return [t for t in self.targets
                if isinstance(t, GNSSCampaignMisfitTarget)]

    @property
    def waveform_targets(self):
        return [t for t in self.targets
                if isinstance(t, WaveformMisfitTarget)]

    @property
    def has_satellite(self):
        if self.satellite_targets:
            return True
        return False

    @property
    def has_waveforms(self):
        if self.waveform_targets:
            return True
        return False

    def set_engine(self, engine):
        self._engine = engine

    def get_engine(self):
        return self._engine

    def get_gf_store(self, target):
        if self.get_engine() is None:
            raise GrondError('Cannot get GF Store, modelling is not set up.')
        return self.get_engine().get_store(target.store_id)

    def random_uniform(self, xbounds, rstate, fixed_magnitude=None):
        if fixed_magnitude is not None:
            raise GrondError(
                'Setting fixed magnitude in random model generation not '
                'supported for this type of problem.')

        x = rstate.uniform(0., 1., self.nparameters)
        x *= (xbounds[:, 1] - xbounds[:, 0])
        x += xbounds[:, 0]
        return x

    def preconstrain(self, x):
        return x

    def extract(self, xs, i):
        if xs.ndim == 1:
            return self.extract(xs[num.newaxis, :], i)[0]

        if i < self.nparameters:
            return xs[:, i]
        else:
            return self.make_dependant(
                xs, self.dependants[i-self.nparameters].name)

    def get_target_weights(self):
        if self._target_weights is None:
            self._target_weights = num.concatenate(
               [target.get_combined_weight() for target in self.targets])

        return self._target_weights

    def get_target_residuals(self):
        pass

    def inter_family_weights(self, ns):
        exp, root = self.get_norm_functions()

        family, nfamilies = self.get_family_mask()

        ws = num.zeros(self.nmisfits)
        for ifamily in range(nfamilies):
            mask = family == ifamily
            ws[mask] = 1.0 / root(num.nansum(exp(ns[mask])))

        return ws

    def inter_family_weights2(self, ns):
        '''
        :param ns: 2D array with normalization factors ``ns[imodel, itarget]``
        :returns: 2D array ``weights[imodel, itarget]``
        '''

        exp, root = self.get_norm_functions()
        family, nfamilies = self.get_family_mask()

        ws = num.zeros(ns.shape)
        for ifamily in range(nfamilies):
            mask = family == ifamily
            ws[:, mask] = (1.0 / root(
                num.nansum(exp(ns[:, mask]), axis=1)))[:, num.newaxis]

        return ws

    def get_reference_model(self):
        model = num.zeros(self.nparameters)
        model_source_params = self.pack(self.base_source)
        model[:model_source_params.size] = model_source_params
        return model

    def get_parameter_bounds(self):
        out = []
        for p in self.problem_parameters:
            r = self.ranges[p.name]
            out.append((r.start, r.stop))

        for target in self.targets:
            for p in target.target_parameters:
                r = target.target_ranges[p.name_nogroups]
                out.append((r.start, r.stop))

        return num.array(out, dtype=num.float)

    def get_dependant_bounds(self):
        return num.zeros((0, 2))

    def get_combined_bounds(self):
        return num.vstack((
            self.get_parameter_bounds(),
            self.get_dependant_bounds()))

    def raise_invalid_norm_exponent(self):
        raise GrondError('Invalid norm exponent: %f' % self.norm_exponent)

    def get_norm_functions(self):
        if self.norm_exponent == 2:
            def sqr(x):
                return x**2

            return sqr, num.sqrt

        elif self.norm_exponent == 1:
            def noop(x):
                return x

            return noop, num.abs

        else:
            self.raise_invalid_norm_exponent()

    def combine_misfits(
            self, misfits,
            extra_weights=None,
            extra_residuals=None,
            extra_correlated_weights=dict(),
            get_contributions=False):

        '''
        Combine misfit contributions (residuals) to global or bootstrap misfits

        :param misfits: 3D array ``misfits[imodel, iresidual, 0]`` are the
            misfit contributions (residuals) ``misfits[imodel, iresidual, 1]``
            are the normalisation contributions. It is also possible to give
            the misfit and normalisation contributions for a single model as
            ``misfits[iresidual, 0]`` and misfits[iresidual, 1]`` in which
            case, the first dimension (imodel) of the result will be stipped
            off.

        :param extra_weights: if given, 2D array of extra weights to be applied
            to the contributions, indexed as
            ``extra_weights[ibootstrap, iresidual]``.

        :param extra_residuals: if given, 2D array of perturbations to be added
            to the residuals, indexed as
            ``extra_residuals[ibootstrap, iresidual]``.

        :param extra_correlated_weights: if a dictionary of
            ``imisfit: correlated weight matrix`` is passed a correlated
            weight matrix is applied to the misfit and normalisation values.
            `imisfit` is the starting index in the misfits vector the
            correlated weight matrix applies to.

        :param get_contributions: get the weighted and perturbed contributions
            (don't do the sum).

        :returns: if no *extra_weights* or *extra_residuals* are given, a 1D
            array indexed as ``misfits[imodel]`` containing the global misfit
            for each model is returned, otherwise a 2D array
            ``misfits[imodel, ibootstrap]`` with the misfit for every model and
            weighting/residual set is returned.
        '''
        if misfits.ndim == 2:
            misfits = misfits[num.newaxis, :, :]
            return self.combine_misfits(
                misfits, extra_weights, extra_residuals,
                extra_correlated_weights, get_contributions)[0, ...]

        if extra_weights is None and extra_residuals is None:
            return self.combine_misfits(
                misfits, False, False,
                extra_correlated_weights, get_contributions)[:, 0]

        assert misfits.ndim == 3
        assert not num.any(extra_weights) or extra_weights.ndim == 2
        assert not num.any(extra_residuals) or extra_residuals.ndim == 2

        if self.norm_exponent != 2 and extra_correlated_weights:
            raise GrondError('Correlated weights can only be used '
                             ' with norm_exponent=2')

        exp, root = self.get_norm_functions()

        nmodels = misfits.shape[0]
        nmisfits = misfits.shape[1]  # noqa

        mf = misfits[:, num.newaxis, :, :].copy()

        if num.any(extra_residuals):
            mf = mf + extra_residuals[num.newaxis, :, :, num.newaxis]

        res = mf[..., 0]
        norms = mf[..., 1]

        for imisfit, corr_weight_mat in extra_correlated_weights.items():

            jmisfit = imisfit + corr_weight_mat.shape[0]

            for imodel in range(nmodels):
                corr_res = res[imodel, :, imisfit:jmisfit]
                corr_norms = norms[imodel, :, imisfit:jmisfit]

                res[imodel, :, imisfit:jmisfit] = \
                    correlated_weights(corr_res, corr_weight_mat)

                norms[imodel, :, imisfit:jmisfit] = \
                    correlated_weights(corr_norms, corr_weight_mat)

        # Apply normalization family weights (these weights depend on
        # on just calculated correlated norms!)
        weights_fam = \
            self.inter_family_weights2(norms[:, 0, :])[:, num.newaxis, :]

        weights_fam = exp(weights_fam)

        res = exp(res)
        norms = exp(norms)

        res *= weights_fam
        norms *= weights_fam

        weights_tar = self.get_target_weights()[num.newaxis, num.newaxis, :]
        if num.any(extra_weights):
            weights_tar = weights_tar * extra_weights[num.newaxis, :, :]

        weights_tar = exp(weights_tar)

        res = res * weights_tar
        norms = norms * weights_tar

        if get_contributions:
            return res / num.nansum(norms, axis=2)[:, :, num.newaxis]

        result = root(
            num.nansum(res, axis=2) /
            num.nansum(norms, axis=2))

        assert result[result < 0].size == 0
        return result

    def make_family_mask(self):
        family_names = set()
        families = num.zeros(self.nmisfits, dtype=num.int)

        idx = 0
        for itarget, target in enumerate(self.targets):
            family_names.add(target.normalisation_family)
            families[idx:idx + target.nmisfits] = len(family_names) - 1
            idx += target.nmisfits

        return families, len(family_names)

    def get_family_mask(self):
        if self._family_mask is None:
            self._family_mask = self.make_family_mask()

        return self._family_mask

    def evaluate(self, x, mask=None, result_mode='full', targets=None,
                 nthreads=1):
        source = self.get_source(x)
        engine = self.get_engine()

        self.set_target_parameter_values(x)

        if mask is not None and targets is not None:
            raise ValueError('Mask cannot be defined with targets set.')
        targets = targets if targets is not None else self.targets

        for target in targets:
            target.set_result_mode(result_mode)

        modelling_targets = []
        t2m_map = {}
        for itarget, target in enumerate(targets):
            t2m_map[target] = target.prepare_modelling(engine, source, targets)
            if mask is None or mask[itarget]:
                modelling_targets.extend(t2m_map[target])

        u2m_map = {}
        for imtarget, mtarget in enumerate(modelling_targets):
            if mtarget not in u2m_map:
                u2m_map[mtarget] = []

            u2m_map[mtarget].append(imtarget)

        modelling_targets_unique = list(u2m_map.keys())

        resp = engine.process(source, modelling_targets_unique,
                              nthreads=nthreads)
        modelling_results_unique = list(resp.results_list[0])

        modelling_results = [None] * len(modelling_targets)

        for mtarget, mresult in zip(
                modelling_targets_unique, modelling_results_unique):

            for itarget in u2m_map[mtarget]:
                modelling_results[itarget] = mresult

        imt = 0
        results = []
        for itarget, target in enumerate(targets):
            nmt_this = len(t2m_map[target])
            if mask is None or mask[itarget]:
                result = target.finalize_modelling(
                    engine, source,
                    t2m_map[target],
                    modelling_results[imt:imt+nmt_this])

                imt += nmt_this
            else:
                result = gf.SeismosizerError(
                    'target was excluded from modelling')

            results.append(result)

        return results

    def misfits(self, x, mask=None, nthreads=1):
        results = self.evaluate(
            x, mask=mask, result_mode='sparse', nthreads=nthreads)
        misfits = num.full((self.nmisfits, 2), num.nan)

        imisfit = 0
        for target, result in zip(self.targets, results):
            if isinstance(result, MisfitResult):
                misfits[imisfit:imisfit+target.nmisfits, :] = result.misfits

            imisfit += target.nmisfits

        return misfits

    def forward(self, x):
        source = self.get_source(x)
        engine = self.get_engine()

        plain_targets = []
        for target in self.targets:
            plain_targets.extend(target.get_plain_targets(engine, source))

        resp = engine.process(source, plain_targets)

        results = []
        for target, result in zip(plain_targets, resp.results_list[0]):
            if isinstance(result, gf.SeismosizerError):
                logger.debug(
                    '%s.%s.%s.%s: %s' % (target.codes + (str(result),)))
            else:
                results.append(result)

        return results

    def get_random_model(self, ntries_limit=100):
        xbounds = self.get_parameter_bounds()

        for _ in range(ntries_limit):
            x = self.random_uniform(xbounds, rstate=g_rstate)
            try:
                return self.preconstrain(x)

            except Forbidden:
                pass

        raise GrondError(
            'Could not find any suitable candidate sample within %i tries' % (
                ntries_limit))


class ProblemInfoNotAvailable(GrondError):
    pass


class ProblemDataNotAvailable(GrondError):
    pass


class NoSuchAttribute(GrondError):
    pass


class InvalidAttributeName(GrondError):
    pass


class ModelHistory(object):
    '''
    Write, read and follow sequences of models produced in an optimisation run.

    :param problem: :class:`grond.Problem` instance
    :param path: path to rundir, defaults to None
    :type path: str, optional
    :param mode: open mode, 'r': read, 'w': write
    :type mode: str, optional
    '''

    nmodels_capacity_min = 1024

    def __init__(self, problem, nchains=None, path=None, mode='r'):
        self.mode = mode

        self.problem = problem
        self.path = path
        self.nchains = nchains

        self._models_buffer = None
        self._misfits_buffer = None
        self._bootstraps_buffer = None
        self._sample_contexts_buffer = None

        self._sorted_misfit_idx = {}

        self.models = None
        self.misfits = None
        self.bootstrap_misfits = None

        self.sampler_contexts = None

        self.nmodels_capacity = self.nmodels_capacity_min
        self.listeners = []

        self._attributes = {}

        if mode == 'r':
            self.load()

    @staticmethod
    def verify_rundir(rundir):
        _rundir_files = ('misfits', 'models')

        if not op.exists(rundir):
            raise ProblemDataNotAvailable(
                'Directory does not exist: %s' % rundir)
        for f in _rundir_files:
            if not op.exists(op.join(rundir, f)):
                raise ProblemDataNotAvailable('File not found: %s' % f)

    @classmethod
    def follow(cls, path, nchains=None, wait=20.):
        '''
        Start following a rundir (constructor).

        :param path: the path to follow, a grond rundir
        :type path: str, optional
        :param wait: wait time until the folder become alive
        :type wait: number in seconds, optional
        :returns: A :py:class:`ModelHistory` instance
        '''
        start_watch = time.time()
        while (time.time() - start_watch) < wait:
            try:
                cls.verify_rundir(path)
                problem = load_problem_info(path)
                return cls(problem, nchains=nchains, path=path, mode='r')
            except (ProblemDataNotAvailable, OSError):
                time.sleep(.25)

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
        if self.nchains is not None:
            self.bootstrap_misfits = self._bootstraps_buffer[:nmodels_new, :, :]  # noqa
        if self._sample_contexts_buffer is not None:
            self.sampler_contexts = self._sample_contexts_buffer[:nmodels_new, :]  # noqa

    @property
    def nmodels_capacity(self):
        if self._models_buffer is None:
            return 0
        else:
            return self._models_buffer.shape[0]

    @nmodels_capacity.setter
    def nmodels_capacity(self, nmodels_capacity_new):
        if self.nmodels_capacity != nmodels_capacity_new:

            models_buffer = num.zeros(
                (nmodels_capacity_new, self.problem.nparameters),
                dtype=num.float)
            misfits_buffer = num.zeros(
                (nmodels_capacity_new, self.problem.nmisfits, 2),
                dtype=num.float)
            sample_contexts_buffer = num.zeros(
                (nmodels_capacity_new, 4),
                dtype=num.int)
            sample_contexts_buffer.fill(-1)

            if self.nchains is not None:
                bootstraps_buffer = num.zeros(
                    (nmodels_capacity_new, self.nchains),
                    dtype=num.float)

            ncopy = min(self.nmodels, nmodels_capacity_new)

            if self._models_buffer is not None:
                models_buffer[:ncopy, :] = \
                    self._models_buffer[:ncopy, :]
                misfits_buffer[:ncopy, :, :] = \
                    self._misfits_buffer[:ncopy, :, :]
                sample_contexts_buffer[:ncopy, :] = \
                    self._sample_contexts_buffer[:ncopy, :]

            self._models_buffer = models_buffer
            self._misfits_buffer = misfits_buffer
            self._sample_contexts_buffer = sample_contexts_buffer

            if self.nchains is not None:
                if self._bootstraps_buffer is not None:
                    bootstraps_buffer[:ncopy, :] = \
                        self._bootstraps_buffer[:ncopy, :]
                self._bootstraps_buffer = bootstraps_buffer

    def clear(self):
        assert self.mode != 'r', 'History is read-only, cannot clear.'
        self.nmodels = 0
        self.nmodels_capacity = self.nmodels_capacity_min

    def extend(
            self, models, misfits,
            bootstrap_misfits=None,
            sampler_contexts=None):

        nmodels = self.nmodels
        n = models.shape[0]

        nmodels_capacity_want = max(
            self.nmodels_capacity_min, nextpow2(nmodels + n))

        if nmodels_capacity_want != self.nmodels_capacity:
            self.nmodels_capacity = nmodels_capacity_want

        self._models_buffer[nmodels:nmodels+n, :] = models
        self._misfits_buffer[nmodels:nmodels+n, :, :] = misfits

        self.models = self._models_buffer[:nmodels+n, :]
        self.misfits = self._misfits_buffer[:nmodels+n, :, :]

        if bootstrap_misfits is not None:
            self._bootstraps_buffer[nmodels:nmodels+n, :] = bootstrap_misfits
            self.bootstrap_misfits = self._bootstraps_buffer[:nmodels+n, :]

        if sampler_contexts is not None:
            self._sample_contexts_buffer[nmodels:nmodels+n, :] \
                = sampler_contexts
            self.sampler_contexts = self._sample_contexts_buffer[:nmodels+n, :]

        if self.path and self.mode == 'w':
            for i in range(n):
                self.problem.dump_problem_data(
                    self.path, models[i, :], misfits[i, :, :],
                    bootstrap_misfits[i, :]
                    if bootstrap_misfits is not None else None,
                    sampler_contexts[i, :]
                    if sampler_contexts is not None else None)

        self._sorted_misfit_idx.clear()

        self.emit('extend', nmodels, n, models, misfits, sampler_contexts)

    def append(
            self, model, misfits,
            bootstrap_misfits=None,
            sampler_context=None):

        if bootstrap_misfits is not None:
            bootstrap_misfits = bootstrap_misfits[num.newaxis, :]

        if sampler_context is not None:
            sampler_context = sampler_context[num.newaxis, :]

        return self.extend(
            model[num.newaxis, :], misfits[num.newaxis, :, :],
            bootstrap_misfits, sampler_context)

    def load(self):
        self.mode = 'r'
        self.verify_rundir(self.path)
        models, misfits, bootstraps, sampler_contexts = load_problem_data(
            self.path, self.problem, nchains=self.nchains)
        self.extend(models, misfits, bootstraps, sampler_contexts)

    def update(self):
        ''' Update history from path '''
        nmodels_available = get_nmodels(self.path, self.problem)
        if self.nmodels == nmodels_available:
            return

        try:
            new_models, new_misfits, new_bootstraps, new_sampler_contexts = \
                load_problem_data(
                    self.path,
                    self.problem,
                    nmodels_skip=self.nmodels,
                    nchains=self.nchains)

        except ValueError:
            return

        self.extend(
            new_models,
            new_misfits,
            new_bootstraps,
            new_sampler_contexts)

    def add_listener(self, listener):
        ''' Add a listener to the history

        The listening class can implement the following methods:
        * ``extend``
        '''
        self.listeners.append(listener)

    def emit(self, event_name, *args, **kwargs):
        for listener in self.listeners:
            slot = getattr(listener, event_name, None)
            if callable(slot):
                slot(*args, **kwargs)

    @property
    def attribute_names(self):
        apath = op.join(self.path, 'attributes')
        if not os.path.exists(apath):
            return []

        return [fn for fn in os.listdir(apath)
                if StringID.regex.match(fn)]

    def get_attribute(self, name):
        if name not in self._attributes:
            if name not in self.attribute_names:
                raise NoSuchAttribute(name)

            path = op.join(self.path, 'attributes', name)

            with open(path, 'rb') as f:
                self._attributes[name] = num.fromfile(
                    f, dtype='<i4',
                    count=self.nmodels).astype(num.int)

            assert self._attributes[name].shape == (self.nmodels,)

        return self._attributes[name]

    def set_attribute(self, name, attribute):
        if not StringID.regex.match(name):
            raise InvalidAttributeName(name)

        attribute = attribute.astype(num.int)
        assert attribute.shape == (self.nmodels,)

        apath = op.join(self.path, 'attributes')

        if not os.path.exists(apath):
            os.mkdir(apath)

        path = op.join(apath, name)

        with open(path, 'wb') as f:
            attribute.astype('<i4').tofile(f)

        self._attributes[name] = attribute

    def ensure_bootstrap_misfits(self, optimiser):
        if self.bootstrap_misfits is None:
            problem = self.problem
            self.bootstrap_misfits = problem.combine_misfits(
                self.misfits,
                extra_weights=optimiser.get_bootstrap_weights(problem),
                extra_residuals=optimiser.get_bootstrap_residuals(problem))

    def imodels_by_cluster(self, cluster_attribute):
        if cluster_attribute is None:
            return [(-1, 100.0, num.arange(self.nmodels))]

        by_cluster = []
        try:
            iclusters = self.get_attribute(cluster_attribute)
            iclusters_avail = num.unique(iclusters)

            for icluster in iclusters_avail:
                imodels = num.where(iclusters == icluster)[0]
                by_cluster.append(
                    (icluster,
                     (100.0 * imodels.size) / self.nmodels,
                     imodels))

            if by_cluster and by_cluster[0][0] == -1:
                by_cluster.append(by_cluster.pop(0))

        except NoSuchAttribute:
            logger.warn(
                'Attribute %s not set in run %s.\n'
                '  Skipping model retrieval by clusters.' % (
                    cluster_attribute, self.problem.name))

        return by_cluster

    def models_by_cluster(self, cluster_attribute):
        if cluster_attribute is None:
            return [(-1, 100.0, self.models)]

        return [
            (icluster, percentage, self.models[imodels])
            for (icluster, percentage, imodels)
            in self.imodels_by_cluster(cluster_attribute)]

    def mean_sources_by_cluster(self, cluster_attribute):
        return [
            (icluster, percentage, stats.get_mean_source(self.problem, models))
            for (icluster, percentage, models)
            in self.models_by_cluster(cluster_attribute)]

    def get_sorted_misfits_idx(self, chain=0):
        if chain not in self._sorted_misfit_idx.keys():
            self._sorted_misfit_idx[chain] = num.argsort(
                self.bootstrap_misfits[:, chain])

        return self._sorted_misfit_idx[chain]

    def get_sorted_misfits(self, chain=0):
        isort = self.get_sorted_misfits_idx(chain)
        return self.bootstrap_misfits[:, chain][isort]

    def get_sorted_models(self, chain=0):
        isort = self.get_sorted_misfits_idx(chain=0)
        return self.models[isort, :]

    def get_sorted_primary_misfits(self):
        return self.get_sorted_misfits(chain=0)

    def get_sorted_primary_models(self):
        return self.get_sorted_models(chain=0)

    def get_best_model(self, chain=0):
        return self.get_sorted_models(chain)[0, ...]

    def get_best_misfit(self, chain=0):
        return self.get_sorted_misfits(chain)[0]

    def get_mean_model(self):
        return num.mean(self.models, axis=0)

    def get_mean_misfit(self, chain=0):
        return num.mean(self.bootstrap_misfits[:, chain])

    def get_best_source(self, chain=0):
        return self.problem.get_source(self.get_best_model(chain))

    def get_mean_source(self, chain=0):
        return self.problem.get_source(self.get_mean_model())

    def get_chain_misfits(self, chain=0):
        return self.bootstrap_misfits[:, chain]

    def get_primary_chain_misfits(self):
        return self.get_chain_misfits(chain=0)


def get_nmodels(dirname, problem):
    fn = op.join(dirname, 'models')
    with open(fn, 'r') as f:
        nmodels1 = os.fstat(f.fileno()).st_size // (problem.nparameters * 8)

    fn = op.join(dirname, 'misfits')
    with open(fn, 'r') as f:
        nmodels2 = os.fstat(f.fileno()).st_size // (problem.nmisfits * 2 * 8)

    return min(nmodels1, nmodels2)


def load_problem_info_and_data(dirname, subset=None, nchains=None):
    problem = load_problem_info(dirname)
    models, misfits, bootstraps, sampler_contexts = load_problem_data(
        xjoin(dirname, subset), problem, nchains=nchains)
    return problem, models, misfits, bootstraps, sampler_contexts


def load_optimiser_info(dirname):
    fn = op.join(dirname, 'optimiser.yaml')
    return guts.load(filename=fn)


def load_problem_info(dirname):
    try:
        fn = op.join(dirname, 'problem.yaml')
        return guts.load(filename=fn)
    except OSError as e:
        logger.debug(e)
        raise ProblemInfoNotAvailable(
            'No problem info available (%s).' % dirname)


def load_problem_data(dirname, problem, nmodels_skip=0, nchains=None):

    def get_chains_fn():
        for fn in (op.join(dirname, 'bootstraps'),
                   op.join(dirname, 'chains')):
            if op.exists(fn):
                return fn
        return False

    try:
        nmodels = get_nmodels(dirname, problem) - nmodels_skip

        fn = op.join(dirname, 'models')
        with open(fn, 'r') as f:
            f.seek(nmodels_skip * problem.nparameters * 8)
            models = num.fromfile(
                    f, dtype='<f8',
                    count=nmodels * problem.nparameters)\
                .astype(num.float)

        models = models.reshape((nmodels, problem.nparameters))

        fn = op.join(dirname, 'misfits')
        with open(fn, 'r') as f:
            f.seek(nmodels_skip * problem.nmisfits * 2 * 8)
            misfits = num.fromfile(
                    f, dtype='<f8',
                    count=nmodels*problem.nmisfits*2)\
                .astype(num.float)
        misfits = misfits.reshape((nmodels, problem.nmisfits, 2))

        chains = None
        fn = get_chains_fn()
        if fn and nchains is not None:
            with open(fn, 'r') as f:
                f.seek(nmodels_skip * nchains * 8)
                chains = num.fromfile(
                        f, dtype='<f8',
                        count=nmodels*nchains)\
                    .astype(num.float)

            chains = chains.reshape((nmodels, nchains))

        sampler_contexts = None
        fn = op.join(dirname, 'choices')
        if op.exists(fn):
            with open(fn, 'r') as f:
                f.seek(nmodels_skip * 4 * 8)
                sampler_contexts = num.fromfile(
                        f, dtype='<i8',
                        count=nmodels*4).astype(num.int)

            sampler_contexts = sampler_contexts.reshape((nmodels, 4))

    except OSError as e:
        logger.debug(str(e))
        raise ProblemDataNotAvailable(
            'No problem data available (%s).' % dirname)

    return models, misfits, chains, sampler_contexts


__all__ = '''
    ProblemConfig
    Problem
    ModelHistory
    ProblemInfoNotAvailable
    ProblemDataNotAvailable
    load_problem_info
    load_problem_info_and_data
    InvalidAttributeName
    NoSuchAttribute
'''.split()
