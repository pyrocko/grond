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

from ..meta import ADict, Parameter, GrondError, xjoin, Forbidden
from ..targets import MisfitResult, MisfitTarget, TargetGroup, \
    WaveformMisfitTarget, SatelliteMisfitTarget, GNSSCampaignMisfitTarget


guts_prefix = 'grond'
logger = logging.getLogger('grond.problems.base')
km = 1e3
as_km = dict(scale_factor=km, scale_unit='km')


def nextpow2(i):
    return 2**int(math.ceil(math.log(i)/math.log(2.)))


class ProblemConfig(Object):
    '''
    Base class for config section defining the objective function setup.

    Factory for :py:class:`Problem` objects.
    '''
    name_template = String.T()
    norm_exponent = Int.T(default=2)

    def get_problem(self, event, target_groups, targets):
        '''
        Instantiate the problem with a given event and targets.

        :returns: :py:class:`Problem` object
        '''
        raise NotImplementedError


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

    def __init__(self, **kwargs):
        Object.__init__(self, **kwargs)

        self._target_weights = None
        self._engine = None
        self._family_mask = None

        if hasattr(self, 'problem_waveform_parameters') and self.has_waveforms:
            self.problem_parameters =\
                self.problem_parameters + self.problem_waveform_parameters

        logger.name = self.__class__.__name__

    @classmethod
    def get_plot_classes(cls):
        from . import plot
        return plot.get_plot_classes()

    def get_engine(self):
        return self._engine

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
            self, dirname, x, misfits,
            accept=None, ibootstrap_choice=None, ibase=None):

        fn = op.join(dirname, 'models')
        if not isinstance(x, num.ndarray):
            x = num.array(x)
        with open(fn, 'ab') as f:
            x.astype('<f8').tofile(f)

        fn = op.join(dirname, 'misfits')
        with open(fn, 'ab') as f:
            misfits.astype('<f8').tofile(f)

        if None not in (ibootstrap_choice, ibase):
            fn = op.join(dirname, 'choices')
            with open(fn, 'ab') as f:
                num.array((ibootstrap_choice, ibase), dtype='<i8').tofile(f)

        if accept is not None:
            fn = op.join(dirname, 'accepted')
            with open(fn, 'ab') as f:
                accept.astype('<i1').tofile(f)

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

    def random_uniform(self, xbounds):
        x = num.random.uniform(0., 1., self.nparameters)
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
                [target.get_combined_weight()
                 for target in self.targets])

        return self._target_weights

    def inter_family_weights(self, ns):
        exp, root = self.get_norm_functions()

        family, nfamilies = self.get_family_mask()

        ws = num.zeros(self.ntargets)
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
        return self.pack(self.base_source)

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
        raise GrondError('invalid norm exponent' % self.norm_exponent)

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

        :param get_contributions: get the weighted and perturbed contributions
            (don't do the sum).

        :returns: if no *extra_weights* or *extra_residuals* are given, a 1D
            array indexed as ``misfits[imodel]`` containing the global misfit
            for each model is returned, otherwise a 2D array
            ``misfits[imodel, ibootstrap]`` with the misfit for every model and
            weighting/residual set is returned.
        '''

        exp, root = self.get_norm_functions()

        if misfits.ndim == 2:
            misfits = misfits[num.newaxis, :, :]
            return self.combine_misfits(
                misfits,
                extra_weights,
                extra_residuals,
                get_contributions)[0, ...]

        assert misfits.ndim == 3
        assert extra_weights is None or extra_weights.ndim == 2
        assert extra_residuals is None or extra_residuals.ndim == 2

        if extra_weights is not None or extra_residuals is not None:
            if extra_weights is not None:
                w = extra_weights[num.newaxis, :, :] \
                    * self.get_target_weights()[num.newaxis, num.newaxis, :] \
                    * self.inter_family_weights2(
                        misfits[:, :, 1])[:, num.newaxis, :]
            else:
                w = 1.0

            if extra_residuals is not None:
                r = extra_residuals[num.newaxis, :, :]
            else:
                r = 0.0

            if get_contributions:
                return exp(w*(misfits[:, num.newaxis, :, 0]+r)) \
                    / num.nansum(
                        exp(w*misfits[:, num.newaxis, :, 1]),
                        axis=2)[:, :, num.newaxis]

            return root(
                num.nansum(exp(w*(misfits[:, num.newaxis, :, 0]+r)), axis=2) /
                num.nansum(exp(w*(misfits[:, num.newaxis, :, 1])), axis=2))
        else:
            w = self.get_target_weights()[num.newaxis, :] \
                * self.inter_family_weights2(misfits[:, :, 1])

            if get_contributions:
                return exp(w*misfits[:, :, 0]) \
                    / num.nansum(
                        exp(w*misfits[:, :, 1]),
                        axis=1)[:, num.newaxis]

            return root(
                num.nansum(exp(w*misfits[:, :, 0]), axis=1) /
                num.nansum(exp(w*misfits[:, :, 1]), axis=1))

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

    def evaluate(self, x, mask=None, result_mode='full', targets=None):
        source = self.get_source(x)
        engine = self.get_engine()

        self.set_target_parameter_values(x)

        if mask is not None and targets is not None:
            raise ValueError('mask cannot be defined with targets set')
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

        resp = engine.process(source, modelling_targets_unique)
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

    def misfits(self, x, mask=None):
        results = self.evaluate(x, mask=mask, result_mode='sparse')
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

    def get_random_model(self):
        xbounds = self.get_parameter_bounds()

        while True:
            x = self.random_uniform(xbounds)
            try:
                return self.preconstrain(x)

            except Forbidden:
                pass


class ProblemInfoNotAvailable(GrondError):
    pass


class ProblemDataNotAvailable(GrondError):
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

    def __init__(self, problem, path=None, mode='r'):

        self.problem = problem
        self.path = path
        self._models_buffer = None
        self._misfits_buffer = None
        self.models = None
        self.misfits = None
        self.nmodels_capacity = self.nmodels_capacity_min
        self.listeners = []
        self.mode = mode

        if mode == 'r':
            self.verify_rundir(self.path)
            models, misfits = load_problem_data(path, problem)
            self.extend(models, misfits)

    @staticmethod
    def verify_rundir(rundir):
        _rundir_files = ['misfits', 'models']

        if not op.exists(rundir):
            raise ProblemDataNotAvailable(
                'Directory %s does not exist!' % rundir)
        for f in _rundir_files:
            if not op.exists(op.join(rundir, f)):
                raise ProblemDataNotAvailable('File %s not found!' % f)

    @classmethod
    def follow(cls, path, wait=20.):
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
                return cls(problem, path, mode='r')
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

            ncopy = min(self.nmodels, nmodels_capacity_new)

            if self._models_buffer is not None:
                models_buffer[:ncopy, :] = \
                    self._models_buffer[:ncopy, :]
                misfits_buffer[:ncopy, :, :] = \
                    self._misfits_buffer[:ncopy, :, :]

            self._models_buffer = models_buffer
            self._misfits_buffer = misfits_buffer

    def clear(self):
        self.nmodels = 0
        self.nmodels_capacity = self.nmodels_capacity_min

    def extend(self, models, misfits):
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

        if self.path and self.mode == 'w':
            for i in range(n):
                self.problem.dump_problem_data(
                        self.path, models[i, :], misfits[i, :, :])

        self.emit('extend', nmodels, n, models, misfits)

    def append(self, model, misfits):
        nmodels = self.nmodels

        nmodels_capacity_want = max(
            self.nmodels_capacity_min, nextpow2(nmodels + 1))

        if nmodels_capacity_want != self.nmodels_capacity:
            self.nmodels_capacity = nmodels_capacity_want

        self._models_buffer[nmodels, :] = model
        self._misfits_buffer[nmodels, :, :] = misfits

        self.models = self._models_buffer[:nmodels+1, :]
        self.misfits = self._misfits_buffer[:nmodels+1, :, :]

        if self.path and self.mode == 'w':
            self.problem.dump_problem_data(
                self.path, model, misfits)

        self.emit(
            'extend', nmodels, 1,
            model[num.newaxis, :], misfits[num.newaxis, :, :])

    def update(self):
        ''' Update history from path '''
        nmodels_available = get_nmodels(self.path, self.problem)
        if self.nmodels == nmodels_available:
            return
        new_models, new_misfits = load_problem_data(
            self.path, self.problem, nmodels_skip=self.nmodels)

        self.extend(new_models, new_misfits)

    def add_listener(self, listener):
        self.listeners.append(listener)

    def emit(self, event_name, *args, **kwargs):
        for listener in self.listeners:
            slot = getattr(listener, event_name, None)
            if callable(slot):
                slot(*args, **kwargs)


def get_nmodels(dirname, problem):
    fn = op.join(dirname, 'models')
    with open(fn, 'r') as f:
        nmodels1 = os.fstat(f.fileno()).st_size // (problem.nparameters * 8)

    fn = op.join(dirname, 'misfits')
    with open(fn, 'r') as f:
        nmodels2 = os.fstat(f.fileno()).st_size // (problem.nmisfits * 2 * 8)

    return min(nmodels1, nmodels2)


def load_problem_info_and_data(dirname, subset=None):
    problem = load_problem_info(dirname)
    models, misfits = load_problem_data(xjoin(dirname, subset), problem)
    return problem, models, misfits


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
            'no problem info available (%s)' % dirname)


def load_problem_data(dirname, problem, nmodels_skip=0):

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

    except OSError as e:
        logger.debug(str(e))
        raise ProblemDataNotAvailable(
            'no problem data available (%s)' % dirname)

    return models, misfits


__all__ = '''
    ProblemConfig
    Problem
    ModelHistory
    ProblemInfoNotAvailable
    ProblemDataNotAvailable
    load_problem_info
    load_problem_info_and_data
'''.split()
