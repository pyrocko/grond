import numpy as num
import math
import copy
import logging
import os.path as op
import os

from pyrocko import gf, util, guts
from pyrocko.guts import Object, String, Bool, List, Dict, Int

from ..meta import ADict, Parameter, GrondError, xjoin
from ..targets import WaveformMisfitTarget, SatelliteMisfitTarget


guts_prefix = 'grond'
logger = logging.getLogger('grond')
km = 1e3
as_km = dict(scale_factor=km, scale_unit='km')


def nextpow2(i):
    return 2**int(math.ceil(math.log(i)/math.log(2.)))


class ProblemConfig(Object):
    name_template = String.T()
    apply_balancing_weights = Bool.T(default=True)
    norm_exponent = Int.T(default=2)


class Problem(Object):
    name = String.T()
    ranges = Dict.T(String.T(), gf.Range.T())
    dependants = List.T(Parameter.T())
    apply_balancing_weights = Bool.T(default=True)
    norm_exponent = Int.T(default=2)
    base_source = gf.Source.T(optional=True)

    targets = List.T()

    def __init__(self, **kwargs):
        Object.__init__(self, **kwargs)

        self._bootstrap_weights = None
        self._target_weights = None
        self._engine = None
        self._group_mask = None

        if hasattr(self, 'problem_waveform_parameters') and self.has_waveforms:
            self.problem_parameters =\
                self.problem_parameters + self.problem_waveform_parameters

        logger.name = self.__class__.__name__

    def get_engine(self):
        return self._engine

    def copy(self):
        o = copy.copy(self)
        o._bootstrap_weights = None
        o._target_weights = None
        return o

    def set_target_parameter_values(self, x):
        i = len(self.problem_parameters)
        for target in self.targets:
            n = len(target.target_parameters)
            target.set_parameter_values(x[i:i+n])
            i += n

    def get_parameter_dict(self, model, group=None):
        params = []
        for ip, p in enumerate(self.parameters):
            if group in p.groups:
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
    def nparameters(self):
        return len(self.parameters)

    @property
    def ntargets(self):
        return len(self.targets)

    @property
    def ntargets_waveform(self):
        return len(self.waveform_targets)

    @property
    def ntargets_static(self):
        return len(self.satellite_targets)

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
    def waveform_targets(self):
        return [t for t in self.targets
                if isinstance(t, WaveformMisfitTarget)]

    @property
    def has_statics(self):
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
        raise NotImplementedError()

    def extract(self, xs, i):
        if xs.ndim == 1:
            return self.extract(xs[num.newaxis, :], i)[0]

        if i < self.nparameters:
            return xs[:, i]
        else:
            return self.make_dependant(
                xs, self.dependants[i-self.nparameters].name)

    def make_bootstrap_weights(self, nbootstrap, type='classic'):
        ntargets = self.ntargets
        ws = num.zeros((nbootstrap, ntargets))
        rstate = num.random.RandomState(23)
        for ibootstrap in range(nbootstrap):
            if type == 'classic':
                ii = rstate.randint(0, ntargets, size=self.ntargets)
                ws[ibootstrap, :] = num.histogram(
                    ii, ntargets, (-0.5, ntargets - 0.5))[0]
            elif type == 'bayesian':
                f = rstate.uniform(0., 1., size=self.ntargets+1)
                f[0] = 0.
                f[-1] = 1.
                f = num.sort(f)
                g = f[1:] - f[:-1]
                ws[ibootstrap, :] = g * ntargets
            else:
                assert False

        return ws

    def get_bootstrap_weights(self, nbootstrap, ibootstrap=None):
        if self._bootstrap_weights is None:
            self._bootstrap_weights = self.make_bootstrap_weights(
                nbootstrap, type='bayesian')

        if ibootstrap is None:
            return self._bootstrap_weights
        else:
            return self._bootstrap_weights[ibootstrap, :]

    def get_target_weights(self):
        if self._target_weights is None:
            self._target_weights = num.array(
                [target.get_combined_weight(
                    apply_balancing_weights=self.apply_balancing_weights)
                 for target in self.targets], dtype=num.float)

        return self._target_weights

    def inter_group_weights(self, ns):
        exp, root = self.get_norm_functions()

        group, ngroups = self.get_group_mask()

        ws = num.zeros(self.ntargets)
        for igroup in range(ngroups):
            mask = group == igroup
            ws[mask] = 1.0 / root(num.nansum(exp(ns[mask])))

        return ws

    def inter_group_weights2(self, ns):
        exp, root = self.get_norm_functions()

        group, ngroups = self.get_group_mask()
        ws = num.zeros(ns.shape)
        for igroup in range(ngroups):
            mask = group == igroup
            ws[:, mask] = (1.0 / root(
                num.nansum(exp(ns[:, mask]), axis=1)))[:, num.newaxis]

        return ws

    def get_xref(self):
        return self.pack(self.base_source)

    def get_parameter_bounds(self):
        out = []
        for p in self.problem_parameters:
            r = self.ranges[p.name]
            out.append((r.start, r.stop))

        for target in self.targets:
            for p in target.target_parameters:
                r = target.target_ranges[p.name]
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

    def bootstrap_misfit(self, misfits, nbootstrap, ibootstrap=None):
        exp, root = self.get_norm_functions()

        ms = misfits[:, 0]
        ns = misfits[:, 1]

        w = self.get_bootstrap_weights(nbootstrap, ibootstrap) * \
            self.get_target_weights() * self.inter_group_weights(ns)

        if ibootstrap is None:
            return root(
                num.nansum(exp(w*ms[num.newaxis, :]), axis=1) /
                num.nansum(exp(w*ns[num.newaxis, :]), axis=1))

        return root(num.nansum(exp(w*ms)) / num.nansum(exp(w*ns)))

    def bootstrap_misfits(self, misfits, nbootstrap, ibootstrap=None):
        exp, root = self.get_norm_functions()

        w = self.get_bootstrap_weights(
                nbootstrap, ibootstrap)[num.newaxis, :] * \
            self.get_target_weights()[num.newaxis, :] * \
            self.inter_group_weights2(misfits[:, :, 1])

        bms = root(num.nansum(exp(w*misfits[:, :, 0]), axis=1) /
                   num.nansum(exp(w*misfits[:, :, 1]), axis=1))

        # From Henriette
        # w = self.get_target_weights()[num.newaxis, :] * \
        #     self.inter_group_weights2(misfits[:, :, 1])
        # #w = self.get_bootstrap_weights(ibootstrap)[num.newaxis, :] * \
        # #    self.get_target_weights()[num.newaxis, :] * \
        # #    self.inter_group_weights2(misfits[:, :, 1])

        # bms = num.sqrt(num.nansum((w*misfits[:, :, 0])**2, axis=1))
        return bms

    def global_misfit(self, misfits):
        exp, root = self.get_norm_functions()

        ws = self.get_target_weights() * \
            self.inter_group_weights(misfits[:, 1])

        return root(num.nansum(exp(ws*misfits[:, 0])) /
                    num.nansum(exp(ws*misfits[:, 1])))

    def global_misfits(self, misfits):
        exp, root = self.get_norm_functions()
        ws = self.get_target_weights()[num.newaxis, :] * \
            self.inter_group_weights2(misfits[:, :, 1])
        return root(num.nansum(exp(ws*misfits[:, :, 0]), axis=1) /
                    num.nansum(exp(ws*misfits[:, :, 1]), axis=1))

    def global_contributions(self, misfits):
        exp, root = self.get_norm_functions()
        ws = self.get_target_weights()[num.newaxis, :] * \
            self.inter_group_weights2(misfits[:, :, 1])

        gcms = exp(ws*misfits[:, :, 0]) / \
            num.nansum(exp(ws*misfits[:, :, 1]), axis=1)[:, num.newaxis]

        return gcms

    def make_group_mask(self):
        family_names = set()
        families = num.zeros(len(self.targets), dtype=num.int)

        for itarget, target in enumerate(self.targets):
            family_names.add(target.normalisation_family)
            families[itarget] = len(family_names) - 1

        return families, len(family_names)

    def get_group_mask(self):
        if self._group_mask is None:
            self._group_mask = self.make_group_mask()

        return self._group_mask

    def evaluate(self, x, dump_data=False):
        raise NotImplementedError()

    def forward(self, x):
        source = self.get_source(x)
        engine = self.get_engine()
        plain_targets = [target.get_plain_target() for target in self.targets]

        resp = engine.process(source, plain_targets)
        results = []
        for target, result in zip(self.targets, resp.results_list[0]):
            if isinstance(result, gf.SeismosizerError):
                logger.debug(
                    '%s.%s.%s.%s: %s' % (target.codes + (str(result),)))

                results.append(None)
            else:
                results.append(result)

        return results


class ModelHistory(object):

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
            models, misfits = load_problem_data(path, problem)
            self.extend(models, misfits)

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
                (nmodels_capacity_new, self.problem.ntargets, 2),
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
            for i in xrange(n):
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

    def add_listener(self, listener):
        self.listeners.append(listener)

    def emit(self, event_name, *args, **kwargs):
        for listener in self.listeners:
            getattr(listener, event_name)(*args, **kwargs)


def load_problem_info_and_data(dirname, subset=None):
    problem = load_problem_info(dirname)
    xs, misfits = load_problem_data(xjoin(dirname, subset), problem)
    return problem, xs, misfits


def load_problem_info(dirname):
    fn = op.join(dirname, 'problem.yaml')
    return guts.load(filename=fn)


def load_problem_data(dirname, problem, skip_models=0):
    fn = op.join(dirname, 'models')
    with open(fn, 'r') as f:
        nmodels = os.fstat(f.fileno()).st_size // (problem.nparameters * 8)
        nmodels -= skip_models
        f.seek(skip_models * problem.nparameters * 8)
        data1 = num.fromfile(
            f, dtype='<f8',
            count=nmodels * problem.nparameters)\
            .astype(num.float)

    nmodels = data1.size // problem.nparameters - skip_models
    xs = data1.reshape((nmodels, problem.nparameters))

    fn = op.join(dirname, 'misfits')
    with open(fn, 'r') as f:
        f.seek(skip_models * problem.ntargets * 2 * 8)
        misfits = num.fromfile(
            f, dtype='<f8', count=nmodels*problem.ntargets*2).astype(num.float)

    misfits = misfits.reshape((nmodels, problem.ntargets, 2))

    return xs, misfits


__all__ = '''
    Problem
    ModelHistory
    ProblemConfig
    load_problem_info
    load_problem_info_and_data
'''.split()
