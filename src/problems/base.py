import numpy as num
import copy
import logging
import os.path as op

from pyrocko import gf, util, guts
from pyrocko.guts import Object, String, Bool, List, Dict, Int

from ..meta import ADict, Parameter, GrondError
from ..targets import MisfitTarget, MisfitSatelliteTarget


guts_prefix = 'grond'
logger = logging.getLogger('grond')
km = 1e3
as_km = dict(scale_factor=km, scale_unit='km')


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
            self, dirname, x, ms, ns,
            accept=None, ibootstrap_choice=None, ibase=None):

        fn = op.join(dirname, 'models')
        if not isinstance(x, num.ndarray):
            x = num.array(x)
        with open(fn, 'ab') as f:
            x.astype('<f8').tofile(f)

        fn = op.join(dirname, 'misfits')
        with open(fn, 'ab') as f:
            ms.astype('<f8').tofile(f)
            ns.astype('<f8').tofile(f)

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
                if isinstance(t, MisfitSatelliteTarget)]

    @property
    def waveform_targets(self):
        return [t for t in self.targets
                if isinstance(t, MisfitTarget)]

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

    def make_bootstrap_weights(self, nbootstrap, type='classic'):
        ntargets = self.ntargets
        ws = num.zeros((nbootstrap, ntargets))
        rstate = num.random.RandomState(23)
        for ibootstrap in xrange(nbootstrap):
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

    def get_bootstrap_weights(self, ibootstrap=None):
        if self._bootstrap_weights is None:
            self._bootstrap_weights = self.make_bootstrap_weights(
                self.nbootstrap, type='bayesian')

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
        for igroup in xrange(ngroups):
            mask = group == igroup
            ws[mask] = 1.0 / root(num.nansum(exp(ns[mask])))

        return ws

    def inter_group_weights2(self, ns):
        exp, root = self.get_norm_functions()

        group, ngroups = self.get_group_mask()
        ws = num.zeros(ns.shape)
        for igroup in xrange(ngroups):
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
        return out

    def get_dependant_bounds(self):
        return []

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

    def bootstrap_misfit(self, ms, ns, ibootstrap=None):
        exp, root = self.get_norm_functions()

        w = self.get_bootstrap_weights(ibootstrap) * \
            self.get_target_weights() * self.inter_group_weights(ns)

        if ibootstrap is None:
            return root(
                num.nansum(exp(w*ms[num.newaxis, :]), axis=1) /
                num.nansum(exp(w*ns[num.newaxis, :]), axis=1))
        else:
            return root(num.nansum(exp(w*ms)) / num.nansum(exp(w*ns)))

    def bootstrap_misfits(self, misfits, ibootstrap):
        exp, root = self.get_norm_functions()

        w = self.get_bootstrap_weights(ibootstrap)[num.newaxis, :] * \
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

    def global_misfit(self, ms, ns):
        exp, root = self.get_norm_functions()

        ws = self.get_target_weights() * self.inter_group_weights(ns)
        m = root(num.nansum(exp(ws*ms)) / num.nansum(exp(ws*ns)))
        return m

    def global_misfits(self, misfits):
        exp, root = self.get_norm_functions()
        ws = self.get_target_weights()[num.newaxis, :] * \
            self.inter_group_weights2(misfits[:, :, 1])
        gms = root(num.nansum(exp(ws*misfits[:, :, 0]), axis=1) /
                   num.nansum(exp(ws*misfits[:, :, 1]), axis=1))
        return gms

    def global_contributions(self, misfits):
        exp, root = self.get_norm_functions()
        ws = self.get_target_weights()[num.newaxis, :] * \
            self.inter_group_weights2(misfits[:, :, 1])

        gcms = exp(ws*misfits[:, :, 0]) / \
            num.nansum(exp(ws*misfits[:, :, 1]), axis=1)[:, num.newaxis]

        return gcms

    def make_group_mask(self):
        super_group_names = set()
        groups = num.zeros(len(self.targets), dtype=num.int)
        ngroups = 0
        for itarget, target in enumerate(self.targets):
            if target.super_group not in super_group_names:
                super_group_names.add(target.super_group)
                ngroups += 1

            groups[itarget] = ngroups - 1

        return groups, ngroups

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
