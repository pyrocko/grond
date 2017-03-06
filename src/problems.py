import numpy as num
import copy
import logging
import os.path as op
import math

from meta import Forbidden, expand_template
from pyrocko import gf, util, guts, moment_tensor as mtm
from pyrocko.guts import (Object, String, Bool, List, Float, Dict, Int,
                          StringChoice)

from .targets import MisfitTarget, MisfitSatelliteTarget


guts_prefix = 'grond'
logger = logging.getLogger('grond')
km = 1e3
as_km = dict(scale_factor=km, scale_unit='km')


class Parameter(Object):
    name = String.T()
    unit = String.T(optional=True)
    scale_factor = Float.T(default=1., optional=True)
    scale_unit = String.T(optional=True)
    label = String.T(optional=True)

    def __init__(self, *args, **kwargs):
        if len(args) >= 1:
            kwargs['name'] = args[0]
        if len(args) >= 2:
            kwargs['unit'] = args[1]

        self.target = kwargs.pop('target', None)

        Object.__init__(self, **kwargs)

    def get_label(self, with_unit=True):
        l = [self.label or self.name]
        if with_unit:
            unit = self.get_unit_label()
            if unit:
                l.append('[%s]' % unit)

        return ' '.join(l)

    def get_value_label(self, value, format='%(value)g%(unit)s'):
        value = self.scaled(value)
        unit = self.get_unit_suffix()
        return format % dict(value=value, unit=unit)

    def get_unit_label(self):
        if self.scale_unit is not None:
            return self.scale_unit
        elif self.unit:
            return self.unit
        else:
            return None

    def get_unit_suffix(self):
        unit = self.get_unit_label()
        if not unit:
            return ''
        else:
            return ' %s' % unit

    def scaled(self, x):
        if isinstance(x, tuple):
            return tuple(v/self.scale_factor for v in x)
        if isinstance(x, list):
            return list(v/self.scale_factor for v in x)
        else:
            return x/self.scale_factor


class ADict(dict):
    def __getattr__(self, k):
        return self[k]

    def __setattr__(self, k, v):
        self[k] = v


class ProblemConfig(Object):
    name_template = String.T()
    apply_balancing_weights = Bool.T(default=True)



class Problem(Object):
    name = String.T()
    ranges = Dict.T(String.T(), gf.Range.T())
    dependants = List.T(Parameter.T())
    apply_balancing_weights = Bool.T(default=True)
    base_source = gf.Source.T()

    targets = List.T()

    def __init__(self, **kwargs):

        Object.__init__(self, **kwargs)
        self.init_satellite_target_leveling()

        self._bootstrap_weights = None
        self._target_weights = None
        self._engine = None
        self._group_mask = None
        logger.name = self.__class__.__name__

    def get_engine(self):
        return self._engine

    def copy(self):
        o = copy.copy(self)
        o._bootstrap_weights = None
        o._target_weights = None
        return o

    @property
    def parameters(self):
        ps = self.problem_parameters
        for target in self.targets:
            ps.extend(target.target_parameters)

        return ps

    def set_target_parameter_values(self, x):
        i = len(self.problem_parameters)
        for target in self.targets:
            n = len(target.target_parameters)
            target.set_parameter_values(x[i:i+n])
            i += n

    def get_parameter_dict(self, x):
        return ADict((p.name, v) for for p, v in zip(self.parameters, x))

    def get_parameter_array(self, d):
        return num.array([d[p.name] for p in self.parameters], dtype=num.float)

    @property
    def parameter_names(self):
        return [p.name for p in self.combined]

    def dump_problem_info(self, dirname):
        fn = op.join(dirname, 'problem.yaml')
        util.ensuredirs(fn)
        guts.dump(self, filename=fn)

    def dump_problem_data(self, dirname, x, ms, ns):
        fn = op.join(dirname, 'models')
        if not isinstance(x, num.ndarray):
            x = num.array(x)
        with open(fn, 'ab') as f:
            x.astype('<f8').tofile(f)

        fn = op.join(dirname, 'misfits')
        with open(fn, 'ab') as f:
            ms.astype('<f8').tofile(f)
            ns.astype('<f8').tofile(f)

    def name_to_index(self, name):
        pnames = [p.name for p in self.combined]
        return pnames.index(name)

    @property
    def nparameters(self):
        return len(self.parameters)

    @property
    def ntargets(self):
        return len(self.targets)

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

    def make_bootstrap_weights(self, nbootstrap):
        ntargets = len(self.targets)
        ws = num.zeros((nbootstrap, ntargets))
        rstate = num.random.RandomState(23)
        for ibootstrap in xrange(nbootstrap):
            ii = rstate.randint(0, ntargets, size=self.ntargets)
            ws[ibootstrap, :] = num.histogram(
                ii, ntargets, (-0.5, ntargets - 0.5))[0]

        return ws

    def get_bootstrap_weights(self, ibootstrap=None):
        if self._bootstrap_weights is None:
            self._bootstrap_weights = self.make_bootstrap_weights(
                self.nbootstrap)

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
        group, ngroups = self.get_group_mask()

        ws = num.zeros(self.ntargets)
        for igroup in xrange(ngroups):
            mask = group == igroup
            ws[mask] = 1. / num.sqrt(num.nansum(ns[mask]**2))

        return ws

    def inter_group_weights2(self, ns):
        group, ngroups = self.get_group_mask()
        ws = num.zeros(ns.shape)
        for igroup in xrange(ngroups):
            mask = group == igroup
            ws[:, mask] = (1.0 / num.sqrt(
                num.nansum(ns[:, mask]**2, axis=1)))[:, num.newaxis]

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
                r = target.ranges[p.name]
                out.append((r.start, r.stop))

        return out

    def get_dependant_bounds(self):
        return []

    def bootstrap_misfit(self, ms, ns, ibootstrap=None):
        w = self.get_bootstrap_weights(ibootstrap) * \
            self.get_target_weights() * self.inter_group_weights(ns)

        if ibootstrap is None:
            return num.sqrt(
                num.nansum((w*ms[num.newaxis, :])**2, axis=1) /
                num.nansum((w*ns[num.newaxis, :])**2, axis=1))
        else:
            return num.sqrt(num.nansum((w*ms)**2) / num.nansum((w*ns)**2))

    def bootstrap_misfits(self, misfits, ibootstrap):
        w = self.get_bootstrap_weights(ibootstrap)[num.newaxis, :] * \
            self.get_target_weights()[num.newaxis, :] * \
            self.inter_group_weights2(misfits[:, :, 1])

        bms = num.sqrt(num.nansum((w*misfits[:, :, 0])**2, axis=1) /
                       num.nansum((w*misfits[:, :, 1])**2, axis=1))
        return bms

    def global_misfit(self, ms, ns):
        ws = self.get_target_weights() * self.inter_group_weights(ns)
        m = num.sqrt(num.nansum((ws*ms)**2) / num.nansum((ws*ns)**2))
        return m

    def global_misfits(self, misfits):
        ws = self.get_target_weights()[num.newaxis, :] * \
            self.inter_group_weights2(misfits[:, :, 1])
        gms = num.sqrt(num.nansum((ws*misfits[:, :, 0])**2, axis=1) /
                       num.nansum((ws*misfits[:, :, 1])**2, axis=1))
        return gms

    def global_contributions(self, misfits):
        ws = self.get_target_weights()[num.newaxis, :] * \
            self.inter_group_weights2(misfits[:, :, 1])

        gcms = (ws*misfits[:, :, 0])**2 / \
            num.nansum((ws*misfits[:, :, 1])**2, axis=1)[:, num.newaxis]

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


class CMTProblem(Problem):

    parameters = [
        Parameter('time', 's', label='Time'),
        Parameter('north_shift', 'm', label='Northing', **as_km),
        Parameter('east_shift', 'm', label='Easting', **as_km),
        Parameter('depth', 'm', label='Depth', **as_km),
        Parameter('magnitude', label='Magnitude'),
        Parameter('rmnn', label='$m_{nn} / M_0$'),
        Parameter('rmee', label='$m_{ee} / M_0$'),
        Parameter('rmdd', label='$m_{dd} / M_0$'),
        Parameter('rmne', label='$m_{ne} / M_0$'),
        Parameter('rmnd', label='$m_{nd} / M_0$'),
        Parameter('rmed', label='$m_{ed} / M_0$'),
        Parameter('duration', 's', label='Duration')]

    dependants = [
        Parameter('rel_moment_iso', label='$M_{0}^{ISO}/M_{0}$'),
        Parameter('rel_moment_clvd', label='$M_{0}^{CLVD}/M_{0}$')]

    distance_min = Float.T(default=0.0)
    nbootstrap = Int.T(default=10)
    mt_type = StringChoice.T(default='full', choices=['full', 'deviatoric'])

    def get_source(self, x):
        d = self.get_parameter_dict(x)
        rm6 = num.array([d.rmnn, d.rmee, d.rmdd, d.rmne, d.rmnd, d.rmed],
                        dtype=num.float)

        m0 = mtm.magnitude_to_moment(d.magnitude)
        m6 = rm6 * m0

        p = {}
        for k in self.base_source.keys():
            if k in d:
                p[k] = float(
                    self.ranges[k].make_relative(self.base_source[k], d[k]))

        stf = gf.HalfSinusoidSTF(duration=float(d.duration))

        source = self.base_source.clone(m6=m6, stf=stf, **p)
        return source

    def make_dependant(self, xs, pname):
        if xs.ndim == 1:
            return self.make_dependant(xs[num.newaxis, :], pname)[0]

        if pname in ('rel_moment_iso', 'rel_moment_clvd'):
            y = num.zeros(xs.shape[0])
            for i, x in enumerate(xs):
                source = self.get_source(x)
                mt = source.pyrocko_moment_tensor()
                res = mt.standard_decomposition()

                if pname == 'rel_moment_iso':
                    ratio_iso, m_iso = res[0][1:3]
                    y[i] = ratio_iso * num.sign(m_iso[0, 0])
                else:
                    ratio_clvd, m_clvd = res[2][1:3]
                    evals, evecs = mtm.eigh_check(m_clvd)
                    ii = num.argmax(num.abs(evals))
                    y[i] = ratio_clvd * num.sign(evals[ii])

            return y

        else:
            raise KeyError(pname)

    def extract(self, xs, i):
        if xs.ndim == 1:
            return self.extract(xs[num.newaxis, :], i)[0]

        if i < self.nparameters:
            return xs[:, i]
        else:
            return self.make_dependant(
                xs, self.dependants[i-self.nparameters].name)

    def pack(self, source):
        m6 = source.m6
        mt = source.pyrocko_moment_tensor()
        rm6 = m6 / mt.scalar_moment()

        x = num.array([
            source.time - self.base_source.time,
            source.north_shift,
            source.east_shift,
            source.depth,
            mt.moment_magnitude(),
            ] + rm6.tolist() + [source.stf.duration], dtype=num.float)

        return x

    def random_uniform(self, xbounds):
        x = num.zeros(self.nparameters)
        for i in [0, 1, 2, 3, 4, 11]:
            x[i] = num.random.uniform(xbounds[i, 0], xbounds[i, 1])

        x[5:11] = mtm.random_m6()

        return x.tolist()

    def preconstrain(self, x):

        d = self.get_parameter_dict(x)
        m6 = num.array([d.rmnn, d.rmee, d.rmdd, d.rmne, d.rmnd, d.rmed],
                       dtype=num.float)

        m9 = mtm.symmat6(*m6)
        if self.mt_type == 'deviatoric':
            trace_m = num.trace(m9)
            m_iso = num.diag([trace_m / 3., trace_m / 3., trace_m / 3.])
            m9 -= m_iso

        m0_unscaled = math.sqrt(num.sum(m9.A**2)) / math.sqrt(2.)

        m9 /= m0_unscaled
        m6 = mtm.to6(m9)
        d.rmnn, d.rmee, d.rmdd, d.rmne, d.rmnd, d.rmed = m6
        x = self.get_parameter_array(d)

        source = self.get_source(x)
        if any(self.distance_min > source.distance_to(t)
               for t in self.targets):
            raise Forbidden()

        return x

    def get_dependant_bounds(self):
        out = [
            (-1., 1.),
            (-1., 1.)]

        return out

    def evaluate(self, x, result_mode='sparse', mask=None):
        source = self.get_source(x)
        engine = self.get_engine()

        for target in self.targets:
            target.set_result_mode(result_mode)

        if mask is not None:
            assert len(mask) == len(self.targets)
            targets_ok = [
                target for (target, ok) in zip(self.targets, mask) if ok]
        else:
            targets_ok = self.targets

        resp = engine.process(source, targets_ok)

        if mask is not None:
            ires_ok = 0
            results = []
            for target, ok in zip(self.targets, mask):
                if ok:
                    results.append(resp.results_list[0][ires_ok])
                    ires_ok += 1
                else:
                    results.append(
                        gf.SeismosizerError(
                            'skipped because of previous failure'))
        else:
            results = list(resp.results_list[0])

        data = []
        for target, result in zip(self.targets, results):
            if isinstance(result, gf.SeismosizerError):
                logger.debug(
                    '%s.%s.%s.%s: %s' % (target.codes + (str(result),)))

                data.append((None, None))
            else:
                data.append((result.misfit_value, result.misfit_norm))

        ms, ns = num.array(data, dtype=num.float).T
        if result_mode == 'full':
            return ms, ns, results
        else:
            return ms, ns

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


class CMTProblemConfig(ProblemConfig):

    ranges = Dict.T(String.T(), gf.Range.T())
    distance_min = Float.T(default=0.0)
    nbootstrap = Int.T(default=10)
    mt_type = StringChoice.T(choices=['full', 'deviatoric'])

    def get_problem(self, event, targets):
        if event.depth is None:
            event.depth = 0.

        base_source = gf.MTSource.from_pyrocko_event(event)
        base_source.stf = gf.HalfSinusoidSTF(duration=event.duration or 0.0)

        subs = dict(
            event_name=event.name,
            event_time=util.time_to_str(event.time))

        problem = CMTProblem(
            name=expand_template(self.name_template, subs),
            apply_balancing_weights=self.apply_balancing_weights,
            base_source=base_source,
            targets=targets,
            ranges=self.ranges,
            distance_min=self.distance_min,
            nbootstrap=self.nbootstrap,
            mt_type=self.mt_type)

        return problem


class DoubleDCProblem(Problem):

    problem_parameters = [
        Parameter('time', 's', label='Time'),
        Parameter('north_shift', 'm', label='Northing', **as_km),
        Parameter('east_shift', 'm', label='Easting', **as_km),
        Parameter('depth', 'm', label='Depth', **as_km),
        Parameter('magnitude', label='Magnitude'),
        Parameter('strike1', 'deg', label='Strike 1'),
        Parameter('dip1', 'deg', label='Dip 1'),
        Parameter('rake1', 'deg', label='Rake 1'),
        Parameter('strike2', 'deg', label='Strike 2'),
        Parameter('dip2', 'deg', label='Dip 2'),
        Parameter('rake2', 'deg', label='Rake 2'),
        Parameter('delta_time', 's', label='$\Delta$ Time'),
        Parameter('delta_depth', 'm', label='$\Delta$ Depth'),
        Parameter('azimuth', 'deg', label='Azimuth'),
        Parameter('distance', 'm', label='Distance'),
        Parameter('mix', label='Mix'),
        Parameter('duration1', 's', label='Duration 1'),
        Parameter('duration2', 's', label='Duration 2')]

    dependants = []
    distance_min = Float.T(default=0.0)
    nbootstrap = Int.T(default=100)

    def get_source(self, x):
        d = self.get_parameter_dict(x)
        p = {}
        for k in self.base_source.keys():
            if k in d:
                p[k] = float(
                    self.ranges[k].make_relative(self.base_source[k], d[k]))

        stf1 = gf.HalfSinusoidSTF(duration=float(d.duration1))
        stf2 = gf.HalfSinusoidSTF(duration=float(d.duration2))

        source = self.base_source.clone(stf1=stf1, stf2=stf2, **p)
        return source

    def make_dependant(self, xs, pname):
        if xs.ndim == 1:
            return self.make_dependant(xs[num.newaxis, :], pname)[0]

        raise KeyError(pname)

    def extract(self, xs, i):
        if xs.ndim == 1:
            return self.extract(xs[num.newaxis, :], i)[0]

        if i < self.nparameters:
            return xs[:, i]
        else:
            return self.make_dependant(
                xs, self.dependants[i-self.nparameters].name)

    def pack(self, source):
        x = num.array([
            source.time - self.base_source.time,
            source.north_shift,
            source.east_shift,
            source.depth,
            source.magnitude,
            source.strike1,
            source.dip1,
            source.rake1,
            source.strike2,
            source.dip2,
            source.rake2,
            source.delta_time,
            source.delta_depth,
            source.azimuth,
            source.distance,
            source.mix,
            source.stf1.duration,
            source.stf2.duration], dtype=num.float)

        return x

    def random_uniform(self, xbounds):
        x = num.zeros(self.nparameters)
        for i in xrange(self.nparameters):
            x[i] = num.random.uniform(xbounds[i, 0], xbounds[i, 1])

        return x.tolist()

    def preconstrain(self, x):
        source = self.get_source(x)
        if any(self.distance_min > source.distance_to(t)
               for t in self.targets):
            raise Forbidden()

        return num.array(x, dtype=num.float)

    def evaluate(self, x, result_mode='sparse'):
        source = self.get_source(x)
        engine = self.get_engine()

        for target in self.targets:
            target.set_result_mode(result_mode)

        resp = engine.process(source, self.targets)

        data = []
        results = []
        for target, result in zip(self.targets, resp.results_list[0]):
            if isinstance(result, gf.SeismosizerError):
                logger.debug(
                    '%s.%s.%s.%s: %s' % (target.codes + (str(result),)))

                data.append((None, None))
                results.append(result)
            else:
                data.append((result.misfit_value, result.misfit_norm))
                results.append(result)

        ms, ns = num.array(data, dtype=num.float).T
        if result_mode == 'full':
            return ms, ns, results
        else:
            return ms, ns

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


class DoubleDCProblemConfig(ProblemConfig):

    ranges = Dict.T(String.T(), gf.Range.T())
    distance_min = Float.T(default=0.0)
    nbootstrap = Int.T(default=100)

    def get_problem(self, event, targets):
        if event.depth is None:
            event.depth = 0.

        base_source = gf.DoubleDCSource.from_pyrocko_event(event)
        base_source.stf = gf.HalfSinusoidSTF(duration=event.duration or 0.0)

        subs = dict(
            event_name=event.name,
            event_time=util.time_to_str(event.time))

        problem = DoubleDCProblem(
            name=expand_template(self.name_template, subs),
            apply_balancing_weights=self.apply_balancing_weights,
            base_source=base_source,
            targets=targets,
            ranges=self.ranges,
            distance_min=self.distance_min,
            nbootstrap=self.nbootstrap)

        return problem


class RectangularProblem(Problem):

    problem_parameters = [
        Parameter('north_shift', 'm', label='Northing', **as_km),
        Parameter('east_shift', 'm', label='Easting', **as_km),
        Parameter('depth', 'm', label='Depth', **as_km),
        Parameter('length', 'm', label='Length', **as_km),
        Parameter('width', 'm', label='Width', **as_km),
        Parameter('dip', 'deg', label='Dip'),
        Parameter('strike', 'deg', label='Strike'),
        Parameter('rake', 'deg', label='Rake'),
        Parameter('slip', 'm', label='Slip'),
        ]

    dependants = []

    nbootstrap = 0

    def pack(self, source):
        return self.get_parameter_array(source)

    def get_source(self, x):
        source = self.base_source.clone(**self.get_parameter_dict(x))
        return source

    def extract(self, xs, i):
        if xs.ndim == 1:
            return self.extract(xs[num.newaxis, :], i)[0]

        return xs[:, i]

    def random_uniform(self, xbounds):
        x = num.zeros(self.nparameters)
        for i in xrange(self.nparameters):
            x[i] = num.random.uniform(xbounds[i, 0], xbounds[i, 1])

        return x.tolist()

    def preconstrain(self, x):
        # source = self.get_source(x)
        # if any(self.distance_min > source.distance_to(t)
        #        for t in self.targets):
            # raise Forbidden()
        return x

    def evaluate(self, x, result_mode='sparse', mask=None, nprocs=0):
        source = self.get_source(x)
        engine = self.get_engine()

        for target in self.targets:
            target.set_result_mode(result_mode)

        self.set_satellite_scene_levels(x)

        if mask is not None:
            assert len(mask) == len(self.targets)
            targets_ok = [
                target for (target, ok) in zip(self.targets, mask) if ok]
        else:
            targets_ok = self.targets

        resp = engine.process(source, targets_ok, nprocs=nprocs)

        if mask is not None:
            ires_ok = 0
            results = []
            for target, ok in zip(self.targets, mask):
                if ok:
                    results.append(resp.results_list[0][ires_ok])
                    ires_ok += 1
                else:
                    results.append(
                        gf.SeismosizerError(
                            'skipped because of previous failure'))
        else:
            results = list(resp.results_list[0])

        data = []
        for target, result in zip(self.targets, results):
            if isinstance(result, gf.SeismosizerError):
                logger.debug(
                    '%s.%s.%s.%s: %s' % (target.codes + (str(result),)))

                data.append((None, None))
            else:
                data.append((result.misfit_value, result.misfit_norm))

        ms, ns = num.array(data, dtype=num.float).T
        if result_mode == 'full':
            return ms, ns, results
        else:
            return ms, ns

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


class RectangularProblemConfig(ProblemConfig):

    ranges = Dict.T(String.T(), gf.Range.T())
    apply_balancing_weights = Bool.T(default=False)

    def get_problem(self, event, targets):
        base_source = gf.RectangularSource(
            lat=event.lat,
            lon=event.lon,
            time=event.time,
            depth=event.depth,
            anchor='top',
            )

        problem = RectangularProblem(
            name=expand_template(self.name_template, event.name),
            apply_balancing_weights=self.apply_balancing_weights,
            base_source=base_source,
            targets=targets,
            ranges=self.ranges,)

        return problem


__all__ = [
    'CMTProblem',
    'CMTProblemConfig',
    'DoubleDCProblem',
    'DoubleDCProblemConfig',
    'RectangularProblem',
    'RectangularProblemConfig',
]
