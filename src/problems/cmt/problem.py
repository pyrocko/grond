import numpy as num
import math
import logging
import copy

from pyrocko import gf, util, moment_tensor as mtm
from pyrocko.guts import String, Float, Dict, StringChoice, Int, clone

from grond.meta import GrondError, Forbidden, expand_template, Parameter, \
    has_get_plot_classes

from ..base import Problem, ProblemConfig

guts_prefix = 'grond'
logger = logging.getLogger('grond.problems.cmt.problem')
km = 1e3
as_km = dict(scale_factor=km, scale_unit='km')


class MTType(StringChoice):
    choices = ['full', 'deviatoric', 'dc']


class STFType(StringChoice):
    choices = ['HalfSinusoidSTF', 'ResonatorSTF']

    cls = {
        'HalfSinusoidSTF': gf.HalfSinusoidSTF,
        'ResonatorSTF': gf.ResonatorSTF}

    @classmethod
    def base_stf(cls, name):
        return cls.cls[name]()


class CMTProblemConfig(ProblemConfig):

    ranges = Dict.T(String.T(), gf.Range.T())
    distance_min = Float.T(default=0.0)
    mt_type = MTType.T(default='full')
    stf_type = STFType.T(default='HalfSinusoidSTF')
    nthreads = Int.T(default=1)

    def get_problem(self, event_group, target_groups, targets):
        if len(event_group.get_events()) != 1:
            raise GrondError('CMTProblem cannot handle multi-events.')

        event = copy.deepcopy(event_group.get_events()[0])

        if event.depth is None:
            event.depth = 0

        base_source = gf.MTSource.from_pyrocko_event(event)

        stf = STFType.base_stf(self.stf_type)
        stf.duration = event.duration or 0.0

        base_source.stf = stf

        subs = dict(
            event_name=event.name,
            event_time=util.time_to_str(event.time))

        problem = CMTProblem(
            name=expand_template(self.name_template, subs),
            base_source=base_source,
            target_groups=target_groups,
            targets=targets,
            ranges=self.ranges,
            distance_min=self.distance_min,
            mt_type=self.mt_type,
            stf_type=self.stf_type,
            norm_exponent=self.norm_exponent,
            nthreads=self.nthreads)

        return problem


@has_get_plot_classes
class CMTProblem(Problem):

    problem_parameters = [
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
        Parameter('rmed', label='$m_{ed} / M_0$')]

    problem_parameters_stf = {
        'HalfSinusoidSTF': [
            Parameter('duration', 's', label='Duration')],
        'ResonatorSTF': [
            Parameter('duration', 's', label='Duration'),
            Parameter('frequency', 'Hz', label='Frequency')]}

    dependants = [
        Parameter('strike1', u'\u00b0', label='Strike 1'),
        Parameter('dip1', u'\u00b0', label='Dip 1'),
        Parameter('rake1', u'\u00b0', label='Rake 1'),
        Parameter('strike2', u'\u00b0', label='Strike 2'),
        Parameter('dip2', u'\u00b0', label='Dip 2'),
        Parameter('rake2', u'\u00b0', label='Rake 2'),
        Parameter('rel_moment_iso', label='$M_{0}^{ISO}/M_{0}$'),
        Parameter('rel_moment_clvd', label='$M_{0}^{CLVD}/M_{0}$')]

    distance_min = Float.T(default=0.0)
    mt_type = MTType.T(default='full')
    stf_type = STFType.T(default='HalfSinusoidSTF')

    def __init__(self, **kwargs):
        Problem.__init__(self, **kwargs)
        self.deps_cache = {}
        self.problem_parameters = self.problem_parameters \
            + self.problem_parameters_stf[self.stf_type]
        self._base_stf = STFType.base_stf(self.stf_type)

    def get_stf(self, d):
        d_stf = {}
        for p in self.problem_parameters_stf[self.stf_type]:
            d_stf[p.name] = float(d[p.name])

        return self._base_stf.clone(**d_stf)

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

        return self.base_source.clone(m6=m6, stf=self.get_stf(d), **p)

    def make_dependant(self, xs, pname):
        cache = self.deps_cache
        if xs.ndim == 1:
            return self.make_dependant(xs[num.newaxis, :], pname)[0]

        if pname not in self.dependant_names:
            raise KeyError(pname)

        mt = self.base_source.pyrocko_moment_tensor()

        sdrs_ref = mt.both_strike_dip_rake()

        y = num.zeros(xs.shape[0])
        for i, x in enumerate(xs):
            k = tuple(x.tolist())
            if k not in cache:
                source = self.get_source(x)
                mt = source.pyrocko_moment_tensor()
                res = mt.standard_decomposition()
                sdrs = mt.both_strike_dip_rake()
                if sdrs_ref:
                    sdrs = mtm.order_like(sdrs, sdrs_ref)

                cache[k] = mt, res, sdrs

            mt, res, sdrs = cache[k]

            if pname == 'rel_moment_iso':
                ratio_iso, m_iso = res[0][1:3]
                y[i] = ratio_iso * num.sign(m_iso[0, 0])
            elif pname == 'rel_moment_clvd':
                ratio_clvd, m_clvd = res[2][1:3]
                evals, evecs = mtm.eigh_check(m_clvd)
                ii = num.argmax(num.abs(evals))
                y[i] = ratio_clvd * num.sign(evals[ii])
            else:
                isdr = {'strike': 0, 'dip': 1, 'rake': 2}[pname[:-1]]
                y[i] = sdrs[int(pname[-1])-1][isdr]

        return y

    def pack_stf(self, stf):
        return [
            stf[p.name] for p in self.problem_parameters_stf[self.stf_type]]

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
            ] + rm6.tolist() + self.pack_stf(source.stf), dtype=num.float)

        return x

    def random_uniform(self, xbounds, rstate):

        x = num.zeros(self.nparameters)
        for i in range(self.nparameters):
            x[i] = rstate.uniform(xbounds[i, 0], xbounds[i, 1])

        x[5:11] = mtm.random_m6(x=rstate.random_sample(6))

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

        elif self.mt_type == 'dc':
            mt = mtm.MomentTensor(m=m9)
            m9 = mt.standard_decomposition()[1][2]

        m0_unscaled = math.sqrt(num.sum(m9.A**2)) / math.sqrt(2.)

        m9 /= m0_unscaled
        m6 = mtm.to6(m9)
        d.rmnn, d.rmee, d.rmdd, d.rmne, d.rmnd, d.rmed = m6
        x = self.get_parameter_array(d)

        source = self.get_source(x)
        for t in self.waveform_targets:
            if t.distance_to(source) < self.distance_min:
                raise Forbidden()

        return x

    def get_dependant_bounds(self):
        out = [
            (0., 360.),
            (0., 90.),
            (-180., 180.),
            (0., 360.),
            (0., 90.),
            (-180., 180.),
            (-1., 1.),
            (-1., 1.)]

        return out

    @classmethod
    def get_plot_classes(cls):
        from .. import plot
        plots = super(CMTProblem, cls).get_plot_classes()
        plots.extend([plot.HudsonPlot, plot.MTDecompositionPlot,
                      plot.MTLocationPlot, plot.MTFuzzyPlot])
        return plots


class MultiCMTProblemConfig(ProblemConfig):

    ranges = Dict.T(String.T(), gf.Range.T())
    distance_min = Float.T(default=0.0)
    mt_type = MTType.T(default='full')
    stf_type = STFType.T(default='HalfSinusoidSTF')
    nthreads = Int.T(default=1)

    def need_event_group(self):
        return True

    def get_problem(self, event_group, target_groups, targets):
        sources = []
        for event in event_group.get_events():
            source = gf.MTSource.from_pyrocko_event(event)
            stf = STFType.base_stf(self.stf_type)
            stf.duration = event.duration or 0.0
            source.stf = stf
            sources.append(source)

        base_source = gf.CombiSource(
            name=event_group.name,
            subsources=sources)

        subs = dict(
            event_name=event_group.name)

        problem = MultiCMTProblem(
            name=expand_template(self.name_template, subs),
            base_source=base_source,
            target_groups=target_groups,
            targets=targets,
            ranges=self.ranges,
            distance_min=self.distance_min,
            mt_type=self.mt_type,
            stf_type=self.stf_type,
            norm_exponent=self.norm_exponent,
            nthreads=self.nthreads)

        return problem


@has_get_plot_classes
class MultiCMTProblem(Problem):

    problem_parameters_single = [
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
        Parameter('rmed', label='$m_{ed} / M_0$')]

    problem_parameters_stf_single = {
        'HalfSinusoidSTF': [
            Parameter('duration', 's', label='Duration')],
        'ResonatorSTF': [
            Parameter('duration', 's', label='Duration'),
            Parameter('frequency', 'Hz', label='Frequency')]}

    dependants_single = [
        Parameter('strike1', u'\u00b0', label='Strike 1'),
        Parameter('dip1', u'\u00b0', label='Dip 1'),
        Parameter('rake1', u'\u00b0', label='Rake 1'),
        Parameter('strike2', u'\u00b0', label='Strike 2'),
        Parameter('dip2', u'\u00b0', label='Dip 2'),
        Parameter('rake2', u'\u00b0', label='Rake 2'),
        Parameter('rel_moment_iso', label='$M_{0}^{ISO}/M_{0}$'),
        Parameter('rel_moment_clvd', label='$M_{0}^{CLVD}/M_{0}$')]

    distance_min = Float.T(default=0.0)
    mt_type = MTType.T(default='full')
    stf_type = STFType.T(default='HalfSinusoidSTF')

    def __init__(self, **kwargs):
        Problem.__init__(self, **kwargs)
        self.deps_cache = {}

        parameters = []
        dependants = []
        for source in self.base_source.subsources:
            for bparameter in self.problem_parameters_single \
                    + self.problem_parameters_stf_single[self.stf_type]:
                parameter = clone(bparameter)
                parameter.name = '.'.join([source.name, bparameter.name])
                parameters.append(parameter)

            for bparameter in self.dependants_single:
                parameter = clone(bparameter)
                parameter.name = '.'.join([source.name, bparameter.name])
                dependants.append(parameter)

        self.nsubsources = len(self.base_source.subsources)

        self._base_stf = STFType.base_stf(self.stf_type)

        self.problem_parameters = parameters
        self.dependants = dependants

    def get_range(self, k):
        try:
            return self.ranges[k]
        except KeyError:
            return self.ranges[k.split('.')[1]]
        except (IndexError, KeyError):
            raise GrondError('Invalid range key: %s' % k)

    def get_stf(self, subsource_name, d):
        d_stf = {}
        for p in self.problem_parameters_stf_single[self.stf_type]:
            d_stf[p.name] = float(d[subsource_name+'.'+p.name])

        return self._base_stf.clone(**d_stf)

    def get_source(self, x):
        d = self.get_parameter_dict(x)
        subsources = []
        for source in self.base_source.subsources:
            n = source.name

            rm6 = num.array([
                d[n+'.rmnn'],
                d[n+'.rmee'],
                d[n+'.rmdd'],
                d[n+'.rmne'],
                d[n+'.rmnd'],
                d[n+'.rmed']], dtype=num.float)

            m0 = mtm.magnitude_to_moment(d[n+'.magnitude'])
            m6 = rm6 * m0

            p = {}
            for k in source.keys():
                if k in d:
                    p[k] = float(
                        self.get_range(n+'.'+k).make_relative(
                            source[k], d[n+'.'+k]))

            csource = source.clone(m6=m6, stf=self.get_stf(n, d), **p)
            subsources.append(csource)

        source = self.base_source.clone(subsources=subsources)
        return source

    def make_dependant(self, xs, pname):
        cache = self.deps_cache
        if xs.ndim == 1:
            return self.make_dependant(xs[num.newaxis, :], pname)[0]

        if pname not in self.dependant_names:
            raise KeyError(pname)

        sdrs_ref_list = []
        for source in self.base_source.subsources:
            mt = source.pyrocko_moment_tensor()
            sdrs_ref_list.append(mt.both_strike_dip_rake())

        pname_ename, pname_par = pname.split('.')
        ename_to_isub = dict(
            (s.name, i) for (i, s) in enumerate(self.base_source.subsources))

        y = num.zeros(xs.shape[0])
        for i, x in enumerate(xs):
            k = tuple(x.tolist())
            if k not in cache:
                source = self.get_source(x)
                cdata = []
                for isub, subsource in enumerate(source.subsources):
                    mt = subsource.pyrocko_moment_tensor()
                    res = mt.standard_decomposition()
                    sdrs = mt.both_strike_dip_rake()
                    if sdrs_ref_list[isub]:
                        sdrs = mtm.order_like(sdrs, sdrs_ref_list[isub])
                    cdata.append((res, sdrs))

                cache[k] = cdata

            res, sdrs = cache[k][ename_to_isub[pname_ename]]

            if pname_par == 'rel_moment_iso':
                ratio_iso, m_iso = res[0][1:3]
                y[i] = ratio_iso * num.sign(m_iso[0, 0])
            elif pname_par == 'rel_moment_clvd':
                ratio_clvd, m_clvd = res[2][1:3]
                evals, evecs = mtm.eigh_check(m_clvd)
                ii = num.argmax(num.abs(evals))
                y[i] = ratio_clvd * num.sign(evals[ii])
            else:
                isdr = {'strike': 0, 'dip': 1, 'rake': 2}[pname_par[:-1]]
                y[i] = sdrs[int(pname_par[-1])-1][isdr]

        return y

    def pack_stf(self, stf):
        return [
            stf[p.name]
            for p in self.problem_parameters_stf_single[self.stf_type]]

    def pack(self, source):

        xs = []
        for isub, subsource in enumerate(source.subsources):
            m6 = subsource.m6
            mt = subsource.pyrocko_moment_tensor()
            rm6 = m6 / mt.scalar_moment()
            xs.append(num.array([
                subsource.time - self.base_source.subsources[isub].time,
                subsource.north_shift,
                subsource.east_shift,
                subsource.depth,
                mt.moment_magnitude()]
                + rm6.tolist()
                + self.pack_stf(subsource.stf), dtype=num.float))

        x = num.concatenate(xs)
        return x

    def random_uniform(self, xbounds, rstate):

        x = num.zeros(self.nparameters)
        for i in range(self.nparameters):
            x[i] = rstate.uniform(xbounds[i, 0], xbounds[i, 1])

        x[5:11] = mtm.random_m6(x=rstate.random_sample(6))

        return x.tolist()

    def preconstrain(self, x):
        d = self.get_parameter_dict(x)
        for source in self.base_source.subsources:
            n = source.name

            m6 = num.array([
                d[n+'.rmnn'],
                d[n+'.rmee'],
                d[n+'.rmdd'],
                d[n+'.rmne'],
                d[n+'.rmnd'],
                d[n+'.rmed']], dtype=num.float)

            m9 = mtm.symmat6(*m6)
            if self.mt_type == 'deviatoric':
                trace_m = num.trace(m9)
                m_iso = num.diag([trace_m / 3., trace_m / 3., trace_m / 3.])
                m9 -= m_iso

            elif self.mt_type == 'dc':
                mt = mtm.MomentTensor(m=m9)
                m9 = mt.standard_decomposition()[1][2]

            m0_unscaled = math.sqrt(num.sum(m9.A**2)) / math.sqrt(2.)

            m9 /= m0_unscaled
            m6 = mtm.to6(m9)

            (d[n+'.rmnn'], d[n+'.rmee'], d[n+'.rmdd'],
                d[n+'.rmne'], d[n+'.rmnd'], d[n+'.rmed']) = m6

        x = self.get_parameter_array(d)

        source = self.get_source(x)
        for t in self.waveform_targets:
            for subsource in source.subsources:
                if t.distance_to(subsource) < self.distance_min:
                    raise Forbidden()

        return x

    def get_dependant_bounds(self):
        out = [
            (0., 360.),
            (0., 90.),
            (-180., 180.),
            (0., 360.),
            (0., 90.),
            (-180., 180.),
            (-1., 1.),
            (-1., 1.)] * self.nsubsources

        return out

    @classmethod
    def get_plot_classes(cls):
        from .. import plot  # noqa
        plots = super(MultiCMTProblem, cls).get_plot_classes()
        return plots


__all__ = '''
    CMTProblem
    CMTProblemConfig
    MultiCMTProblem
    MultiCMTProblemConfig
'''.split()
