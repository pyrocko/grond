import numpy as num
import logging

from pyrocko import gf
from pyrocko.guts import String, Bool, Float, Dict, Int

from .base import Problem, ProblemConfig
from ..meta import expand_template, Parameter

guts_prefix = 'grond'
logger = logging.getLogger('grond.problems').getChild('rectangular')
km = 1e3
as_km = dict(scale_factor=km, scale_unit='km')


class RectangularProblemConfig(ProblemConfig):

    ranges = Dict.T(String.T(), gf.Range.T())
    apply_balancing_weights = Bool.T(default=False)
    decimation_factor = Int.T(default=1)
    distance_min = Float.T(default=0.)

    def get_problem(self, event, targets):
        base_source = gf.RectangularSource(
            lat=event.lat,
            lon=event.lon,
            time=event.time,
            depth=event.depth,
            anchor='top',
            decimation_factor=self.decimation_factor,
            )

        problem = RectangularProblem(
            name=expand_template(self.name_template, event.name),
            apply_balancing_weights=self.apply_balancing_weights,
            base_source=base_source,
            distance_min=self.distance_min,
            targets=targets,
            ranges=self.ranges,
            norm_exponent=self.norm_exponent)

        return problem


class RectangularProblem(Problem):
    # nucleation_x
    # nucleation_y
    # time
    # stf

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

    problem_waveform_parameters = [
        Parameter('nucleation_x', 'offset', label='Nucleation X'),
        Parameter('nucleation_y', 'offset', label='Nucleation Y'),
        Parameter('time', 's', label='Time'),
    ]

    dependants = []

    distance_min = Float.T(default=0.0)

    def pack(self, source):
        arr = self.get_parameter_array(source)
        for ip, p in enumerate(self.parameters):
            if p.name == 'time':
                arr[ip] -= self.base_source.time
        return arr

    def get_source(self, x):
        d = self.get_parameter_dict(x)
        if 'time' in d.keys():
            d['time'] += self.base_source['time']
        source = self.base_source.clone(**d)
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

        if mask is not None:
            assert len(mask) == len(self.targets)
            targets_ok = [
                target for (target, ok) in zip(self.targets, mask) if ok]
        else:
            targets_ok = self.targets

        self.set_target_parameter_values(x)
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

    def forward(self, x, nprocs=0):
        source = self.get_source(x)
        engine = self.get_engine()
        plain_targets = [target.get_plain_target() for target in self.targets]

        resp = engine.process(source, plain_targets, nprocs=nprocs)
        results = []
        for target, result in zip(self.targets, resp.results_list[0]):
            if isinstance(result, gf.SeismosizerError):
                logger.debug(
                    '%s.%s.%s.%s: %s' % (target.codes + (str(result),)))

                results.append(None)
            else:
                results.append(result)

        return results
