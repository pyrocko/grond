import numpy as num
import logging

from pyrocko import trace
from pyrocko.guts import (Object, Dict, String, Float, Bool, Int)

logger = logging.getLogger('grond.synthetic_tests')

guts_prefix = 'grond'


class SyntheticWaveformNotAvailable(Exception):
    pass


class SyntheticTest(Object):
    inject_solution = Bool.T(default=False)
    respect_data_availability = Bool.T(default=False)
    real_noise_scale = Float.T(default=0.0)
    white_noise_scale = Float.T(default=0.0)
    relative_white_noise_scale = Float.T(default=0.0)
    random_response_scale = Float.T(default=0.0)
    real_noise_toffset = Float.T(default=-3600.)
    random_seed = Int.T(optional=True)
    x = Dict.T(String.T(), Float.T())

    def __init__(self, **kwargs):
        Object.__init__(self, **kwargs)
        self._problem = None
        self._synthetics = None

    def set_problem(self, problem):
        self._problem = problem
        self._synthetics = None

    def get_problem(self):
        if self._problem is None:
            raise SyntheticWaveformNotAvailable(
                'SyntheticTest.set_problem() has not been called yet')

        return self._problem

    def get_x(self):
        problem = self.get_problem()
        if self.x:
            x = problem.preconstrain(
                problem.get_parameter_array(self.x))

        else:
            x = problem.preconstrain(
                problem.pack(
                    problem.base_source))

        return x

    def get_synthetics(self):
        problem = self.get_problem()
        if self._synthetics is None:
            x = self.get_x()
            results = problem.forward(x)
            synthetics = {}
            for iresult, result in enumerate(results):
                tr = result.trace.pyrocko_trace()
                tfade = tr.tmax - tr.tmin
                tr_orig = tr.copy()
                tr.extend(tr.tmin - tfade, tr.tmax + tfade)
                rstate = num.random.RandomState(
                    (self.random_seed or 0) + iresult)

                if self.random_response_scale != 0:
                    tf = RandomResponse(scale=self.random_response_scale)
                    tf.set_random_state(rstate)
                    tr = tr.transfer(
                        tfade=tfade,
                        transfer_function=tf)

                if self.white_noise_scale != 0.0:
                    u = rstate.normal(
                        scale=self.white_noise_scale,
                        size=tr.data_len())

                    tr.ydata += u

                if self.relative_white_noise_scale != 0.0:
                    u = rstate.normal(
                        scale=self.relative_white_noise_scale * num.std(
                            tr_orig.ydata),
                        size=tr.data_len())

                    tr.ydata += u

                synthetics[result.trace.codes] = tr

            self._synthetics = synthetics

        return self._synthetics

    def get_waveform(self, nslc, tmin, tmax, tfade=0., freqlimits=None):
        synthetics = self.get_synthetics()
        if nslc not in synthetics:
            s = 'no synthetics for %s available' % '.'.join(nslc)
            logger.warn(s)
            from grond.dataset import NotFound
            raise NotFound(s)

        tr = synthetics[nslc]
        tr.extend(tmin - tfade * 2.0, tmax + tfade * 2.0)

        tr = tr.transfer(
            tfade=tfade,
            freqlimits=freqlimits)

        tr.chop(tmin, tmax)
        return tr


class RandomResponse(trace.FrequencyResponse):

    scale = Float.T(default=0.0)

    def set_random_state(self, rstate):
        self._rstate = rstate

    def evaluate(self, freqs):
        n = freqs.size
        return 1.0 + (
            self._rstate.normal(scale=self.scale, size=n) +
            0.0J * self._rstate.normal(scale=self.scale, size=n))


__all__ = '''
SyntheticTest
SyntheticWaveformNotAvailable
'''.split()
