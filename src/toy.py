
import numpy as num

from pyrocko.guts import Object, Float, Dict, List, String, Int
from pyrocko import gf
from grond import Problem, Parameter, MisfitTarget
from grond.optimisers.highscore.plot import HighScoreOptimiserPlot

guts_prefix = 'grond.toy'


class ToyOptimiserPlot(HighScoreOptimiserPlot):

    def set_source(self, source):
        self._source = source

    def set_targets(self, targets):
        self._targets = targets

    def set_contour_data(self, contour_data):
        self._contour_data = contour_data

    def set_limits(self):
        self.axes.set_xlim(-10., 10.)
        self.axes.set_ylim(-10., 10.)

    def start(self):
        HighScoreOptimiserPlot.start(self)
        x = [getattr(t, self.xpar_name) for t in self._targets]
        y = [getattr(t, self.ypar_name) for t in self._targets]
        self.axes.plot(x, y, '^', color='black')

        for ibootstrap, (xc, yc, zc) in enumerate(
                self._contour_data['east', 'depth']):

            zmin = num.min(zc)

            if self.optimiser.nbootstrap < 5:
                alpha = 1.0
            else:
                alpha = 0.5

            self.axes.contour(
                xc, yc, zc, [zmin + 0.01],
                colors=[self.bcolors[ibootstrap]], alpha=alpha)

        self.axes.plot(
            getattr(self._source, self.xpar_name),
            getattr(self._source, self.ypar_name),
            '*', color='black')


class ToyTarget(MisfitTarget):
    north = Float.T()
    east = Float.T()
    depth = Float.T()
    obs_distance = Float.T()
    nmisfits = Int.T(default=1)


class ToySource(Object):
    north = Float.T()
    east = Float.T()
    depth = Float.T()


class ToyProblem(Problem):
    problem_parameters = [
        Parameter('north', 'm', label='North'),
        Parameter('east', 'm', label='East'),
        Parameter('depth', 'm', label='Depth')]

    ranges = Dict.T(String.T(), gf.Range.T())

    targets = List.T(ToyTarget.T())
    base_source = ToySource.T()

    def __init__(self, **kwargs):
        Problem.__init__(self, **kwargs)
        self._xtargets = None
        self._obs_distances = None

    def pack(self, source):
        return num.array(
            [source.north, source.east, source.depth], dtype=num.float)

    def _setup_modelling(self):
        if self._xtargets is None:
            self._xtargets = num.array(
                [(t.north, t.east, t.depth) for t in self.targets],
                dtype=num.float)

            self._obs_distances = num.array(
                [t.obs_distance for t in self.targets],
                dtype=num.float)

    def evaluate(self, x):
        raise NotImplementedError('Toy problem does not have evaluate()')

    def misfits(self, x, mask=None):
        self._setup_modelling()
        distances = num.sqrt(
            num.sum((x[num.newaxis, :]-self._xtargets)**2, axis=1))

        misfits = num.zeros((self.ntargets, 2))
        misfits[:, 0] = num.abs(distances - self._obs_distances)
        misfits[:, 1] = num.ones(self.ntargets) \
            * num.mean(num.abs(self._obs_distances))
        return misfits

    def misfits_many(self, xs):
        self._setup_modelling()
        distances = num.sqrt(
            num.sum(
                (xs[:, num.newaxis, :]-self._xtargets[num.newaxis, :])**2,
                axis=2))

        misfits = num.zeros((xs.shape[0], self.ntargets, 2))

        misfits[:, :, 0] = num.abs(
            distances - self._obs_distances[num.newaxis, :])

        misfits[:, :, 1] = num.mean(num.abs(self._obs_distances))

        return misfits

    def xref(self):
        base_source = self.base_source
        return num.array([
            base_source.north, base_source.east, base_source.depth])

    def extract(self, xs, i):
        if xs.ndim == 1:
            return self.extract(xs[num.newaxis, :], i)[0]

        if i < self.nparameters:
            return xs[:, i]
        else:
            return self.make_dependant(
                xs, self.dependants[i-self.nparameters].name)


def scenario(station_setup, noise_setup):

    snorth = 0.
    seast = 0.
    sdepth = 5.

    source = ToySource(
        north=snorth,
        east=seast,
        depth=sdepth)

    n = 10
    num.random.seed(10)
    norths = num.random.uniform(-10., 10., n)

    if station_setup == 'wellposed':
        easts = num.random.uniform(-10., 10., n)
    elif station_setup == 'illposed':
        easts = num.random.uniform(0, 10., n)
    else:
        assert False

    depths = num.zeros(n)

    distances = num.sqrt(
        (norths-snorth)**2 +
        (easts-seast)**2 +
        (depths-sdepth)**2)

    if noise_setup == 'noisefree':
        measured_distances = distances
    elif noise_setup == 'lownoise':
        measured_distances = distances + num.random.normal(scale=0.4, size=n)
    elif noise_setup == 'highnoise':
        measured_distances = distances + num.random.normal(scale=0.8, size=n)

    targets = [
        ToyTarget(
            path='t%03i' % i,
            north=float(norths[i]),
            east=float(easts[i]),
            depth=float(depths[i]),
            obs_distance=float(measured_distances[i]))
        for i in range(n)]

    return source, targets


__all__ = '''
    ToyTarget
    ToySource
    ToyProblem
    scenario
'''.split()
