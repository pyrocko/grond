import copy

import numpy as num

from pyrocko import gf
from pyrocko.guts_array import Array
from pyrocko.guts import Object, Float, String, Dict, List, Choice, load, dump

from grond.analysers.base import AnalyserResult
from grond.meta import has_get_plot_classes, GrondError


guts_prefix = 'grond'


class TargetGroup(Object):
    normalisation_family = gf.StringID.T(
        optional=True,
        help='Group with common misfit normalisation')
    path = gf.StringID.T(
        optional=True,
        help='Targets.id will be prefixed with this path')
    weight = Float.T(
        default=1.0,
        help='Additional manual weight of the target group')
    interpolation = gf.InterpolationMethod.T(
        default='nearest_neighbor',
        help='Interpolation from pre-calculated GF store.')
    store_id = gf.StringID.T(
        optional=True,
        help='ID of the Green\'s function store for this TargetGroup.')

    def get_targets(self, ds, event, default_path='none'):
        if not self._targets:
            raise NotImplementedError()


class MisfitResult(Object):
    misfits = Array.T(
        shape=(None, 2),
        dtype=num.float)


class MisfitConfig(Object):
    pass


@has_get_plot_classes
class MisfitTarget(Object):

    manual_weight = Float.T(
        default=1.0,
        help='Relative weight of this target')
    analyser_results = Dict.T(
        gf.StringID.T(),
        AnalyserResult.T(),
        help='Dictionary of analyser results')
    normalisation_family = gf.StringID.T(
        optional=True,
        help='Normalisation family of this misfit target')
    path = gf.StringID.T(
        help='A path identifier used for plotting')
    misfit_config = MisfitConfig.T(
        default=MisfitConfig.D(),
        help='Misfit configuration')
    bootstrap_weights = Array.T(
        dtype=num.float,
        serialize_as='base64',
        optional=True)
    bootstrap_residuals = Array.T(
        dtype=num.float,
        serialize_as='base64',
        optional=True)

    can_bootstrap_weights = False
    can_bootstrap_residuals = False

    plot_misfits_cumulative = True

    def __init__(self, **kwargs):
        Object.__init__(self, **kwargs)
        self.parameters = []

        self._ds = None
        self._result_mode = 'sparse'

        self._combined_weight = None
        self._target_parameters = None
        self._target_ranges = None

        self._combined_weight = None

    @classmethod
    def get_plot_classes(cls):
        return []

    def set_dataset(self, ds):
        self._ds = ds

    def get_dataset(self):
        return self._ds

    def string_id(self):
        return str(self.path)

    def misfits_string_ids(self):
        raise NotImplementedError('%s does not implement misfits_string_id'
                                  % self.__class__.__name__)

    @property
    def nmisfits(self):
        return 1

    def noise_weight_matrix(self):
        return num.array([[1]])

    @property
    def nparameters(self):
        if self._target_parameters is None:
            return 0
        return len(self._target_parameters)

    @property
    def target_parameters(self):
        if self._target_parameters is None:
            self._target_parameters = copy.deepcopy(self.parameters)
            for p in self._target_parameters:
                p.set_groups([self.string_id()])
        return self._target_parameters

    @property
    def target_ranges(self):
        return {}

    def set_parameter_values(self, model):
        for i, p in enumerate(self.parameters):
            self.parameter_values[p.name_nogroups] = model[i]

    def set_result_mode(self, result_mode):
        self._result_mode = result_mode

    def post_process(self, engine, source, statics):
        raise NotImplementedError()

    def get_combined_weight(self):
        if self._combined_weight is None:
            self._combined_weight = num.ones(1, dtype=num.float)
        return self._combined_weight

    def get_correlated_weights(self, nthreads=0):
        pass

    def set_bootstrap_weights(self, weights):
        self.bootstrap_weights = weights

    def get_bootstrap_weights(self):
        if self.bootstrap_weights is None:
            raise Exception('Bootstrap weights have not been set!')
        nbootstraps = self.bootstrap_weights.size // self.nmisfits
        return self.bootstrap_weights.reshape(nbootstraps, self.nmisfits)

    def init_bootstrap_residuals(self, nbootstrap, rstate=None, nthreads=0):
        raise NotImplementedError()

    def set_bootstrap_residuals(self, residuals):
        self.bootstrap_residuals = residuals

    def get_bootstrap_residuals(self):
        if self.bootstrap_residuals is None:
            raise Exception('Bootstrap residuals have not been set!')
        nbootstraps = self.bootstrap_residuals.size // self.nmisfits
        return self.bootstrap_residuals.reshape(nbootstraps, self.nmisfits)

    def prepare_modelling(self, engine, source, targets):
        ''' Prepare modelling target

        This function shall return a list of :class:`pyrocko.gf.Target`
        for forward modelling in the :class:`pyrocko.gf.LocalEngine`.
        '''
        return [self]

    def finalize_modelling(
            self, engine, source, modelling_targets, modelling_results):
        ''' Manipulate modelling before misfit calculation

        This function can be overloaded interact with the modelling results.
        '''
        return modelling_results[0]


class MisfitResultError(Object):
    message = String.T()


class MisfitResultCollection(Object):
    results = List.T(List.T(
        Choice.T([MisfitResult.T(), MisfitResultError.T()])))


def dump_misfit_result_collection(misfit_result_collection, path):
    dump(misfit_result_collection, filename=path)


def load_misfit_result_collection(path):
    try:
        obj = load(filename=path)

    except OSError as e:
        raise GrondError(
            'Failed to read ensemble misfit results from file "%s" (%s)' % (
                path, e))

    if not isinstance(obj, MisfitResultCollection):
        raise GrondError(
            'File "%s" does not contain any misfit result collection.' % path)

    return obj


__all__ = '''
    TargetGroup
    MisfitTarget
    MisfitResult
    MisfitResultError
    dump_misfit_result_collection
    load_misfit_result_collection
    MisfitResultCollection
'''.split()
