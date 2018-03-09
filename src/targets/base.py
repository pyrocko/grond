import copy

import numpy as num

from pyrocko import gf
from pyrocko.guts_array import Array
from pyrocko.guts import Object, Float


guts_prefix = 'grond'


class TargetGroup(Object):
    normalisation_family = gf.StringID.T(optional=True)
    path = gf.StringID.T(optional=True)
    weight = Float.T(default=1.0)

    interpolation = gf.InterpolationMethod.T()
    store_id = gf.StringID.T(optional=True)

    def get_targets(self, ds, event, default_path):
        raise NotImplementedError()


class TargetAnalysisResult(Object):
    class NoResult(Exception):
        pass

    balancing_weight = Float.T()


class MisfitResult(Object):
    misfits = Array.T(
        shape=(None, 2),
        dtype=num.float)


class MisfitTarget(Object):

    manual_weight = Float.T(
        default=1.0,
        help='Relative weight of this target')
    analysis_result = TargetAnalysisResult.T(
        optional=True)
    normalisation_family = gf.StringID.T(
        optional=True,
        help='Normalisation family of this misfit target')
    path = gf.StringID.T(
        help='A path identifier used for plotting')

    def __init__(self, **kwargs):
        Object.__init__(self, **kwargs)
        self.parameters = []

        self._ds = None
        self._result_mode = 'sparse'

        self._target_parameters = None
        self._target_ranges = None

    def set_dataset(self, ds):
        self._ds = ds

    def get_dataset(self):
        return self._ds

    @property
    def nmisfits(self):
        return 1

    @property
    def nparameters(self):
        return len(self._target_parameters)

    @property
    def target_parameters(self):
        if self._target_parameters is None:
            self._target_parameters = copy.deepcopy(self.parameters)
            for p in self._target_parameters:
                p.set_groups([self.id])
        return self._target_parameters

    @property
    def target_ranges(self):
        return {}

    def set_parameter_values(self, model):
        for i, p in enumerate(self.parameters):
            self.parameter_values[p.name_nogroups] = model[i]

    def set_result_mode(self, result_mode):
        self._result_mode = result_mode

    def string_id(self):
        return '.'.join([self.path, self.id])

    def post_process(self, engine, source, statics):
        raise NotImplementedError()

    def get_combined_weight(self, apply_balancing_weights=False):
        return num.ones(1, dtype=num.float)

    def init_modelling(self):
        return []

    def finalize_modelling(self, modelling_results):
        raise NotImplemented('must be overloaded in subclass')


__all__ = '''
    TargetGroup
    MisfitTarget
    MisfitResult
    TargetAnalysisResult
'''.split()
