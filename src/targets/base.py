import copy

from pyrocko import gf
from pyrocko.guts import Object, Float


class MisfitConfig(Object):
    pass


class TargetGroup(Object):
    normalisation_family = gf.StringID.T(default='', optional=True)
    path = gf.StringID.T(optional=True)
    weight = Float.T(default=1.0)

    misfit_config = MisfitConfig.T(optional=True)

    interpolation = gf.InterpolationMethod.T()
    store_id = gf.StringID.T(optional=True)

    def get_targets(self, ds, event, default_path):
        raise NotImplementedError()


class MisfitTarget(object):
    def set_dataset(self, ds):
        self._ds = ds

    def get_dataset(self):
        return self._ds

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
        if self._target_ranges is None:
            self._target_ranges = self.misfit_config.ranges.copy()
            for k in self._target_ranges.keys():
                self._target_ranges['%s:%s' % (self.id, k)] =\
                    self._target_ranges.pop(k)
        return self._target_ranges

    def set_parameter_values(self, model):
        for i, p in enumerate(self.parameters):
            self.parameter_values[p.name_nogroups] = model[i]

    def set_result_mode(self, result_mode):
        self._result_mode = result_mode

    def string_id(self):
        return '.'.join([self.normalisation_family, self.path, self.id])

    def post_process(self, engine, source, statics):
        raise NotImplementedError()

    def get_combined_weight(self, apply_balancing_weights=False):
        raise NotImplementedError()


class MisfitResult(gf.Result):
    misfit_value = Float.T()
    misfit_norm = Float.T()


class TargetAnalysisResult(Object):
    class NoResult(Exception):
        pass

    balancing_weight = Float.T()
