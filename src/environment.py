import os.path as op

from grond.config import read_config
from grond import meta
from grond.problems.base import load_optimiser_info, load_problem_info, \
    ModelHistory


class GrondEnvironmentError(meta.GrondError):
    pass


class NoEventSelected(GrondEnvironmentError):
    pass


class NoRundirAvailable(GrondEnvironmentError):
    pass


class NoPlotCollectionManagerAvailable(GrondEnvironmentError):
    pass


class Environment(object):

    def __init__(self, rundir_or_config_path=None, event_name=None):

        if op.isdir(rundir_or_config_path):
            self._rundir_path = rundir_or_config_path
            self._config_path = op.join(self._rundir_path, 'config.yaml')
        else:
            self._rundir_path = None
            self._config_path = rundir_or_config_path

        self._event_name = event_name
        self._config = None
        self._plot_collection_manager = None

        self.reset()

    def reset(self):
        self._histories = {}
        self._dataset = None
        self._optimiser = None
        self._problem = None

    def get_config(self):
        if self._config is None:
            self._config = read_config(self._config_path)

        return self._config

    def get_events_names(self):
        return self.get_config().get_event_names()

    def set_event_name(self, event_name):
        self._event_name = event_name
        self.reset()

    def get_event_name(self):
        if self._event_name is None:
            try:
                self.get_rundir_path()
                self._event_name = self.get_problem().base_source.name
            except NoRundirAvailable:
                raise NoEventSelected()

        return self._event_name

    def get_dataset(self):
        if self._dataset is None:
            event_name = self.get_event_name()
            self._dataset = self.get_config().get_dataset(event_name)

        return self._dataset

    def get_rundir_path(self):
        if self._rundir_path is None:
            raise NoRundirAvailable()

        return self._rundir_path

    def get_optimiser(self):
        if self._optimiser is None:
            try:
                self._optimiser = load_optimiser_info(self.get_rundir_path())
            except NoRundirAvailable:
                self._optimiser = \
                    self.get_config().optimiser_config.get_optimiser()

        return self._optimiser

    def get_problem(self):
        if self._problem is None:
            try:
                self._problem = load_problem_info(self.get_rundir_path())
            except NoRundirAvailable:
                self._problem = \
                    self.get_config().get_problem(self.get_event_name())

        return self._problem

    def get_history(self, subset=None):
        if subset not in self._histories:
            self._histories[subset] = \
                ModelHistory(
                    self.get_problem(),
                    path=meta.xjoin(self.get_rundir_path(), subset))

        return self._histories[subset]

    def set_plot_collection_manager(self, pcm):
        self._plot_collection_manager = pcm

    def get_plot_collection_manager(self):
        if self._plot_collection_manager is None:
            raise NoPlotCollectionManagerAvailable()

        return self._plot_collection_manager

    def setup_modelling(self):
        '''Must be called before any modelling can be done.'''

        self.get_config().setup_modelling_environment(self.get_problem())
        ds = self.get_dataset()
        for target in self.get_problem().targets:
            target.set_dataset(ds)

    def get_plot_classes(self):
        '''Discover all plot classes relevant for the setup.'''

        plots = set()
        try:
            plots.update(self.get_problem().get_plot_classes())
        except GrondEnvironmentError:
            pass

        try:
            plots.update(self.get_optimiser().get_plot_classes())
        except GrondEnvironmentError:
            pass

        try:
            for target in self.get_problem().targets:
                plots.update(target.get_plot_classes())
        except GrondEnvironmentError:
            pass

        return sorted(list(plots), key=lambda plot: plot.name)

    @property
    def path_plots(self):
        return op.join(self.get_rundir_path(), 'plots')
