import time
import logging
import os.path as op

from grond.config import read_config
from grond import meta, run_info
from grond.problems.base import load_optimiser_info, load_problem_info, \
    ModelHistory

logger = logging.getLogger('grond.environment')


class GrondEnvironmentError(meta.GrondError):
    pass


class EventSelectionFailed(GrondEnvironmentError):
    pass


class NoCurrentEventAvailable(GrondEnvironmentError):
    def __init__(self, message='no current event available'):
        GrondEnvironmentError.__init__(self, message)


class NoEventSelectionAvailable(GrondEnvironmentError):
    def __init__(self, message='no event selection available'):
        GrondEnvironmentError.__init__(self, message)


class NoRundirAvailable(GrondEnvironmentError):
    def __init__(self, message='no rundir available'):
        GrondEnvironmentError.__init__(self, message)


class NoPlotCollectionManagerAvailable(GrondEnvironmentError):
    def __init__(self, message='no plot collection manager available'):
        GrondEnvironmentError.__init__(self, message)


class Environment(object):

    def __init__(self, args):

        self._current_name = None
        self._selected_names = None
        self._config = None
        self._plot_collection_manager = None
        if isinstance(args, str):
            args = [args]

        self.reset()

        if not args:
            raise GrondEnvironmentError('missing arguments')

        if op.isdir(args[0]):
            self._rundir_path = args[0]
            self._config_path = op.join(self._rundir_path, 'config.yaml')

        else:
            self._rundir_path = None
            self._config_path = args[0]
            self.set_selected_names(args[1:])

    def __str__(self):
        try:
            rundir_path = self.get_rundir_path()
        except NoRundirAvailable:
            rundir_path = 'n/a'

        try:
            current_name = self.get_current_name()
        except NoCurrentEventAvailable:
            current_name = None

        return '''Grond Runtime Environment
current name:   %s
selected names: %s
rundir:         %s
config path:    %s
''' % (
            current_name or 'n/a',
            ', '.join(self._selected_names) if self._selected_names else 'n/a',
            rundir_path,
            self.get_config_path())

    @classmethod
    def discover(cls, rundir, wait=20.):
        start_watch = time.time()
        while (time.time() - start_watch) < wait:
            try:
                cls.verify_rundir(rundir)
                return cls([rundir])
            except GrondEnvironmentError:
                time.sleep(.25)

    @staticmethod
    def verify_rundir(rundir_path):
        files = [
            'config.yaml',
            'problem.yaml',
            'optimiser.yaml',
            'misfits'
            ]
        for fn in files:
            if not op.exists(op.join(rundir_path, fn)):
                raise GrondEnvironmentError('inconsistent rundir')

    def reset(self):
        self._histories = {}
        self._dataset = None
        self._optimiser = None
        self._problem = None
        self._rundir_path = None

    def get_config(self):
        if self._config is None:
            self._config = read_config(self._config_path)

        return self._config

    def get_available_names(self):
        return self.get_config().get_names()

    def names_are_groups(self):
        return self.get_config().need_event_group()

    def set_current_name(self, name):
        self._current_name = name
        self.reset()

    def get_current_name(self):
        if self._current_name is None:
            try:
                self.get_rundir_path()
                self._current_name = self.get_problem().base_source.name
            except NoRundirAvailable:
                try:
                    names = self.get_selected_names()
                    if len(names) == 1:
                        self._current_name = names[0]
                    else:
                        raise NoCurrentEventAvailable()

                except NoEventSelectionAvailable:
                    raise NoCurrentEventAvailable()

        return self._current_name

    def set_selected_names(self, args):
        names = self.get_available_names()
        sevent = 'event group' if self.names_are_groups() else 'event'

        def savail(sevent, names):
            return 'Select from available %ss:' \
                '\n    %s\n  or \'all\' to use all available %ss.' % (
                    sevent, '\n    '.join(names), sevent)

        if len(args) == 0:
            if len(names) == 1:
                self._selected_names = names
            else:
                if not names:
                    raise EventSelectionFailed(
                        'No %s, check your config!' % sevent)

                raise EventSelectionFailed(
                    'Ambiguous %s selection. %s' % (
                        (sevent, savail(sevent, names))))

        elif len(args) == 1 and args[0] == 'all':
            self._selected_names = names

        else:
            self._selected_names = []
            for name in args:
                if name not in names:
                    self._selected_names = None
                    raise EventSelectionFailed(
                        'No such %s: %s. %s' % (
                            (sevent, name, savail(sevent, names))))

                self._selected_names.append(name)

    @property
    def nselected(self):
        return len(self.get_selected_names())

    def get_selected_names(self):
        if self._selected_names is None:
            raise NoEventSelectionAvailable()

        return self._selected_names

    def get_dataset(self):
        if self._dataset is None:
            name = self.get_current_name()
            self._dataset = self.get_config().get_dataset(name)

        return self._dataset

    def set_rundir_path(self, path):
        self._rundir_path = path

    def get_rundir_path(self):
        if self._rundir_path is None:
            raise NoRundirAvailable()

        return self._rundir_path

    def get_run_info_path(self):
        return op.join(self.get_rundir_path(), 'run_info.yaml')

    def get_run_info(self):
        run_info_path = self.get_run_info_path()
        if not op.exists(run_info_path):
            info = run_info.RunInfo()
            return info
        else:
            return run_info.read_info(run_info_path)

    def set_run_info(self, info):
        run_info_path = self.get_run_info_path()
        run_info.write_info(info, run_info_path)

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
                conf = self.get_config()
                ds = self.get_dataset()
                self._problem = conf.get_problem(ds.get_event_group())

        return self._problem

    def get_history(self, subset=None):
        if subset not in self._histories:
            self._histories[subset] = \
                ModelHistory(
                    self.get_problem(),
                    nchains=self.get_optimiser().nchains,
                    path=meta.xjoin(self.get_rundir_path(), subset))

            self._histories[subset].ensure_bootstrap_misfits(
                self.get_optimiser())

        return self._histories[subset]

    def set_plot_collection_manager(self, pcm):
        self._plot_collection_manager = pcm

    def get_plot_collection_manager(self):
        if self._plot_collection_manager is None:
            raise NoPlotCollectionManagerAvailable()

        return self._plot_collection_manager

    def setup_modelling(self):
        '''Must be called before any modelling can be done.'''
        logger.debug('Setting up modelling...')
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

    def get_plots_path(self):
        try:
            return op.join(self.get_rundir_path(), 'plots')
        except NoRundirAvailable:
            return 'plots'

    def get_config_path(self):
        return self._config_path


__all__ = [
    'GrondEnvironmentError',
    'EventSelectionFailed',
    'NoCurrentEventAvailable',
    'NoRundirAvailable',
    'NoPlotCollectionManagerAvailable',
    'Environment',
]
