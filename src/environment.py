import time
import logging
import shutil
import os

from grond.config import read_config, write_config
from grond import meta, run_info
from grond.problems.base import load_optimiser_info, load_problem_info, \
    ModelHistory

op = os.path

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

    def __init__(self, args=None, config=None, event_names=None):

        self._current_event_name = None
        self._selected_event_names = None
        self._config = None
        self._plot_collection_manager = None

        if isinstance(args, str):
            args = [args]

        if not args and not config:
            raise GrondEnvironmentError('missing arguments')

        if config and event_names:
            self._config_path = None
            self._rundir_path = None
            self._config = config

            if isinstance(event_names, str):
                event_names = [event_names]
            self.set_selected_event_names(event_names)

        elif op.isdir(args[0]):
            self._rundir_path = args[0]
            self._config_path = op.join(self._rundir_path, 'config.yaml')

        else:
            self._rundir_path = None
            self._config_path = args[0]
            self.set_selected_event_names(args[1:])

        self.reset()

    @classmethod
    def discover(cls, rundir):
        running_fn = op.join(rundir, '.running')
        while op.exists(running_fn):
            try:
                cls.verify_rundir(rundir)
                return cls([rundir])
            except GrondEnvironmentError:
                time.sleep(.25)
        raise GrondEnvironmentError('could not discover rundir')

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

    def copy(self, destination, force=False):
        ''' Copy the environment and return it '''
        files = [
            'config.yaml',
            'problem.yaml',
            'optimiser.yaml',
            'misfits',
            'models',
            'choices',
            'chains'
            ]

        if op.exists(destination) and not force:
            raise OSError('Directory %s already exists' % destination)

        destination = op.abspath(destination)
        os.makedirs(destination, exist_ok=True)

        for file in files:
            src = op.join(self._rundir_path, file)
            dest = op.join(destination, file)

            if not op.isfile(src):
                logger.debug('Cannot find file %s', src)
                continue
            logger.debug('Copying %s to %s', src, dest)

            shutil.copy(src, dest)

        cls = self.__class__
        return cls(destination)

    def reset(self):
        self._histories = {}
        self._dataset = None
        self._optimiser = None
        self._problem = None

    def get_config(self):
        if self._config is None:
            self._config = read_config(self._config_path)

        return self._config

    def write_config(self):
        write_config(self.get_config(), self.get_config_path())

    def get_available_event_names(self):
        return self.get_config().get_event_names()

    def set_current_event_name(self, event_name):
        self._current_event_name = event_name
        self.reset()

    def get_current_event_name(self):
        if self._current_event_name is None:
            try:
                self.get_rundir_path()
                self._current_event_name = self.get_problem().base_source.name
            except NoRundirAvailable:
                try:
                    event_names = self.get_selected_event_names()
                    if len(event_names) == 1:
                        self._current_event_name = event_names[0]
                    else:
                        raise NoCurrentEventAvailable()

                except NoEventSelectionAvailable:
                    raise NoCurrentEventAvailable()

        return self._current_event_name

    def set_selected_event_names(self, args):
        event_names = self.get_available_event_names()
        if len(args) == 0:
            if len(event_names) == 1:
                self._selected_event_names = event_names
            else:
                if not event_names:
                    raise EventSelectionFailed(
                        'No event file found, check your config!')
                raise EventSelectionFailed(
                    'Ambiguous event selection. Select from available events:'
                    '\n    %s\n  or \'all\' to use all available events'
                    % '\n    '.join(event_names))

        elif len(args) == 1 and args[0] == 'all':
            self._selected_event_names = event_names

        else:
            self._selected_event_names = []
            for event_name in args:
                if event_name not in event_names:
                    self._selected_event_names = None
                    raise EventSelectionFailed(
                        'No such event: %s' % event_name)

                self._selected_event_names.append(event_name)

    @property
    def nevents_selected(self):
        return len(self.get_selected_event_names())

    def get_selected_event_names(self):
        if self._selected_event_names is None:
            raise NoEventSelectionAvailable()

        return self._selected_event_names

    def get_dataset(self):
        if self._dataset is None:
            event_name = self.get_current_event_name()
            self._dataset = self.get_config().get_dataset(event_name)

        return self._dataset

    def set_rundir_path(self, path):
        self._rundir_path = path

    def get_rundir_path(self):
        if self._rundir_path is None:
            raise NoRundirAvailable()

        return self._rundir_path

    def have_rundir(self):
        return self._rundir_path is not None

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
                self._problem = \
                    self.get_config().get_problem(
                        self.get_dataset().get_event())

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

    def is_running(self):
        return op.exists(self.get_rundir_path, '.running')


__all__ = [
    'GrondEnvironmentError',
    'EventSelectionFailed',
    'NoCurrentEventAvailable',
    'NoRundirAvailable',
    'NoPlotCollectionManagerAvailable',
    'Environment',
]
