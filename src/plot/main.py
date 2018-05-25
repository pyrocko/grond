import logging
from grond.meta import GrondError
from grond.environment import Environment, GrondEnvironmentError
from grond.plot.collection import PlotCollectionManager
from grond.plot.config import PlotConfigCollection

logger = logging.getLogger('grond.plots')


def get_plot_names(env):
    plot_classes = env.get_plot_classes()
    return [plot_class.name for plot_class in plot_classes]


def get_plot_config_collection(env, plot_names=None):

    plots = []
    for plot_class in env.get_plot_classes():
        if plot_names is None or plot_class.name in plot_names:
            plots.append(plot_class())

    if plot_names is not None:
        if set(plot_names) - set([p.name for p in plots]):
            logger.warning(
                'Plots %s not available!'
                % ', '.join(set(plot_names) - set([p.name for p in plots])))

    collection = PlotConfigCollection(plot_configs=plots)

    return collection


def make_plots(
        env,
        plot_config_collection=None,
        plot_names=None,
        plots_path=None):

    if plot_config_collection is None:
        plot_config_collection = get_plot_config_collection(env, plot_names)

    if plots_path is None:
        plots_path = env.get_plots_path()

    plots = plot_config_collection.plot_configs
    manager = PlotCollectionManager(plots_path)
    env.set_plot_collection_manager(manager)

    for plot in plots:
        try:
            plot.make(env)
        except (GrondEnvironmentError, GrondError) as e:
            logger.warning('cannot create plot %s: %s' % (
                plot.name, str(e)))


def make_movie(dirname, xpar_name, ypar_name, movie_filename):
    env = Environment([dirname])
    optimiser = env.get_optimiser()
    problem = env.get_problem()
    history = env.get_history()
    movie_maker = optimiser.get_movie_maker(
        problem, history, xpar_name, ypar_name, movie_filename)

    movie_maker.render()


__all__ = [
    'get_plot_names',
    'get_plot_config_collection',
    'make_plots',
    'make_movie']
