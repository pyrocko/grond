import os.path as op
from grond.environment import Environment
from grond.plot.collection import PlotCollectionManager
from grond.plot.config import PlotConfigCollection


def get_plot_names(args):
    env = Environment(*args)
    plot_classes = env.get_plots()
    return [plot_class.name for plot_class in plot_classes]


def get_plot_config_collection(args):
    env = Environment(*args)
    plot_classes = env.get_plots()
    collection = PlotConfigCollection()

    for plot_class in plot_classes:
        plot_config = plot_class()
        collection.plot_configs.append(plot_config)

    return collection


def make_plots(plots, args, plots_path=None):
    env = Environment(*args)
    if isinstance(plots, PlotConfigCollection):
        plots = plots.plot_configs

    else:
        plot_classes = env.get_plots()
        for plot_class in plot_classes:
            plots = [
                plot_class()
                for plot_class in plot_classes
                if plot_class.name in plots]

    if plots_path is None:
        plots_path = op.join(env.get_rundir_path(), 'plots')

    manager = PlotCollectionManager(plots_path)
    env.set_plot_collection_manager(manager)

    for plot in plots:
        plot.make(env)


def make_movie(dirname, xpar_name, ypar_name, movie_filename):
    env = Environment(dirname)
    optimizer = env.get_optimizer()
    problem = env.get_problem()
    history = env.get_history()
    movie_maker = optimizer.get_movie_maker(
        problem, history, xpar_name, ypar_name, movie_filename)

    movie_maker.render()


__all__ = [
    'get_plot_names',
    'get_plot_config_collection',
    'make_plots',
    'make_movie']
