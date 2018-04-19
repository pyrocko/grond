
def discover(mediator):
    plots = set() 
    plots.update(mediator.get_problem())
    return plots


def plot_result(dirname, plotnames_want,
                save=False, save_path=None, formats=('pdf',), dpi=None):

    if isinstance(formats, str):
        formats = formats.split(',')

    plotnames_want = set(plotnames_want)
    plotnames_avail = set(plot_dispatch.keys())

    if save_path is None:
        plot_dirname = op.join(dirname, 'plots')
    else:
        plot_dirname = save_path
        save = True

    unavailable = plotnames_want - plotnames_avail
    if unavailable:
        raise core.GrondError(
            'unavailable plotname: %s' % ', '.join(unavailable))

    fontsize = 10.0

    mpl_init(fontsize=fontsize)
    fns = defaultdict(list)
    config = None

    optimizer_fn = op.join(dirname, 'optimizer.yaml')
    optimizer = guts.load(filename=optimizer_fn)

    if 3 != len({'bootstrap', 'sequence', 'contributions'} - plotnames_want):
        problem = load_problem_info(dirname)
        history = ModelHistory(problem, path=dirname)

        for plotname in ['bootstrap', 'sequence', 'contributions']:
            if plotname in plotnames_want:
                figs = plot_dispatch[plotname](history, optimizer, plt)
                if save:
                    fns[plotname].extend(
                        save_figs(figs, plot_dirname, plotname, formats, dpi))

                    for fig in figs:
                        plt.close(fig)

    if 8 != len({
            'fits',
            'fits_statics',
            'fits_ensemble',
            'jointpar',
            'histogram',
            'hudson',
            'solution',
            'location'} - plotnames_want):

        problem = load_problem_info(dirname)
        history = ModelHistory(problem, path=meta.xjoin(dirname, 'harvest'))

        for plotname in ['fits', 'fits_ensemble', 'fits_statics']:
            if plotname in plotnames_want:
                event_name = problem.base_source.name

                if config is None:
                    config = guts.load(
                        filename=op.join(dirname, 'config.yaml'))
                    config.set_basepath(dirname)
                    config.setup_modelling_environment(problem)

                ds = config.get_dataset(event_name)
                figs = plot_dispatch[plotname](ds, history, optimizer, plt)
                if save:
                    fns[plotname].extend(
                        save_figs(figs, plot_dirname, plotname, formats, dpi))

                    for fig in figs:
                        plt.close(fig)

        for plotname in [
                'jointpar',
                'histogram',
                'hudson',
                'solution',
                'location']:

            if plotname in plotnames_want:
                figs = plot_dispatch[plotname](history, optimizer, plt)
                if save:
                    fns[plotname].extend(
                        save_figs(figs, plot_dirname, plotname, formats, dpi))

                    for fig in figs:
                        plt.close(fig)

    if not save:
        plt.show()

    return fns


def make_movie(dirname, xpar_name, ypar_name, movie_filename):
    optimizer_fn = op.join(dirname, 'optimizer.yaml')
    optimizer = guts.load(filename=optimizer_fn)
    problem = load_problem_info(dirname)
    history = ModelHistory(problem, path=dirname)
    movie_maker = optimizer.get_movie_maker(
        problem, history, xpar_name, ypar_name, movie_filename)

    movie_maker.render()


def draw_target_check_figures(sources, target, results):
    try:
        plotter_class = target.get_plotter_class()
        return plotter_class.draw_check_figures(sources, target, results)

    except plotter.NoPlotterClassAvailable:
        return []


__all__ = []
