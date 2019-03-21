import numpy as num
from matplotlib import pyplot as plt
from matplotlib import cm

from pyrocko.guts import Tuple, Float, String, Int
from pyrocko.plot import mpl_init, mpl_margins
from pyrocko import gf

from grond import core, meta
from grond.plot.config import PlotConfig
from grond.plot.collection import PlotItem

from .target import PhaseRatioTarget

guts_prefix = 'grond'


def S(text):
    return ' '.join(text.split())


class FitsPhaseRatioPlot(PlotConfig):
    '''Plot showing the phase ratio fits for the best model.'''

    name = 'fits_phase_ratio'

    size_cm = Tuple.T(
        2, Float.T(),
        default=(15., 7.5),
        help='Width and length of the figure in [cm]')

    misfit_cutoff = Float.T(
        optional=True,
        help='Plot fits for models up to this misfit value.')

    color_parameter = String.T(
        default='misfit',
        help='Choice of value to color, options: "misfit" (default), '
             '"dist" or any source parameter name.')

    istride_ensemble = Int.T(
        default=10,
        help='Stride value N to choose every Nth model from the solution '
             'ensemble.')

    font_size_title = Float.T(
        default=10,
        help='Font size of title [pt]')

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        mpl_init(fontsize=self.font_size)

        environ.setup_modelling()
        ds = environ.get_dataset()
        history = environ.get_history(subset='harvest')

        scolor = {
            'misfit': S('''
                The synthetic markers are colored according to their
                respective global (non-bootstrapped) misfit value. Red
                indicates better fit, blue worse.'''),

            'dist': S('''
                The synthetic markers are colored according to their
                Mahalanobis distance from the mean solution.''')

        }.get(self.color_parameter, S('''
            The synthetic markers are colored according to source
            parameter "%s".''' % self.color_parameter))

        cm.create_group_mpl(
            self, self.draw_figures(ds, history),
            title=u'Fits Phase Ratios',
            section='fits',
            feather_icon='activity',
            description=u'''
Observed (black markers) and synthetic waveform amplitude phase ratio estimates
(spectral average ratio; colored markers) at different stations for every Nth
model in the bootstrap solution ensemble (N=%i).

%s

The frequency range used to estimate spectral averages is not shown (see config
file). Optimal solutions show good agreement between black and colored markers.
''' % (self.istride_ensemble, scolor))

    def draw_figures(self, ds, history):
        problem = history.problem

        for target in problem.targets:
            target.set_dataset(ds)

        targets = [
            t for t in problem.targets if isinstance(t, PhaseRatioTarget)]

        tpaths = sorted(set(t.path for t in targets))

        for tpath in tpaths:
            for (item, fig) in self.draw_figure(ds, history, tpath):
                yield item, fig

    def draw_figure(self, ds, history, tpath):
        problem = history.problem
        color_parameter = self.color_parameter
        misfit_cutoff = self.misfit_cutoff
        fontsize = self.font_size
        targets = [
            t for t in problem.targets
            if isinstance(t, PhaseRatioTarget) and t.path == tpath]

        gms = history.get_sorted_misfits(chain=0)[::-1]
        models = history.get_sorted_models(chain=0)[::-1]

        if misfit_cutoff is not None:
            ibest = gms < misfit_cutoff
            gms = gms[ibest]
            models = models[ibest]

        gms = gms[::self.istride_ensemble]
        models = models[::self.istride_ensemble]

        nmodels = models.shape[0]
        if color_parameter == 'dist':
            mx = num.mean(models, axis=0)
            cov = num.cov(models.T)
            mdists = core.mahalanobis_distance(models, mx, cov)
            icolor = meta.ordersort(mdists)

        elif color_parameter == 'misfit':
            iorder = num.arange(nmodels)
            icolor = iorder

        elif color_parameter in problem.parameter_names:
            ind = problem.name_to_index(color_parameter)
            icolor = problem.extract(models, ind)

        from matplotlib import colors
        cmap = cm.ScalarMappable(
            norm=colors.Normalize(vmin=num.min(icolor), vmax=num.max(icolor)),
            cmap=plt.get_cmap('coolwarm'))

        imodel_to_color = []
        for imodel in range(nmodels):
            imodel_to_color.append(cmap.to_rgba(icolor[imodel]))

        data = []
        for imodel in range(nmodels):
            model = models[imodel, :]

            # source = problem.get_source(model)
            results = problem.evaluate(model, targets=targets)

            for target, result in zip(targets, results):
                if isinstance(result, gf.SeismosizerError):
                    continue

                if not isinstance(target, PhaseRatioTarget):
                    continue

                a_obs = result.a_obs
                b_obs = result.b_obs
                a_syn = result.a_syn
                b_syn = result.b_syn

                r_obs = a_obs / (a_obs + b_obs)
                r_syn = a_syn / (a_syn + b_syn)

                data.append(('.'.join(target.codes), imodel, r_obs, r_syn))

        fontsize = self.font_size

        item = PlotItem(
            name='fig_%s' % tpath)

        item.attributes['targets'] = [
            t.string_id() for t in targets]

        fig = plt.figure(figsize=self.size_inch)

        labelpos = mpl_margins(
            fig, nw=1, nh=1, left=7., right=1., bottom=10., top=3,
            units=fontsize)

        axes = fig.add_subplot(1, 1, 1)
        labelpos(axes, 2.5, 2.0)

        labels = sorted(set(x[0] for x in data))

        ntargets = len(labels)
        string_id_to_itarget = dict((x, i) for (i, x) in enumerate(labels))

        itargets = num.array([string_id_to_itarget[x[0]] for x in data])

        imodels = num.array([x[1] for x in data], dtype=num.int).T
        r_obs, r_syn = num.array([x[2:] for x in data]).T

        r_obs_median = num.zeros(ntargets)
        for itarget in range(ntargets):
            r_obs_median[itarget] = num.median(r_obs[itargets == itarget])

        iorder = meta.ordersort(r_obs_median)

        for imodel in range(nmodels):
            mask = imodels == imodel
            axes.plot(
                iorder[itargets[mask]], r_obs[mask], '_',
                ms=20.,
                zorder=-10,
                alpha=0.5,
                color='black')
            axes.plot(
                iorder[itargets[mask]], r_syn[mask], '_',
                ms=10.,
                alpha=0.5,
                color=imodel_to_color[imodel])

        axes.set_yscale('log')
        axes.set_ylabel('Ratio')

        axes.set_xticks(
            iorder[num.arange(ntargets)])
        axes.set_xticklabels(labels, rotation='vertical')

        fig.suptitle(tpath, fontsize=self.font_size_title)

        yield item, fig


def get_plot_classes():
    return [
        FitsPhaseRatioPlot,
    ]
