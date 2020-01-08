import numpy as num
from matplotlib import pyplot as plt, cm as mcm, colors as mcolors

from pyrocko.guts import Tuple, Float
from pyrocko.plot import mpl_init, mpl_margins

from grond.plot.config import PlotConfig
from grond.plot.collection import PlotItem

from .target import PhasePickTarget


class PickResidualsPlot(PlotConfig):
    '''Travel time residuals plot.'''

    name = 'pick_residuals'

    size_cm = Tuple.T(
        2, Float.T(),
        default=(15., 7.5),
        help='Width and length of the figure in [cm]')

    def make(self, environ):
        environ.setup_modelling()

        cm = environ.get_plot_collection_manager()
        mpl_init(fontsize=self.font_size)

        ds = environ.get_dataset()
        history = environ.get_history(subset='harvest')
        optimiser = environ.get_optimiser()

        cm.create_group_mpl(
            self,
            self.draw_figures(ds, history, optimiser),
            title=u'Pick Residuals',
            section='fits',
            feather_icon='activity',
            description=u'''
Travel time residuals.
''')

    def draw_figures(self, ds, history, optimiser):

        # gms = history.get_sorted_primary_misfits()[::-1]
        models = history.get_sorted_primary_models()[::-1]
        problem = history.problem

        targets = [
            t for t in problem.targets if
            isinstance(t, PhasePickTarget)]

        targets.sort(key=lambda t: t.distance_to(problem.base_source))

        tpaths = sorted(set(t.path for t in targets))

        ntargets = len(targets)
        nmodels = history.nmodels
        tts = num.zeros((nmodels, ntargets, 2))
        # distances = num.zeros((nmodels, ntargets))

        for imodel in range(nmodels):
            model = models[imodel, :]
            # source = problem.get_source(model)
            results = problem.evaluate(model, targets=targets)
            for itarget, result in enumerate(results):
                result = results[itarget]
                tts[imodel, itarget, :] = result.tobs, result.tsyn
                # distances[imodel, itarget] = source.distance_to(
                #     targets[itarget])

        ok = num.all(num.isfinite(tts), axis=2)
        ok = num.all(ok, axis=0)
        print(ok.size, ntargets)

        iorder = num.arange(nmodels)
        icolor = iorder

        cmap = mcm.ScalarMappable(
            norm=mcolors.Normalize(vmin=num.min(icolor), vmax=num.max(icolor)),
            cmap=plt.get_cmap('coolwarm'))

        imodel_to_color = []
        for imodel in range(nmodels):
            imodel_to_color.append(cmap.to_rgba(icolor[imodel]))

        for tpath in tpaths:
            mask = num.array(
                [t.path == tpath for t in targets], dtype=num.bool)

            mask &= ok

            targets_this = [t for (it, t) in enumerate(targets) if mask[it]]

            fontsize = self.font_size

            item = PlotItem(
                name='fig_%s' % tpath)

            item.attributes['targets'] = [t.string_id() for t in targets_this]
            print(item.attributes['targets'])

            fig = plt.figure(figsize=self.size_inch)

            labelpos = mpl_margins(
                fig, nw=1, nh=1, left=7., right=1., bottom=10., top=3,
                units=fontsize)

            axes = fig.add_subplot(1, 1, 1)
            labelpos(axes, 2.5, 2.0)

            for imodel in range(0, nmodels, 10):
                color = imodel_to_color[imodel]
                axes.plot(
                    num.arange(len(targets_this)),
                    # distances[imodel, mask],
                    # tts[imodel, mask, 0],
                    tts[imodel, mask, 1] - tts[imodel, mask, 0],
                    '_',
                    ms=2.0,
                    color=color)

            yield item, fig


def get_plot_classes():
    return [
        PickResidualsPlot]
