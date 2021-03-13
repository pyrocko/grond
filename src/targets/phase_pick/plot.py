import math

from matplotlib import pyplot as plt, colors as mcolors
import numpy as num

from pyrocko import util
from pyrocko.guts import Tuple, Float
from pyrocko.plot import mpl_init, mpl_margins

from grond import meta
from grond.plot.config import PlotConfig
from grond.plot.collection import PlotItem

from .target import PhasePickTarget
from ..plot import StationDistributionPlot

km = 1000.
d2r = math.pi / 180.


class PickResidualsPlot(PlotConfig):
    '''Travel time residuals plot.'''

    name = 'pick_residuals'

    size_cm = Tuple.T(
        2, Float.T(),
        default=(15.0, 10.0),
        help='Width and height of the figure in [cm]')

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
            feather_icon='watch',
            description=u'''
Travel time residuals.

The difference between observed and synthetic phase arrival times is shown
for a subset of the models in the solution ensemble. Models with good global
performance are shown in red, less good ones in blue. Vertical black bars mark
the classical best model.
''')

    def draw_figures(self, ds, history, optimiser):
        show_mean_residuals = False

        # gms = history.get_sorted_primary_misfits()[::-1]
        models = history.get_sorted_primary_models()[::-1]
        problem = history.problem

        targets = [
            t for t in problem.targets if
            isinstance(t, PhasePickTarget)]

        # targets.sort(key=lambda t: t.distance_to(problem.base_source))
        # tpaths = sorted(set(t.path for t in targets))

        ntargets = len(targets)
        nmodels = history.nmodels
        tts = num.zeros((nmodels, ntargets, 2))
        distances = num.zeros((nmodels, ntargets))
        azimuths = num.zeros((nmodels, ntargets))

        for imodel in range(nmodels):
            model = models[imodel, :]
            source = problem.get_source(model)
            results = problem.evaluate(model, targets=targets)
            for itarget, result in enumerate(results):
                result = results[itarget]
                tts[imodel, itarget, :] = result.tobs, result.tsyn
                distances[imodel, itarget] = source.distance_to(
                    targets[itarget])

                azimuths[imodel, itarget] = source.azibazi_to(
                    targets[itarget])[0]

        ok = num.all(num.isfinite(tts), axis=2)
        ok = num.all(ok, axis=0)

        targets = [
            target for (itarget, target) in enumerate(targets) if ok[itarget]]
        tts = tts[:, ok, :]
        distances = distances[:, ok]
        azimuths = azimuths[:, ok]
        residuals = tts[:, :, 0] - tts[:, :, 1]
        residual_amax = num.max(num.abs(residuals))
        mean_residuals = num.mean(residuals, axis=0)
        mean_distances = num.mean(distances, axis=0)

        isort = num.array(
            [x[2] for x in sorted(
                zip(targets, mean_distances, range(len(targets))),
                key=lambda x: (x[0].path, x[1]))], dtype=num.int)

        distances = distances[:, isort]
        distance_min = num.min(distances)
        distance_max = num.max(distances)
        azimuths = azimuths[:, isort]
        targets = [targets[i] for i in isort]
        ntargets = len(targets)
        residuals = residuals[:, isort]
        mean_residuals = mean_residuals[isort]

        icolor = num.arange(nmodels)

        norm = mcolors.Normalize(vmin=num.min(icolor), vmax=num.max(icolor))
        cmap = plt.get_cmap('coolwarm')

        napprox = max(
            1, int(round((self.size_points[1] / self.font_size - 7))))

        npages = (ntargets-1) // napprox + 1
        ntargets_page = (ntargets-1) // npages + 1

        for ipage in range(npages):

            ilo, ihi = ipage*ntargets_page, (ipage+1)*ntargets_page
            ihi = min(len(targets), ihi)

            targets_page = targets[ilo:ihi]
            item = PlotItem(
                name='fig_%02i' % ipage)

            item.attributes['targets'] = [t.string_id() for t in targets_page]

            fig = plt.figure(figsize=self.size_inch)

            labelpos = mpl_margins(
                fig, nw=2, nh=1, left=12., right=1., bottom=5., top=1.,
                wspace=1., units=self.font_size)

            axes1 = fig.add_subplot(1, 3, 1)
            axes2 = fig.add_subplot(1, 3, 2)
            axes3 = fig.add_subplot(1, 3, 3)

            for axes in (axes1, axes2, axes3):
                axes.set_ylim(-0.5, len(targets_page)-0.5)
                axes.invert_yaxis()
                labelpos(axes, 2.5, 2.0)

            axes1.axvline(0.0, color='black')

            labels = []
            lastpath = None
            for ipos, t in enumerate(targets_page):
                scodes = '.'.join(t.codes)
                if lastpath is None or lastpath != t.path:
                    labels.append(t.path + '.' + scodes)
                    if lastpath is not None:
                        for axes in (axes1, axes2, axes3):
                            axes.axhline(ipos-0.5, color='black')
                else:
                    labels.append(scodes)

                lastpath = t.path

            for ipos, itarget in enumerate(range(ilo, ihi)):
                axes1.plot(
                    residuals[-1, itarget],
                    ipos,
                    '|',
                    ms=self.font_size,
                    color='black')

                if show_mean_residuals:
                    axes1.plot(
                        mean_residuals[itarget],
                        ipos,
                        's',
                        ms=self.font_size,
                        color='none',
                        mec='black')

                axes1.scatter(
                    residuals[::10, itarget],
                    util.num_full(icolor[::10].size, ipos, dtype=num.float),
                    c=icolor[::10],
                    cmap=cmap,
                    norm=norm,
                    alpha=0.5)

                axes2.scatter(
                    distances[::10, itarget]/km,
                    util.num_full(icolor[::10].size, ipos, dtype=num.float),
                    c=icolor[::10],
                    cmap=cmap,
                    norm=norm,
                    alpha=0.5)

                axes3.scatter(
                    azimuths[::10, itarget],
                    util.num_full(icolor[::10].size, ipos, dtype=num.float),
                    c=icolor[::10],
                    cmap=cmap,
                    norm=norm,
                    alpha=0.5)

            axes1.set_yticks(num.arange(len(labels)))
            axes1.set_yticklabels(labels)
            axes1.set_xlabel('$T_{obs} - T_{syn}$ [s]')
            axes1.set_xlim(-residual_amax, residual_amax)

            axes2.set_xlabel('Distance [km]')
            axes2.set_xlim(distance_min/km, distance_max/km)

            axes3.set_xlabel('Azimuth [deg]')
            axes3.set_xlim(-180., 180.)

            axes2.get_yaxis().set_ticks([])
            axes2.get_yaxis().set_ticks([])

            axes3.get_yaxis().set_ticks([])
            axes3.get_yaxis().set_ticks([])

            yield item, fig


class PickResidualsStationDistribution(StationDistributionPlot):
    '''
    Plot showing residuals of seismic phase arrivals by azimuth and distance.
    '''

    name = 'pick_residuals_stations'

    def make(self, environ):
        environ.setup_modelling()

        cm = environ.get_plot_collection_manager()
        mpl_init(fontsize=self.font_size)

        problem = environ.get_problem()
        dataset = environ.get_dataset()
        history = environ.get_history(subset='harvest')

        cm.create_group_mpl(
            self,
            self.draw_figures(problem, dataset, history),
            title=u'Pick Residuals by Location',
            section='fits',
            feather_icon='target',
            description=u'''
Plot showing residuals for seismic phase arrivals by azimuth and distance.

Station locations in dependence of distance and azimuth are shown. The center
of the plot corresponds to the origin of the search space. Optimized source
locations are shown as small black dots (subset of the solution ensemble).
''')

    def draw_figures(self, problem, dataset, history):

        target_index = {}
        i = 0
        for target in problem.targets:
            target_index[target] = i, i+target.nmisfits
            i += target.nmisfits

        misfits = history.misfits[history.get_sorted_misfits_idx(), ...]
        gcms = problem.combine_misfits(
            misfits[:1, :, :], get_contributions=True)[0, :]
        models = history.get_sorted_primary_models()[::-1]

        origin = problem.base_source

        targets = [
            t for t in problem.targets if isinstance(t, PhasePickTarget)]

        ntargets = len(targets)
        nmodels = history.nmodels
        tts = num.zeros((nmodels, ntargets, 2))

        sdata = []
        for imodel in range(nmodels):
            model = models[imodel, :]
            source = problem.get_source(model)
            sdata.append(
                (origin.distance_to(source), origin.azibazi_to(source)[0]))

            results = problem.evaluate(model, targets=targets)
            for itarget, result in enumerate(results):
                result = results[itarget]
                tts[imodel, itarget, :] = result.tobs, result.tsyn

        ok = num.all(num.isfinite(tts), axis=2)
        ok = num.all(ok, axis=0)

        targets_ok = [
            target for (itarget, target) in enumerate(targets) if ok[itarget]]
        tts = tts[:, ok, :]
        residuals = tts[:, :, 0] - tts[:, :, 1]
        mean_residuals = num.mean(residuals, axis=0)

        rlo, rhi = num.percentile(mean_residuals, [10., 90.])
        residual_amax = max(abs(rlo), abs(rhi))

        target_to_residual = dict(
            (target, residual)
            for (target, residual)
            in zip(targets_ok, mean_residuals))

        cg_to_targets = meta.gather(targets, lambda t: (t.path, ))

        cgs = sorted(cg_to_targets.keys())

        for cg in cgs:
            cg_str = '.'.join(cg)

            targets = cg_to_targets[cg]
            if len(targets) == 0:
                continue

            assert all(target_index[target][0] == target_index[target][1] - 1
                       for target in targets)

            itargets = num.array(
                [target_index[target][0] for target in targets])

            labels = ['.'.join(x for x in t.codes[:3] if x) for t in targets]

            azimuths = num.array(
                [origin.azibazi_to(t)[0] for t in targets])
            distances = num.array(
                [t.distance_to(origin) for t in targets])
            residuals = num.array(
                [target_to_residual.get(t, num.nan) for t in targets])

            item = PlotItem(
                name='picks_contributions_%s' % cg_str,
                title=u'Pick residuals and contributions (%s)' % cg_str,
                description=u'\n\nMarkers are scaled according to their '
                            u'misfit contribution for the globally best '
                            u'source model.')

            fig, axes, legend = self.plot_station_distribution(
                azimuths, distances,
                gcms[itargets],
                labels,
                colors=residuals,
                cnorm=(-residual_amax, residual_amax),
                cmap='RdYlBu',
                clabel='$T_{obs} - T_{syn}$ [s]',
                scatter_kwargs=dict(alpha=1.0),
                legend_title='Contribution')

            sources = [
                    problem.get_source(x) for x in models[::10]]

            azimuths = num.array(
                [origin.azibazi_to(s)[0] for s in sources])
            distances = num.array(
                [s.distance_to(origin) for s in sources])

            axes.plot(
                azimuths*d2r, distances, 'o', ms=2.0, color='black', alpha=0.3)

            yield (item, fig)


def get_plot_classes():
    return [
        PickResidualsPlot,
        PickResidualsStationDistribution]
