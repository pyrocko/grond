import logging
import math
from collections import defaultdict

import numpy as num

from matplotlib import pyplot as plt
from matplotlib import cm, patches

from pyrocko import plot, gf, trace
from pyrocko.plot import mpl_init, mpl_color
from pyrocko.guts import Tuple, Float, Int, String

from grond import core, meta
from .target import WaveformMisfitResult, WaveformMisfitTarget
from ..plot import StationDistributionPlot

from grond.plot.config import PlotConfig
from grond.plot.common import light
from grond.plot.collection import PlotItem

guts_prefix = 'grond'
km = 1e3
d2r = math.pi / 180.

logger = logging.getLogger('targets.waveform.plot')


def make_norm_trace(a, b, exponent):
    tmin = max(a.tmin, b.tmin)
    tmax = min(a.tmax, b.tmax)
    c = a.chop(tmin, tmax, inplace=False)
    bc = b.chop(tmin, tmax, inplace=False)
    c.set_ydata(num.abs(c.get_ydata() - bc.get_ydata())**exponent)
    return c


def amp_spec_max(spec_trs, key):
    amaxs = {}
    for spec_tr in spec_trs:
        amax = num.max(num.abs(spec_tr.ydata))
        k = key(spec_tr)
        if k not in amaxs:
            amaxs[k] = amax
        else:
            amaxs[k] = max(amaxs[k], amax)

    return amaxs


def plot_trace(axes, tr, **kwargs):
    return axes.plot(tr.get_xdata(), tr.get_ydata(), **kwargs)


def plot_taper(axes, t, taper, **kwargs):
    y = num.ones(t.size) * 0.9
    taper(y, t[0], t[1] - t[0])
    y2 = num.concatenate((y, -y[::-1]))
    t2 = num.concatenate((t, t[::-1]))
    axes.fill(t2, y2, **kwargs)


def plot_dtrace(axes, tr, space, mi, ma, **kwargs):
    t = tr.get_xdata()
    y = tr.get_ydata()
    y2 = (num.concatenate((y, num.zeros(y.size))) - mi) / \
        (ma - mi) * space - (1.0 + space)
    t2 = num.concatenate((t, t[::-1]))
    axes.fill(
        t2, y2,
        clip_on=False,
        **kwargs)


def plot_spectrum(
        axes, spec_syn, spec_obs, fmin, fmax, space, mi, ma,
        syn_color='red', obs_color='black',
        syn_lw=1.5, obs_lw=1.0, color_vline='gray', fontsize=9.):

    fpad = (fmax - fmin) / 6.

    for spec, color, lw in [
            (spec_syn, syn_color, syn_lw),
            (spec_obs, obs_color, obs_lw)]:

        f = spec.get_xdata()
        mask = num.logical_and(fmin - fpad <= f, f <= fmax + fpad)

        f = f[mask]
        y = num.abs(spec.get_ydata())[mask]

        y2 = (num.concatenate((y, num.zeros(y.size))) - mi) / \
            (ma - mi) * space - (1.0 + space)
        f2 = num.concatenate((f, f[::-1]))
        axes2 = axes.twiny()
        axes2.set_axis_off()

        axes2.set_xlim(fmin - fpad * 5, fmax + fpad * 5)

        axes2.plot(f2, y2, clip_on=False, color=color, lw=lw)
        axes2.fill(f2, y2, alpha=0.1, clip_on=False, color=color)

    axes2.plot([fmin, fmin], [-1.0 - space, -1.0], color=color_vline)
    axes2.plot([fmax, fmax], [-1.0 - space, -1.0], color=color_vline)

    for (text, fx, ha) in [
            ('%.3g Hz' % fmin, fmin, 'right'),
            ('%.3g Hz' % fmax, fmax, 'left')]:

        axes2.annotate(
            text,
            xy=(fx, -1.0),
            xycoords='data',
            xytext=(
                fontsize * 0.4 * [-1, 1][ha == 'left'],
                -fontsize * 0.2),
            textcoords='offset points',
            ha=ha,
            va='top',
            color=color_vline,
            fontsize=fontsize)


def plot_dtrace_vline(axes, t, space, **kwargs):
    axes.plot([t, t], [-1.0 - space, -1.0], **kwargs)


def layout(source, targets, nxmax, nymax):
    if nxmax == 1 and nymax == 1:
        nx = len(targets)
        ny = 1
        nxx = 1
        nyy = 1
        frame_to_target = {}
        for ix, target in enumerate(sorted(targets, key=lambda t: t.path)):
            frame_to_target[0, ix] = target

        return frame_to_target, nx, ny, nxx, nyy

    nframes = len(targets)

    nx = int(math.ceil(math.sqrt(nframes)))
    ny = (nframes - 1) // nx + 1

    nxx = (nx - 1) // nxmax + 1
    nyy = (ny - 1) // nymax + 1

    # nz = nxx * nyy

    xs = num.arange(nx) // ((max(2, nx) - 1.0) / 2.)
    ys = num.arange(ny) // ((max(2, ny) - 1.0) / 2.)

    xs -= num.mean(xs)
    ys -= num.mean(ys)

    fxs = num.tile(xs, ny)
    fys = num.repeat(ys, nx)

    data = []

    for target in targets:
        azi = source.azibazi_to(target)[0]
        dist = source.distance_to(target)
        x = dist * num.sin(num.deg2rad(azi))
        y = dist * num.cos(num.deg2rad(azi))
        data.append((x, y, dist))

    gxs, gys, dists = num.array(data, dtype=num.float).T

    iorder = num.argsort(dists)

    gxs = gxs[iorder]
    gys = gys[iorder]
    targets_sorted = [targets[ii] for ii in iorder]

    gxs -= num.mean(gxs)
    gys -= num.mean(gys)

    gmax = max(num.max(num.abs(gys)), num.max(num.abs(gxs)))
    if gmax == 0.:
        gmax = 1.

    gxs /= gmax
    gys /= gmax

    dists = num.sqrt(
        (fxs[num.newaxis, :] - gxs[:, num.newaxis])**2 +
        (fys[num.newaxis, :] - gys[:, num.newaxis])**2)

    distmax = num.max(dists)

    availmask = num.ones(dists.shape[1], dtype=num.bool)
    frame_to_target = {}
    for itarget, target in enumerate(targets_sorted):
        iframe = num.argmin(
            num.where(availmask, dists[itarget], distmax + 1.))
        availmask[iframe] = False
        iy, ix = num.unravel_index(iframe, (ny, nx))
        frame_to_target[iy, ix] = target

    return frame_to_target, nx, ny, nxx, nyy


class CheckWaveformsPlot(PlotConfig):
    ''' Plot for checking the waveforms fit with a number of synthetics '''
    name = 'check_waveform'

    size_cm = Tuple.T(
        2, Float.T(),
        default=(9., 7.5),
        help='width and length of the figure in cm')
    n_random_synthetics = Int.T(
        default=10,
        help='Number of Synthetics to generate')

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        mpl_init(fontsize=self.font_size)

        environ.setup_modelling()

        problem = environ.get_problem()
        results_list = []
        sources = []
        if self.n_random_synthetics == 0:
            x = problem.preconstrain(problem.get_reference_model())
            sources.append(problem.base_source)
            results = problem.evaluate(x)
            results_list.append(results)

        else:
            for _ in range(self.n_random_synthetics):
                x = problem.get_random_model()
                sources.append(problem.get_source(x))
                results = problem.evaluate(x)
                results_list.append(results)

        cm.create_group_mpl(self, self.draw_figures(
            sources, problem.targets, results_list),
            title=u'Waveform Check',
            section='checks',
            feather_icon='activity',
            description=u'''
Plot to judge waveform time window settings and source model parameter ranges.

For each waveform target, observed and synthetic waveforms are shown. For the
latter, models are randomly drawn from the configured parameter search space.

The top panel shows the observed waveform; filtered (faint gray), and filtered
and tapered (black). The colored outline around the observed trace shows the
taper position for each drawn model in a different color. The middle panel
shows the filtered synthetic waveforms of the drawn models and the bottom plot
shows the corresponding filtered and tapered synthetic waveforms. The colors of
taper and synthetic traces are consistent for each random model. The given time
is relative to the reference event origin time.
''')

    def draw_figures(self, sources, targets, results_list):
        results_list = list(zip(*results_list))
        for itarget, target, results in zip(
                range(len(targets)), targets, results_list):

            if isinstance(target, WaveformMisfitTarget) and results:
                item = PlotItem(name='t%i' % itarget)
                item.attributes['targets'] = [target.string_id()]
                fig = self.draw_figure(sources, target, results)
                if fig is not None:
                    yield item, fig

    def draw_figure(self, sources, target, results):
        t0_mean = num.mean([s.time for s in sources])

        # distances = [
        #    s.distance_to(target) for s in sources]

        # distance_min = num.min(distances)
        # distance_max = num.max(distances)

        yabsmaxs = []

        for result in results:
            if isinstance(result, WaveformMisfitResult):
                yabsmaxs.append(
                    num.max(num.abs(
                        result.filtered_obs.get_ydata())))

        if yabsmaxs:
            yabsmax = max(yabsmaxs) or 1.0
        else:
            yabsmax = None

        fontsize = self.font_size

        fig = plt.figure(figsize=self.size_inch)

        labelpos = plot.mpl_margins(
            fig, nw=1, nh=1, w=1., h=5.,
            units=fontsize)

        axes = fig.add_subplot(1, 1, 1)

        labelpos(axes, 2.5, 2.0)
        axes.set_frame_on(False)
        axes.set_ylim(1., 4.)
        axes.get_yaxis().set_visible(False)
        axes.set_title('%s' % target.string_id())
        axes.set_xlabel('Time [s]')
        ii = 0

        for source, result in zip(sources, results):
            if not isinstance(result, WaveformMisfitResult):
                continue

            if result.tobs_shift != 0.0:
                t0 = result.tsyn_pick
            else:
                t0 = t0_mean

            t = result.filtered_obs.get_xdata()
            ydata = result.filtered_obs.get_ydata() / yabsmax
            axes.plot(
                t-t0, ydata*0.40 + 3.5, color='black', lw=1.0)

            color = plot.mpl_graph_color(ii)

            t = result.filtered_syn.get_xdata()
            ydata = result.filtered_syn.get_ydata()
            ydata = ydata / (num.max(num.abs(ydata)) or 1.0)

            axes.plot(t-t0, ydata*0.47 + 2.5, color=color, alpha=0.5, lw=1.0)

            t = result.processed_syn.get_xdata()
            ydata = result.processed_syn.get_ydata()
            ydata = ydata / (num.max(num.abs(ydata)) or 1.0)

            axes.plot(t-t0, ydata*0.47 + 1.5, color=color, alpha=0.5, lw=1.0)
            if result.tobs_shift != 0.0:
                axes.axvline(
                    result.tsyn_pick - t0,
                    color=(0.7, 0.7, 0.7),
                    zorder=2)

            t = result.processed_syn.get_xdata()
            taper = result.taper

            y = num.ones(t.size) * 0.9
            taper(y, t[0], t[1] - t[0])
            y2 = num.concatenate((y, -y[::-1]))
            t2 = num.concatenate((t, t[::-1]))
            axes.plot(t2-t0, y2 * 0.47 + 3.5, color=color, alpha=0.2, lw=1.0)
            ii += 1

        return fig


class FitsWaveformEnsemblePlot(PlotConfig):
    ''' Plot showing all waveform fits for the ensemble of solutions'''

    name = 'fits_waveform_ensemble'
    size_cm = Tuple.T(
        2, Float.T(),
        default=(9., 5.),
        help='width and length of the figure in cm')
    nx = Int.T(
        default=1,
        help='horizontal number of subplots on every page')
    ny = Int.T(
        default=1,
        help='vertical number of subplots on every page')
    misfit_cutoff = Float.T(
        optional=True,
        help='Plot fits for models up to this misfit value')
    color_parameter = String.T(
        default='misfit',
        help='Choice of value to color, options: dist and misfit')
    font_size = Float.T(
        default=8,
        help='Font Size of all fonts, except title')
    font_size_title = Float.T(
        default=10,
        help='Font size of title')

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        mpl_init(fontsize=self.font_size)
        environ.setup_modelling()
        ds = environ.get_dataset()
        history = environ.get_history(subset='harvest')
        optimiser = environ.get_optimiser()

        cm.create_group_mpl(
            self,
            self.draw_figures(ds, history, optimiser),
            title=u'Waveform fits for the ensemble',
            section='fits',
            feather_icon='activity',
            description=u'''
Plot showing waveform (attribute) fits for the ensemble of solutions.

Waveform fits for every nth model in the ensemble of bootstrap solutions.
Depending on the target configuration different types of comparisons are
possible: (i) time domain waveform differences, (ii) amplitude spectra, (iii)
envelopes, (iv) cross correlation functions. Each waveform plot gives a number
of details:

1) Target information (left side, from top to bottom) gives station name with
component, distance to source, azimuth of station with respect to source,
target weight, target misfit and starting time of the waveform relative to the
origin time.

2) The background gray area shows the applied taper function.

3) The waveforms shown are: the restituted and filtered observed trace without
tapering (light grey) and the same trace with tapering and processing (dark
gray), the synthetic trace (light red) and the filtered, tapered and (if
enabled) shifted and processed synthetic trace (colored). The colors of the
synthetic traces indicate how well the corresponding models fit in the global
weighting scheme (when all bootstrap weights are equal), from better fit (red)
to worse fit (blue). The amplitudes of the traces are scaled according to the
target weight (small weight, small amplitude) and normed relative to the
maximum amplitude of the targets of the corresponding normalisation family.

4) The bottom panel shows, depending on the type of comparison, sample-wise
residuals for time domain comparisons (red filled), spectra of observed and
synthetic traces for amplitude spectrum comparisons, or cross correlation
traces.''')

    def draw_figures(self, ds, history, optimiser):

        color_parameter = self.color_parameter
        misfit_cutoff = self.misfit_cutoff
        fontsize = self.font_size
        fontsize_title = self.font_size_title

        nxmax = self.nx
        nymax = self.ny

        problem = history.problem

        for target in problem.targets:
            target.set_dataset(ds)

        target_index = {}
        i = 0
        for target in problem.targets:
            target_index[target] = i, i+target.nmisfits
            i += target.nmisfits

        gms = history.get_sorted_primary_misfits()[::-1]
        models = history.get_sorted_primary_models()[::-1]

        if misfit_cutoff is not None:
            ibest = gms < misfit_cutoff
            gms = gms[ibest]
            models = models[ibest]

        gms = gms[::10]
        models = models[::10]

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

        target_to_results = defaultdict(list)
        all_syn_trs = []

        dtraces = []
        for imodel in range(nmodels):
            model = models[imodel, :]

            source = problem.get_source(model)
            results = problem.evaluate(model)

            dtraces.append([])

            for target, result in zip(problem.targets, results):
                w = target.get_combined_weight()

                if isinstance(result, gf.SeismosizerError) or \
                        not isinstance(target, WaveformMisfitTarget) or \
                        not num.all(num.isfinite(w)):

                    dtraces[-1].extend([None] * target.nmisfits)
                    continue

                itarget, itarget_end = target_index[target]
                assert itarget_end == itarget + 1

                if target.misfit_config.domain == 'cc_max_norm':
                    tref = (
                        result.filtered_obs.tmin + result.filtered_obs.tmax) \
                        * 0.5

                    for tr_filt, tr_proc, tshift in (
                            (result.filtered_obs,
                             result.processed_obs,
                             0.),
                            (result.filtered_syn,
                             result.processed_syn,
                             result.tshift)):

                        norm = num.sum(num.abs(tr_proc.ydata)) \
                            / tr_proc.data_len()
                        tr_filt.ydata /= norm
                        tr_proc.ydata /= norm

                        tr_filt.shift(tshift)
                        tr_proc.shift(tshift)

                    ctr = result.cc
                    ctr.shift(tref)

                    dtrace = ctr

                else:
                    for tr in (
                            result.filtered_obs,
                            result.filtered_syn,
                            result.processed_obs,
                            result.processed_syn):

                        tr.ydata *= w

                    if result.tshift is not None and result.tshift != 0.0:
                        # result.filtered_syn.shift(result.tshift)
                        result.processed_syn.shift(result.tshift)

                    dtrace = make_norm_trace(
                        result.processed_syn, result.processed_obs,
                        problem.norm_exponent)

                target_to_results[target].append(result)

                dtrace.meta = dict(
                    normalisation_family=target.normalisation_family,
                    path=target.path)

                dtraces[-1].append(dtrace)

                result.processed_syn.meta = dict(
                    normalisation_family=target.normalisation_family,
                    path=target.path)

                all_syn_trs.append(result.processed_syn)

        if not all_syn_trs:
            logger.warn('No traces to show!')
            return

        def skey(tr):
            return tr.meta['normalisation_family'], tr.meta['path']

        trace_minmaxs = trace.minmax(all_syn_trs, skey)

        dtraces_all = []
        for dtraces_group in dtraces:
            dtraces_all.extend(dtraces_group)

        dminmaxs = trace.minmax([
            dtrace_ for dtrace_ in dtraces_all if dtrace_ is not None], skey)

        for tr in dtraces_all:
            if tr:
                dmin, dmax = dminmaxs[skey(tr)]
                tr.ydata /= max(abs(dmin), abs(dmax))

        cg_to_targets = meta.gather(
            problem.waveform_targets,
            lambda t: (t.path, t.codes[3]),
            filter=lambda t: t in target_to_results)

        cgs = sorted(cg_to_targets.keys())

        from matplotlib import colors
        cmap = cm.ScalarMappable(
            norm=colors.Normalize(vmin=num.min(icolor), vmax=num.max(icolor)),
            cmap=plt.get_cmap('coolwarm'))

        imodel_to_color = []
        for imodel in range(nmodels):
            imodel_to_color.append(cmap.to_rgba(icolor[imodel]))

        for cg in cgs:
            targets = cg_to_targets[cg]

            frame_to_target, nx, ny, nxx, nyy = layout(
                source, targets, nxmax, nymax)

            figures = {}
            for iy in range(ny):
                for ix in range(nx):
                    if (iy, ix) not in frame_to_target:
                        continue

                    ixx = ix // nxmax
                    iyy = iy // nymax
                    if (iyy, ixx) not in figures:
                        title = '_'.join(x for x in cg if x)
                        item = PlotItem(
                            name='fig_%s_%i_%i' % (title, ixx, iyy))
                        item.attributes['targets'] = []
                        figures[iyy, ixx] = (
                            item, plt.figure(figsize=self.size_inch))

                        figures[iyy, ixx][1].subplots_adjust(
                            left=0.03,
                            right=1.0 - 0.03,
                            bottom=0.03,
                            top=1.0 - 0.06,
                            wspace=0.2,
                            hspace=0.2)

                    item, fig = figures[iyy, ixx]

                    target = frame_to_target[iy, ix]

                    item.attributes['targets'].append(target.string_id())

                    amin, amax = trace_minmaxs[
                        target.normalisation_family, target.path]
                    absmax = max(abs(amin), abs(amax))

                    ny_this = nymax  # min(ny, nymax)
                    nx_this = nxmax  # min(nx, nxmax)
                    i_this = (iy % ny_this) * nx_this + (ix % nx_this) + 1

                    axes2 = fig.add_subplot(ny_this, nx_this, i_this)

                    space = 0.5
                    space_factor = 1.0 + space
                    axes2.set_axis_off()
                    axes2.set_ylim(-1.05 * space_factor, 1.05)

                    axes = axes2.twinx()
                    axes.set_axis_off()

                    if target.misfit_config.domain == 'cc_max_norm':
                        axes.set_ylim(-10. * space_factor, 10.)
                    else:
                        axes.set_ylim(-absmax*1.33 * space_factor, absmax*1.33)

                    itarget, itarget_end = target_index[target]
                    assert itarget_end == itarget + 1

                    for imodel, result in enumerate(target_to_results[target]):

                        syn_color = imodel_to_color[imodel]

                        dtrace = dtraces[imodel][itarget]

                        tap_color_annot = (0.35, 0.35, 0.25)
                        tap_color_edge = (0.85, 0.85, 0.80)
                        tap_color_fill = (0.95, 0.95, 0.90)

                        plot_taper(
                            axes2,
                            result.processed_obs.get_xdata(),
                            result.taper,
                            fc=tap_color_fill, ec=tap_color_edge, alpha=0.2)

                        obs_color = mpl_color('aluminium5')
                        obs_color_light = light(obs_color, 0.5)

                        plot_dtrace(
                            axes2, dtrace, space, 0., 1.,
                            fc='none',
                            ec=syn_color)

                        # plot_trace(
                        #     axes, result.filtered_syn,
                        #     color=syn_color_light, lw=1.0)

                        if imodel == 0:
                            plot_trace(
                                axes, result.filtered_obs,
                                color=obs_color_light, lw=0.75)

                        plot_trace(
                            axes, result.processed_syn,
                            color=syn_color, lw=1.0, alpha=0.3)

                        plot_trace(
                            axes, result.processed_obs,
                            color=obs_color, lw=0.75, alpha=0.3)

                        if imodel != 0:
                            continue
                        xdata = result.filtered_obs.get_xdata()
                        axes.set_xlim(xdata[0], xdata[-1])

                        tmarks = [
                            result.processed_obs.tmin,
                            result.processed_obs.tmax]

                        for tmark in tmarks:
                            axes2.plot(
                                [tmark, tmark], [-0.9, 0.1],
                                color=tap_color_annot)

                        dur = tmarks[1] - tmarks[0]
                        for tmark, text, ha in [
                                (tmarks[0],
                                 '$\\,$ ' + meta.str_duration(
                                    tmarks[0] - source.time),
                                 'left'),
                                (tmarks[1],
                                 '$\\Delta$ ' + meta.str_duration(
                                    dur),
                                 'right')]:

                            axes2.annotate(
                                text,
                                xy=(tmark, -0.9),
                                xycoords='data',
                                xytext=(
                                    fontsize*0.4 * [-1, 1][ha == 'left'],
                                    fontsize*0.2),
                                textcoords='offset points',
                                ha=ha,
                                va='bottom',
                                color=tap_color_annot,
                                fontsize=fontsize)

                        axes2.set_xlim(
                            tmarks[0] - dur*0.1, tmarks[1] + dur*0.1)

                    scale_string = None

                    if target.misfit_config.domain == 'cc_max_norm':
                        scale_string = 'Syn/obs scales differ!'

                    infos = []
                    if scale_string:
                        infos.append(scale_string)

                    if self.nx == 1 and self.ny == 1:
                        infos.append(target.string_id())
                    else:
                        infos.append('.'.join(x for x in target.codes if x))
                    dist = source.distance_to(target)
                    azi = source.azibazi_to(target)[0]
                    infos.append(meta.str_dist(dist))
                    infos.append(u'%.0f\u00B0' % azi)
                    axes2.annotate(
                        '\n'.join(infos),
                        xy=(0., 1.),
                        xycoords='axes fraction',
                        xytext=(2., 2.),
                        textcoords='offset points',
                        ha='left',
                        va='top',
                        fontsize=fontsize,
                        fontstyle='normal')

                    if (self.nx == 1 and self.ny == 1):
                        yield item, fig
                        del figures[iyy, ixx]

            if not (self.nx == 1 and self.ny == 1):
                for (iyy, ixx), (_, fig) in figures.items():
                    title = '.'.join(x for x in cg if x)
                    if len(figures) > 1:
                        title += ' (%i/%i, %i/%i)' % (iyy+1, nyy, ixx+1, nxx)

                    fig.suptitle(title, fontsize=fontsize_title)

            for item, fig in figures.values():
                yield item, fig


class FitsWaveformPlot(PlotConfig):
    ''' Plot showing the waveform fits for the best model '''
    name = 'fits_waveform'
    size_cm = Tuple.T(
        2, Float.T(),
        default=(9., 5.),
        help='width and length of the figure in cm')
    nx = Int.T(
        default=1,
        help='horizontal number of subplots on every page')
    ny = Int.T(
        default=1,
        help='vertical number of subplots on every page')
    font_size = Float.T(
        default=8,
        help='Font Size of all fonts, except title')
    font_size_title = Float.T(
        default=10,
        help='Font Size of title')

    def make(self, environ):
        cm = environ.get_plot_collection_manager()
        mpl_init(fontsize=self.font_size)
        environ.setup_modelling()
        ds = environ.get_dataset()
        optimiser = environ.get_optimiser()

        environ.setup_modelling()

        history = environ.get_history(subset='harvest')
        cm.create_group_mpl(
            self,
            self.draw_figures(ds, history, optimiser),
            title=u'Waveform fits for best model',
            section='fits',
            feather_icon='activity',
            description=u'''
Plot showing observed and synthetic waveform (attributes) for the best fitting
model.

Best model's waveform fits for all targets. Depending on the target
configurations different types of comparisons are possible: (i) time domain
waveform differences, (ii) amplitude spectra, (iii) envelopes, (iv) cross
correlation functions. Each waveform plot gives a number of details:

1) Target information (left side, from top to bottom) gives station name with
component, distance to source, azimuth of station with respect to source,
target weight, target misfit and starting time of the waveform relative to the
origin time.

2) The background gray area shows the applied taper function.

3) The waveforms shown are: the restituted and filtered observed trace without
tapering (light grey) and the same trace with tapering and processing (dark
gray), the synthetic trace (light red) and the filtered, tapered and (if
enabled) shifted and processed synthetic target trace (red). The traces are
scaled according to the target weight (small weight, small amplitude) and
normed relative to the maximum amplitude of the targets of the corresponding
normalisation family.

4) The bottom panel shows, depending on the type of comparison, sample-wise
residuals for time domain comparisons (red filled), spectra of observed and
synthetic traces for amplitude spectrum comparisons, or cross correlation
traces.

5) Colored boxes on the upper right show the relative weight of the target
within the entire dataset of the optimisation (top box, orange) and the
relative misfit contribution to the global misfit of the optimisation (bottom
box, red).
''')

    def draw_figures(self, ds, history, optimiser):

        fontsize = self.font_size
        fontsize_title = self.font_size_title

        nxmax = self.nx
        nymax = self.ny

        problem = history.problem

        for target in problem.targets:
            target.set_dataset(ds)

        target_index = {}
        i = 0
        for target in problem.targets:
            target_index[target] = i, i+target.nmisfits
            i += target.nmisfits

        xbest = history.get_best_model()
        misfits = history.misfits[history.get_sorted_misfits_idx(chain=0), ...]

        ws = problem.get_target_weights()

        gcms = problem.combine_misfits(
            misfits[:1, :, :],
            extra_correlated_weights=optimiser.get_correlated_weights(problem),
            get_contributions=True)[0, :]

        w_max = num.nanmax(ws)
        gcm_max = num.nanmax(gcms)

        source = problem.get_source(xbest)

        target_to_result = {}
        all_syn_trs = []
        all_syn_specs = []
        results = problem.evaluate(xbest)

        dtraces = []
        for target, result in zip(problem.targets, results):
            if not isinstance(result, WaveformMisfitResult):
                dtraces.extend([None] * target.nmisfits)
                continue

            itarget, itarget_end = target_index[target]
            assert itarget_end == itarget + 1

            w = target.get_combined_weight()

            if target.misfit_config.domain == 'cc_max_norm':
                tref = (
                    result.filtered_obs.tmin + result.filtered_obs.tmax) * 0.5
                for tr_filt, tr_proc, tshift in (
                        (result.filtered_obs,
                         result.processed_obs,
                         0.),
                        (result.filtered_syn,
                         result.processed_syn,
                         result.tshift)):

                    norm = num.sum(num.abs(tr_proc.ydata)) / tr_proc.data_len()
                    tr_filt.ydata /= norm
                    tr_proc.ydata /= norm

                    tr_filt.shift(tshift)
                    tr_proc.shift(tshift)

                ctr = result.cc
                ctr.shift(tref)

                dtrace = ctr

            else:
                for tr in (
                        result.filtered_obs,
                        result.filtered_syn,
                        result.processed_obs,
                        result.processed_syn):

                    tr.ydata *= w

                for spec in (
                        result.spectrum_obs,
                        result.spectrum_syn):

                    if spec is not None:
                        spec.ydata *= w

                if result.tshift is not None and result.tshift != 0.0:
                    # result.filtered_syn.shift(result.tshift)
                    result.processed_syn.shift(result.tshift)

                dtrace = make_norm_trace(
                    result.processed_syn, result.processed_obs,
                    problem.norm_exponent)

            target_to_result[target] = result

            dtrace.meta = dict(
                normalisation_family=target.normalisation_family,
                path=target.path)
            dtraces.append(dtrace)

            result.processed_syn.meta = dict(
                normalisation_family=target.normalisation_family,
                path=target.path)

            all_syn_trs.append(result.processed_syn)

            if result.spectrum_syn:
                result.spectrum_syn.meta = dict(
                    normalisation_family=target.normalisation_family,
                    path=target.path)

                all_syn_specs.append(result.spectrum_syn)

        if not all_syn_trs:
            logger.warn('No traces to show!')
            return

        def skey(tr):
            return tr.meta['normalisation_family'], tr.meta['path']

        trace_minmaxs = trace.minmax(all_syn_trs, skey)

        amp_spec_maxs = amp_spec_max(all_syn_specs, skey)

        dminmaxs = trace.minmax([x for x in dtraces if x is not None], skey)

        for tr in dtraces:
            if tr:
                dmin, dmax = dminmaxs[skey(tr)]
                tr.ydata /= max(abs(dmin), abs(dmax))

        cg_to_targets = meta.gather(
            problem.waveform_targets,
            lambda t: (t.path, t.codes[3]),
            filter=lambda t: t in target_to_result)

        cgs = sorted(cg_to_targets.keys())

        for cg in cgs:
            targets = cg_to_targets[cg]

            frame_to_target, nx, ny, nxx, nyy = layout(
                source, targets, nxmax, nymax)

            figures = {}
            for iy in range(ny):
                for ix in range(nx):
                    if (iy, ix) not in frame_to_target:
                        continue

                    ixx = ix // nxmax
                    iyy = iy // nymax
                    if (iyy, ixx) not in figures:
                        title = '_'.join(x for x in cg if x)
                        item = PlotItem(
                            name='fig_%s_%i_%i' % (title, ixx, iyy))
                        item.attributes['targets'] = []
                        figures[iyy, ixx] = (
                            item, plt.figure(figsize=self.size_inch))

                        figures[iyy, ixx][1].subplots_adjust(
                            left=0.03,
                            right=1.0 - 0.03,
                            bottom=0.03,
                            top=1.0 - 0.06,
                            wspace=0.2,
                            hspace=0.2)

                    item, fig = figures[iyy, ixx]

                    target = frame_to_target[iy, ix]

                    item.attributes['targets'].append(target.string_id())

                    amin, amax = trace_minmaxs[
                        target.normalisation_family, target.path]
                    absmax = max(abs(amin), abs(amax))

                    ny_this = nymax  # min(ny, nymax)
                    nx_this = nxmax  # min(nx, nxmax)
                    i_this = (iy % ny_this) * nx_this + (ix % nx_this) + 1

                    axes2 = fig.add_subplot(ny_this, nx_this, i_this)

                    space = 0.5
                    space_factor = 1.0 + space
                    axes2.set_axis_off()
                    axes2.set_ylim(-1.05 * space_factor, 1.05)

                    axes = axes2.twinx()
                    axes.set_axis_off()

                    if target.misfit_config.domain == 'cc_max_norm':
                        axes.set_ylim(-10. * space_factor, 10.)
                    else:
                        axes.set_ylim(
                            -absmax * 1.33 * space_factor, absmax * 1.33)

                    itarget, itarget_end = target_index[target]
                    assert itarget_end == itarget + 1

                    result = target_to_result[target]

                    dtrace = dtraces[itarget]

                    tap_color_annot = (0.35, 0.35, 0.25)
                    tap_color_edge = (0.85, 0.85, 0.80)
                    tap_color_fill = (0.95, 0.95, 0.90)

                    plot_taper(
                        axes2, result.processed_obs.get_xdata(), result.taper,
                        fc=tap_color_fill, ec=tap_color_edge)

                    obs_color = mpl_color('aluminium5')
                    obs_color_light = light(obs_color, 0.5)

                    syn_color = mpl_color('scarletred2')
                    syn_color_light = light(syn_color, 0.5)

                    misfit_color = mpl_color('scarletred2')
                    weight_color = mpl_color('chocolate2')

                    cc_color = mpl_color('aluminium5')

                    if target.misfit_config.domain == 'cc_max_norm':
                        tref = (result.filtered_obs.tmin +
                                result.filtered_obs.tmax) * 0.5

                        plot_dtrace(
                            axes2, dtrace, space, -1., 1.,
                            fc=light(cc_color, 0.5),
                            ec=cc_color)

                        plot_dtrace_vline(
                            axes2, tref, space, color=tap_color_annot)

                    elif target.misfit_config.domain == 'frequency_domain':

                        asmax = amp_spec_maxs[
                            target.normalisation_family, target.path]
                        fmin, fmax = \
                            target.misfit_config.get_full_frequency_range()

                        plot_spectrum(
                            axes2,
                            result.spectrum_syn,
                            result.spectrum_obs,
                            fmin, fmax,
                            space, 0., asmax,
                            syn_color=syn_color,
                            obs_color=obs_color,
                            syn_lw=1.0,
                            obs_lw=0.75,
                            color_vline=tap_color_annot,
                            fontsize=fontsize)

                    else:
                        plot_dtrace(
                            axes2, dtrace, space, 0., 1.,
                            fc=light(misfit_color, 0.3),
                            ec=misfit_color)

                    plot_trace(
                        axes, result.filtered_syn,
                        color=syn_color_light, lw=1.0)

                    plot_trace(
                        axes, result.filtered_obs,
                        color=obs_color_light, lw=0.75)

                    plot_trace(
                        axes, result.processed_syn,
                        color=syn_color, lw=1.0)

                    plot_trace(
                        axes, result.processed_obs,
                        color=obs_color, lw=0.75)

                    # xdata = result.filtered_obs.get_xdata()

                    tmarks = [
                        result.processed_obs.tmin,
                        result.processed_obs.tmax]

                    for tmark in tmarks:
                        axes2.plot(
                            [tmark, tmark], [-0.9, 0.1], color=tap_color_annot)

                    dur = tmarks[1] - tmarks[0]
                    for tmark, text, ha in [
                            (tmarks[0],
                             '$\\,$ ' + meta.str_duration(
                                 tmarks[0] - source.time),
                             'left'),
                            (tmarks[1],
                             '$\\Delta$ ' + meta.str_duration(dur),
                             'right')]:

                        axes2.annotate(
                            text,
                            xy=(tmark, -0.9),
                            xycoords='data',
                            xytext=(
                                fontsize * 0.4 * [-1, 1][ha == 'left'],
                                fontsize * 0.2),
                            textcoords='offset points',
                            ha=ha,
                            va='bottom',
                            color=tap_color_annot,
                            fontsize=fontsize)

                    axes2.set_xlim(tmarks[0] - dur*0.1, tmarks[1] + dur*0.1)

                    rel_w = ws[itarget] / w_max
                    rel_c = gcms[itarget] / gcm_max

                    sw = 0.25
                    sh = 0.1
                    ph = 0.01

                    for (ih, rw, facecolor, edgecolor) in [
                            (0, rel_w, light(weight_color, 0.5),
                             weight_color),
                            (1, rel_c, light(misfit_color, 0.5),
                             misfit_color)]:

                        bar = patches.Rectangle(
                            (1.0 - rw * sw, 1.0 - (ih + 1) * sh + ph),
                            rw * sw,
                            sh - 2 * ph,
                            facecolor=facecolor, edgecolor=edgecolor,
                            zorder=10,
                            transform=axes.transAxes, clip_on=False)

                        axes.add_patch(bar)

                    scale_string = None

                    if target.misfit_config.domain == 'cc_max_norm':
                        scale_string = 'Syn/obs scales differ!'

                    infos = []
                    if scale_string:
                        infos.append(scale_string)

                    if self.nx == 1 and self.ny == 1:
                        infos.append(target.string_id())
                    else:
                        infos.append('.'.join(x for x in target.codes if x))

                    dist = source.distance_to(target)
                    azi = source.azibazi_to(target)[0]
                    infos.append(meta.str_dist(dist))
                    infos.append('%.0f\u00B0' % azi)
                    infos.append('%.3g' % ws[itarget])
                    infos.append('%.3g' % gcms[itarget])
                    axes2.annotate(
                        '\n'.join(infos),
                        xy=(0., 1.),
                        xycoords='axes fraction',
                        xytext=(2., 2.),
                        textcoords='offset points',
                        ha='left',
                        va='top',
                        fontsize=fontsize,
                        fontstyle='normal')

                    if (self.nx == 1 and self.ny == 1):
                        yield item, fig
                        del figures[iyy, ixx]

            if not (self.nx == 1 and self.ny == 1):
                for (iyy, ixx), (_, fig) in figures.items():
                    title = '.'.join(x for x in cg if x)
                    if len(figures) > 1:
                        title += ' (%i/%i, %i/%i)' % (
                            iyy + 1, nyy, ixx + 1, nxx)

                    fig.suptitle(title, fontsize=fontsize_title)

            for item, fig in figures.values():
                yield item, fig


class WaveformStationDistribution(StationDistributionPlot):
    ''' Plot showing all waveform fits for the ensemble of solutions'''

    name = 'seismic_stations'

    def make(self, environ):
        from grond.problems.base import ProblemDataNotAvailable
        from grond.environment import NoRundirAvailable

        cm = environ.get_plot_collection_manager()
        mpl_init(fontsize=self.font_size)

        problem = environ.get_problem()
        dataset = environ.get_dataset()
        try:
            history = environ.get_history(subset='harvest')
        except (NoRundirAvailable, ProblemDataNotAvailable):
            history = None

        cm.create_group_mpl(
            self,
            self.draw_figures(problem, dataset, history),
            title=u'Seismic station locations',
            section='checks',
            feather_icon='target',
            description=u'''
Plot showing seismic station locations and attributes.

Station locations in dependence of distance and azimuth are shown. The center
of the plot corresponds to the origin of the search space, not to the optimised
location of the source.
''')

    def draw_figures(self, problem, dataset, history):

        target_index = {}
        i = 0
        for target in problem.targets:
            target_index[target] = i, i+target.nmisfits
            i += target.nmisfits

        ws = problem.get_target_weights()

        if history:
            misfits = history.misfits[history.get_sorted_misfits_idx(), ...]
            gcms = problem.combine_misfits(
                misfits[:1, :, :], get_contributions=True)[0, :]

        event = problem.base_source

        cg_to_targets = meta.gather(
            problem.waveform_targets,
            lambda t: (t.path, t.codes[3]))

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

            azimuths = num.array([event.azibazi_to(t)[0] for t in targets])
            distances = num.array([t.distance_to(event) for t in targets])

            item = PlotItem(
                name='seismic_stations_weights_%s' % cg_str,
                title=u'Station weights (%s)' % cg_str,
                description=u'\n\nMarkers are scaled according to the '
                            u'weighting factor of the corresponding target\'s '
                            u'contribution in the misfit function.')
            fig, ax, legend = self.plot_station_distribution(
                azimuths, distances, ws[itargets], labels)
            legend.set_title(
                'Weight',
                prop=dict(size=self.font_size))

            yield (item, fig)

            if history:
                item = PlotItem(
                    name='seismic_stations_contributions_%s' % cg_str,
                    title=u'Station misfit contributions (%s)' % cg_str,
                    description=u'\n\nMarkers are scaled according to their '
                                u'misfit contribution for the globally best '
                                u'source model.')
                fig, ax, legend = self.plot_station_distribution(
                    azimuths, distances, gcms[itargets], labels)
                legend.set_title(
                    'Contribution',
                    prop=dict(size=self.font_size))

                yield (item, fig)


def get_plot_classes():
    return [
        WaveformStationDistribution,
        CheckWaveformsPlot,
        FitsWaveformPlot,
        FitsWaveformEnsemblePlot
    ]
