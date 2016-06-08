import math
import random
import logging
import os.path as op
import numpy as num
from scipy import signal
from pyrocko import automap, beachball, guts, trace, util
from pyrocko import hudson
from grond import core
from matplotlib import pyplot as plt
from matplotlib import cm, patches
from pyrocko.cake_plot import mpl_init, labelspace, colors, \
    str_to_mpl_color as scolor, light

logger = logging.getLogger('grond.plot')

km = 1000.


def ordersort(x):
    isort = num.argsort(x)
    iorder = num.empty(isort.size)
    iorder[isort] = num.arange(isort.size)
    return iorder


def nextpow2(i):
    return 2**int(math.ceil(math.log(i)/math.log(2.)))


def fixlim(lo, hi):
    if lo == hi:
        return lo - 1.0, hi + 1.0
    else:
        return lo, hi


def str_dist(dist):
    if dist < 10.0:
        return '%g m' % dist
    elif 10. <= dist < 1.*km:
        return '%.0f m' % dist
    elif 1.*km <= dist < 10.*km:
        return '%.1f km' % (dist / km)
    else:
        return '%.0f km' % (dist / km)


def str_duration(t):
    if t < 1.0:
        return '%g s' % t
    elif 1.0 <= t < 10.0:
        return '%.1f s' % t
    elif 10.0 <= t < 3600.:
        return util.time_to_str(t, format='%M:%S min')
    elif 3600. <= t < 24*3600.:
        return util.time_to_str(t, format='%H:%M h')
    else:
        return '%.1f d' % (t / (24.*3600.))


def str_duration2(t):
    if t < 1.0:
        return '%g s' % t
    elif 1.0 <= t < 10.0:
        return '%.1f s' % t
    elif 10.0 <= t < 3600.:
        return '%.0f s' % t
    elif 10.0 <= t:
        return '%.1f h' % (t / 3600.)


def eigh_sorted(mat):
    evals, evecs = num.linalg.eigh(mat)
    iorder = num.argsort(evals)
    return evals[iorder], evecs[:, iorder]


def plot(stations, center_lat, center_lon, radius, output_path,
         width=25., height=25.,
         show_station_labels=False):

    station_lats = num.array([s.lat for s in stations])
    station_lons = num.array([s.lon for s in stations])

    map = automap.Map(
        width=width,
        height=height,
        lat=center_lat,
        lon=center_lon,
        radius=radius,
        show_rivers=False,
        show_topo=False,
        illuminate_factor_land=0.35,
        color_dry=(240, 240, 235),
        topo_cpt_wet='white_sea_land',
        topo_cpt_dry='white_sea_land')

    map.gmt.psxy(
        in_columns=(station_lons, station_lats),
        S='t8p',
        G='black',
        *map.jxyr)

    if show_station_labels:
        for s in stations:
            map.add_label(s.lat, s.lon, '%s' % s.station)

    map.save(output_path)


def map_geometry(config, output_path):
    stations = config.get_dataset().get_stations()

    lat0, lon0, radius = core.stations_mean_latlondist(stations)

    radius *= 1.5

    plot(stations, lat0, lon0, radius, output_path,
         show_station_labels=True)


class GrondModel(object):
    def __init__(self, **kwargs):
        self.listeners = []
        self.set_problem(None)

    def add_listener(self, listener):
        self.listeners.append(listener)

    def set_problem(self, problem):

        self.problem = problem
        if problem:
            nparameters = problem.nparameters
            ntargets = problem.ntargets
        else:
            nparameters = 0
            ntargets = 0

        nmodels = 0
        nmodels_capacity = 1024

        self._xs_buffer = num.zeros(
            (nmodels_capacity, nparameters), dtype=num.float)
        self._misfits_buffer = num.zeros(
            (nmodels_capacity, ntargets, 2), dtype=num.float)

        self.xs = self._xs_buffer[:nmodels, :]
        self.misfits = self._misfits_buffer[:nmodels, :, :]

        self.data_changed()

    @property
    def nmodels(self):
        return self.xs.shape[0]

    @property
    def nmodels_capacity(self):
        return self._xs_buffer.shape[0]

    def append(self, xs, misfits):
        assert xs.shape[0] == misfits.shape[0]

        nmodels_add = xs.shape[0]

        nmodels = self.nmodels
        nmodels_new = nmodels + nmodels_add
        nmodels_capacity_new = max(1024, nextpow2(nmodels_new))

        nmodels_capacity = self.nmodels_capacity
        if nmodels_capacity_new > nmodels_capacity:
            xs_buffer = num.zeros(
                (nmodels_capacity_new, self.problem.nparameters),
                dtype=num.float)

            misfits_buffer = num.zeros(
                (nmodels_capacity_new, self.problem.ntargets, 2),
                dtype=num.float)

            xs_buffer[:nmodels, :] = self._xs_buffer[:nmodels]
            misfits_buffer[:nmodels, :] = self._misfits_buffer[:nmodels]
            self._xs_buffer = xs_buffer
            self._misfits_buffer = misfits_buffer

        self._xs_buffer[nmodels:nmodels+nmodels_add, :] = xs
        self._misfits_buffer[nmodels:nmodels+nmodels_add, :, :] = misfits

        nmodels = nmodels_new

        self.xs = self._xs_buffer[:nmodels, :]
        self.misfits = self._misfits_buffer[:nmodels, :, :]

        self.data_changed()

    def data_changed(self):
        for listener in self.listeners:
            listener()


def draw_sequence_figures(model, plt, misfit_cutoff=None):
    problem = model.problem
    if not problem:
        return

    imodels = num.arange(model.nmodels)
    bounds = problem.bounds() + problem.dependant_bounds()

    xref = problem.pack(problem.base_source)

    xs = model.xs

    npar = problem.nparameters
    ndep = problem.ndependants

    gms = problem.global_misfits(model.misfits)
    gms_softclip = num.where(gms > 1.0, 0.2 * num.log10(gms) + 1.0, gms)

    isort = num.argsort(gms)[::-1]

    imodels = imodels[isort]
    gms = gms[isort]
    gms_softclip = gms_softclip[isort]
    xs = xs[isort, :]

    iorder = num.empty_like(isort)
    iorder = num.arange(iorder.size)

    if misfit_cutoff is None:
        ibest = num.ones(gms.size, dtype=num.bool)
    else:
        ibest = gms < misfit_cutoff

    nfx = 2
    nfy = 4
    # nfz = (npar + ndep + 1 - 1) / (nfx*nfy) + 1
    cmap = cm.YlOrRd
    cmap = cm.jet
    axes = None
    fig = None
    alpha = 0.5
    for ipar in xrange(npar):
        impl = ipar % (nfx*nfy) + 1

        if impl == 1:
            fig = plt.figure()

        par = problem.parameters[ipar]

        axes = fig.add_subplot(nfy, nfx, impl, sharex=axes)
        axes.set_ylabel(par.get_label())
        axes.get_yaxis().set_major_locator(plt.MaxNLocator(4))
        if impl < (nfx*nfy-1):
            axes.get_xaxis().set_visible(False)
        else:
            axes.set_xlabel('Iteration')

        axes.set_ylim(*fixlim(*par.scaled(bounds[ipar])))
        axes.set_xlim(0, model.nmodels)
        axes.axhline(par.scaled(xref[ipar]), color='black', alpha=0.3)

        axes.scatter(
            imodels[ibest], par.scaled(xs[ibest, ipar]), s=3, c=iorder[ibest],
            lw=0, cmap=cmap, alpha=alpha)

    for idep in xrange(ndep):
        # ifz, ify, ifx = num.unravel_index(ipar, (nfz, nfy, nfx))
        impl = (npar+idep) % (nfx*nfy) + 1

        if impl == 1:
            fig = plt.figure()

        par = problem.dependants[idep]

        axes = fig.add_subplot(nfy, nfx, impl, sharex=axes)
        axes.set_ylabel(par.get_label())
        axes.get_yaxis().set_major_locator(plt.MaxNLocator(4))
        if impl < (nfx*nfy-1):
            axes.get_xaxis().set_visible(False)
        else:
            axes.set_xlabel('Iteration')
        axes.set_ylim(*fixlim(*par.scaled(bounds[npar+idep])))
        axes.set_xlim(0, model.nmodels)

        y = problem.make_dependant(xref, par.name)
        axes.axhline(par.scaled(y), color='black', alpha=0.3)

        ys = problem.make_dependant(xs[ibest, :], par.name)
        axes.scatter(
            imodels[ibest], par.scaled(ys), s=3, c=iorder[ibest],
            lw=0, cmap=cmap, alpha=alpha)

    impl = (npar+ndep) % (nfx*nfy) + 1
    if impl == 1:
        fig = plt.figure()

    axes = fig.add_subplot(nfy, nfx, impl, sharex=axes)

    axes.set_ylim(0., 1.5)
    axes.axhspan(1.0, 1.5, color=(0.8, 0.8, 0.8), alpha=0.2)
    axes.axhline(1.0, color=(0.5, 0.5, 0.5), zorder=2)
    axes.set_yticks([0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
    axes.set_yticklabels(['0.0', '0.2', '0.4', '0.6', '0.8', '1', '10', '100'])

    axes.scatter(
        imodels[ibest], gms_softclip[ibest], c=iorder[ibest],
        s=3, lw=0, cmap=cmap, alpha=alpha)

    axes.set_xlim(0, model.nmodels)
    axes.set_xlabel('Iteration')

    axes.set_ylabel('Misfit')

    fig.canvas.draw()


def draw_jointpar_figures(
        model, plt, misfit_cutoff=None, ibootstrap=None, color=None,
        exclude=None, include=None):

    color = 'misfit'
    # exclude = ['duration']
    # include = ['magnitude', 'rel_moment_iso', 'rel_moment_clvd', 'depth']
    neach = 6
    figsize = (8, 8)
    # cmap = cm.YlOrRd
    # cmap = cm.jet
    cmap = cm.coolwarm

    problem = model.problem
    if not problem:
        return

    xs = model.xs

    bounds = problem.bounds() + problem.dependant_bounds()
    for ipar in xrange(problem.ncombined):
        par = problem.combined[ipar]
        lo, hi = bounds[ipar]
        if lo == hi:
            if exclude is None:
                exclude = []

            exclude.append(par.name)

    xref = problem.pack(problem.base_source)

    if ibootstrap is not None:
        gms = problem.bootstrap_misfits(model.misfits, ibootstrap)
    else:
        gms = problem.global_misfits(model.misfits)

    isort = num.argsort(gms)[::-1]

    gms = gms[isort]
    xs = xs[isort, :]

    if misfit_cutoff is not None:
        ibest = gms < misfit_cutoff
        gms = gms[ibest]
        xs = xs[ibest]

    nmodels = xs.shape[0]

    if color == 'dist':
        mx = num.mean(xs, axis=0)
        cov = num.cov(xs.T)
        mdists = core.mahalanobis_distance(xs, mx, cov)
        color = ordersort(mdists)

    elif color == 'misfit':
        iorder = num.arange(nmodels)
        color = iorder

    elif color in problem.parameter_names:
        ind = problem.name_to_index(color)
        color = ordersort(problem.extract(xs, ind))

    smap = {}
    iselected = 0
    for ipar in xrange(problem.ncombined):
        par = problem.combined[ipar]
        if exclude and par.name in exclude or \
                include and par.name not in include:
            continue

        smap[iselected] = ipar
        iselected += 1

    nselected = iselected

    if nselected == 0:
        return

    nfig = (nselected-2) / neach + 1

    figs = []
    for ifig in xrange(nfig):
        figs_row = []
        for jfig in xrange(nfig):
            if ifig >= jfig:
                figs_row.append(plt.figure(figsize=figsize))
            else:
                figs_row.append(None)

        figs.append(figs_row)

    for iselected in xrange(nselected):
        ipar = smap[iselected]
        ypar = problem.combined[ipar]
        for jselected in xrange(iselected):
            jpar = smap[jselected]
            xpar = problem.combined[jpar]

            ixg = (iselected - 1)
            iyg = jselected

            ix = ixg % neach
            iy = iyg % neach

            ifig = ixg/neach
            jfig = iyg/neach

            aind = (neach, neach, (ix * neach) + iy + 1)

            fig = figs[ifig][jfig]

            axes = fig.add_subplot(*aind)

            axes.axvline(0., color=scolor('aluminium3'), lw=0.5)
            axes.axhline(0., color=scolor('aluminium3'), lw=0.5)
            for spine in axes.spines.values():
                spine.set_edgecolor(scolor('aluminium5'))
                spine.set_linewidth(0.5)

            xmin, xmax = fixlim(*xpar.scaled(bounds[jpar]))
            ymin, ymax = fixlim(*ypar.scaled(bounds[ipar]))

            if ix == 0 or jselected + 1 == iselected:
                for (xpos, xoff, x) in [(0.0, 10., xmin), (1.0, -10., xmax)]:
                    axes.annotate(
                        '%.2g%s' % (x, xpar.get_unit_suffix()),
                        xy=(xpos, 1.05),
                        xycoords='axes fraction',
                        xytext=(xoff, 5.),
                        textcoords='offset points',
                        verticalalignment='bottom',
                        horizontalalignment='left',
                        rotation=45.)

            if iy == neach - 1 or jselected + 1 == iselected:
                for (ypos, yoff, y) in [(0., 10., ymin), (1.0, -10., ymax)]:
                    axes.annotate(
                        '%.2g%s' % (y, ypar.get_unit_suffix()),
                        xy=(1.0, ypos),
                        xycoords='axes fraction',
                        xytext=(5., yoff),
                        textcoords='offset points',
                        verticalalignment='bottom',
                        horizontalalignment='left',
                        rotation=45.)

            axes.set_xlim(xmin, xmax)
            axes.set_ylim(ymin, ymax)

            axes.get_xaxis().set_ticks([])
            axes.get_yaxis().set_ticks([])

            if iselected == nselected - 1 or ix == neach - 1:
                axes.annotate(
                    xpar.get_label(with_unit=False),
                    xy=(0.5, -0.05),
                    xycoords='axes fraction',
                    verticalalignment='top',
                    horizontalalignment='right',
                    rotation=45.)

            if iy == 0:
                axes.annotate(
                    ypar.get_label(with_unit=False),
                    xy=(-0.05, 0.5),
                    xycoords='axes fraction',
                    verticalalignment='top',
                    horizontalalignment='right',
                    rotation=45.)

            fx = problem.extract(xs, jpar)
            fy = problem.extract(xs, ipar)

            axes.scatter(
                xpar.scaled(fx),
                ypar.scaled(fy),
                c=color,
                s=3, alpha=0.5, cmap=cmap, edgecolors='none')

            cov = num.cov((xpar.scaled(fx), ypar.scaled(fy)))
            evals, evecs = eigh_sorted(cov)
            evals = num.sqrt(evals)
            ell = patches.Ellipse(
                xy=(num.mean(xpar.scaled(fx)), num.mean(ypar.scaled(fy))),
                width=evals[0]*2,
                height=evals[1]*2,
                angle=num.rad2deg(num.arctan2(evecs[1][0], evecs[0][0])))

            ell.set_facecolor('none')
            axes.add_artist(ell)

            fx = problem.extract(xref, jpar)
            fy = problem.extract(xref, ipar)

            ref_color = scolor('aluminium6')
            ref_color_light = 'none'
            axes.plot(
                xpar.scaled(fx), ypar.scaled(fy), 's',
                mew=1.5, ms=5, color=ref_color_light, mec=ref_color)


def draw_solution_figure(
        model, plt, misfit_cutoff=None, beachball_type='full'):

    fontsize = 10.

    fig = plt.figure(figsize=(6, 2))
    axes = fig.add_subplot(1, 1, 1, aspect=1.0)
    fig.subplots_adjust(left=0., right=1., bottom=0., top=1.)

    problem = model.problem
    if not problem:
        return

    xs = model.xs

    if xs.size == 0:
        return

    gms = problem.global_misfits(model.misfits)
    isort = num.argsort(gms)
    iorder = num.empty_like(isort)
    iorder[isort] = num.arange(iorder.size)[::-1]

    mean_source = core.get_mean_source(problem, model.xs)
    best_source = core.get_best_source(problem, model.xs, model.misfits)
    ref_source = problem.base_source

    for xpos, label in [
            (0., 'Full'),
            (2., 'Isotropic'),
            (4., 'Deviatoric'),
            (6., 'CLVD'),
            (8., 'DC')]:

        axes.annotate(
            label,
            xy=(1+xpos, 3),
            xycoords='data',
            xytext=(0., 0.),
            textcoords='offset points',
            ha='center',
            va='center',
            color='black',
            fontsize=fontsize)

    decos = []
    for source in [best_source, mean_source, ref_source]:
        mt = source.pyrocko_moment_tensor()
        deco = mt.standard_decomposition()
        decos.append(deco)

    moment_full_max = max(deco[-1][0] for deco in decos)

    for ypos, label, deco, color_t in [
            (2., 'Ensemble best', decos[0], scolor('aluminium5')),
            (1., 'Ensemble mean', decos[1], scolor('scarletred1')),
            (0., 'Reference', decos[2], scolor('aluminium3'))]:

        [(moment_iso, ratio_iso, m_iso),
         (moment_dc, ratio_dc, m_dc),
         (moment_clvd, ratio_clvd, m_clvd),
         (moment_devi, ratio_devi, m_devi),
         (moment_full, ratio_full, m_full)] = deco

        size0 = moment_full / moment_full_max

        axes.annotate(
            label,
            xy=(-2., ypos),
            xycoords='data',
            xytext=(0., 0.),
            textcoords='offset points',
            ha='left',
            va='center',
            color='black',
            fontsize=fontsize)

        for xpos, mt_part, ratio, ops in [
                (0., m_full, ratio_full, '-'),
                (2., m_iso, ratio_iso, '='),
                (4., m_devi, ratio_devi, '='),
                (6., m_clvd, ratio_clvd, '+'),
                (8., m_dc, ratio_dc, None)]:

            if ratio != 0.0:
                beachball.plot_beachball_mpl(
                    mt_part, axes,
                    beachball_type='full',
                    position=(1.+xpos, ypos),
                    size=0.9*size0*math.sqrt(ratio),
                    size_units='data',
                    color_t=color_t,
                    linewidth=1.0)

            else:
                axes.annotate(
                    'N/A',
                    xy=(1.+xpos, ypos),
                    ha='center',
                    va='center',
                    color='black',
                    fontsize=fontsize)

            if ops is not None:
                axes.annotate(
                    ops,
                    xy=(2. + xpos, ypos),
                    ha='center',
                    va='center',
                    color='black',
                    fontsize=fontsize)

    axes.axison = False
    axes.set_xlim(-2.25, 9.75)
    axes.set_ylim(-0.5, 3.5)
    fig.savefig('test.pdf')


def draw_contributions_figure(model, plt):

    fig = plt.figure()

    problem = model.problem
    if not problem:
        return

    xs = model.xs

    if xs.size == 0:
        return

    imodels = num.arange(model.nmodels)

    gms = problem.global_misfits(model.misfits)**2

    isort = num.argsort(gms)[::-1]

    gms = gms[isort]

    gms_softclip = num.where(gms > 1.0, 0.1 * num.log10(gms) + 1.0, gms)

    gcms = problem.global_contributions(model.misfits)
    gcms = gcms[isort, :]

    jsort = num.argsort(gcms[-1, :])[::-1]

    # ncols = 4
    # nrows = ((problem.ntargets + 1) - 1) / ncols + 1

    axes = fig.add_subplot(2, 2, 1)
    labelspace(axes)
    axes.set_ylabel('Relative contribution (smoothed)')
    axes.set_ylim(0.0, 1.0)

    axes2 = fig.add_subplot(2, 2, 3, sharex=axes)
    labelspace(axes2)
    axes2.set_xlabel('Tested model, sorted descending by global misfit value')

    axes2.set_ylabel('Square of misfit')

    axes2.set_ylim(0., 1.5)
    axes2.axhspan(1.0, 1.5, color=(0.8, 0.8, 0.8))
    axes2.set_yticks([0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5])
    axes2.set_yticklabels(
        ['0.0', '0.2', '0.4', '0.6', '0.8', '1', '10', '100', '1000', '10000',
         '100000'])

    axes2.set_xlim(imodels[0], imodels[-1])

    rel_ms_sum = num.zeros(model.nmodels)
    rel_ms_smooth_sum = num.zeros(model.nmodels)
    ms_smooth_sum = num.zeros(model.nmodels)
    b = num.hanning(100)
    b /= num.sum(b)
    a = [1]
    ii = 0

    for itarget in jsort:
        target = problem.targets[itarget]
        ms = gcms[:, itarget]
        ms = num.where(num.isfinite(ms), ms, 0.0)
        if num.all(ms == 0.0):
            continue

        rel_ms = ms / gms

        rel_ms_smooth = signal.filtfilt(b, a, rel_ms)

        ms_smooth = rel_ms_smooth * gms_softclip

        rel_poly_y = num.concatenate(
            [rel_ms_smooth_sum[::-1], rel_ms_smooth_sum + rel_ms_smooth])
        poly_x = num.concatenate([imodels[::-1], imodels])

        axes.fill(
            poly_x, rel_poly_y,
            alpha=0.5,
            color=colors[ii % len(colors)],
            label='%s (%.2g)' % (target.string_id(), num.mean(rel_ms[-1])))

        poly_y = num.concatenate(
            [ms_smooth_sum[::-1], ms_smooth_sum + ms_smooth])

        axes2.fill(poly_x, poly_y, alpha=0.5, color=colors[ii % len(colors)])

        rel_ms_sum += rel_ms

        # axes.plot(imodels, rel_ms_sum, color='black', alpha=0.1, zorder=-1)

        ms_smooth_sum += ms_smooth
        rel_ms_smooth_sum += rel_ms_smooth
        ii += 1

    axes.legend(
        title='Contributions (large to small at minimal global misfit)',
        bbox_to_anchor=(1.05, 0.0, 1.0, 1.0),
        loc='upper left',
        ncol=2, borderaxespad=0., prop={'size': 12})

    axes2.plot(imodels, gms_softclip, color='black')
    axes2.axhline(1.0, color=(0.5, 0.5, 0.5))
    fig.tight_layout()


def draw_bootstrap_figure(model, plt):

    fig = plt.figure()

    problem = model.problem
    gms = problem.global_misfits(model.misfits)

    imodels = num.arange(model.nmodels)

    axes = fig.add_subplot(1, 1, 1)

    gms_softclip = num.where(gms > 1.0, 0.1 * num.log10(gms) + 1.0, gms)

    ibests = []
    for ibootstrap in xrange(problem.nbootstrap):
        bms = problem.bootstrap_misfits(model.misfits, ibootstrap)
        isort_bms = num.argsort(bms)[::-1]

        ibests.append(isort_bms[-1])
        print num.argmin(bms), isort_bms[-1]

        bms_softclip = num.where(bms > 1.0, 0.1 * num.log10(bms) + 1.0, bms)
        axes.plot(imodels, bms_softclip[isort_bms], color='red', alpha=0.2)

    isort = num.argsort(gms)[::-1]
    iorder = num.empty(isort.size)
    iorder[isort] = imodels

    axes.plot(iorder[ibests], gms_softclip[ibests], 'x', color='black')

    m = num.median(gms[ibests])
    s = num.std(gms[ibests])

    axes.axhline(m+s, color='black', alpha=0.5)
    axes.axhline(m, color='black')
    axes.axhline(m-s, color='black', alpha=0.5)

    axes.plot(imodels, gms_softclip[isort], color='black')

    axes.set_xlim(imodels[0], imodels[-1])
    axes.set_xlabel('Tested model, sorted descending by global misfit value')


def gather(l, key, sort=None, filter=None):
    d = {}
    for x in l:
        if filter is not None and not filter(x):
            continue

        k = key(x)
        if k not in d:
            d[k] = []

        d[k].append(x)

    if sort is not None:
        for v in d.itervalues():
            v.sort(key=sort)

    return d


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
        (ma-mi) * space - (1.0 + space)
    t2 = num.concatenate((t, t[::-1]))
    axes.fill(
        t2, y2,
        clip_on=False,
        **kwargs)


def plot_dtrace_vline(axes, t, space, **kwargs):
    axes.plot([t, t], [-1.0 - space, -1.0], **kwargs)


def draw_fits_figures(ds, model, plt):
    fontsize = 10

    problem = model.problem

    for target in problem.targets:
        target.set_dataset(ds)

    target_index = dict(
        (target, i) for (i, target) in enumerate(problem.targets))

    gms = problem.global_misfits(model.misfits)
    isort = num.argsort(gms)
    gms = gms[isort]
    xs = model.xs[isort, :]
    misfits = model.misfits[isort, :]

    xbest = xs[0, :]

    ws = problem.get_target_weights()
    gcms = problem.global_contributions(misfits[:1])[0]

    w_max = num.nanmax(ws)
    gcm_max = num.nanmax(gcms)

    source = problem.unpack(xbest)

    target_to_result = {}
    all_syn_trs = []
    ms, ns, results = problem.evaluate(xbest, result_mode='full')

    dtraces = []
    for target, result in zip(problem.targets, results):
        if result is None:
            dtraces.append(None)
            continue

        itarget = target_index[target]
        w = target.get_combined_weight(problem.apply_balancing_weights)

        if target.misfit_config.domain == 'cc_max_norm':
            tref = (result.filtered_obs.tmin + result.filtered_obs.tmax) * 0.5
            for tr_filt, tr_proc, tshift in (
                    (result.filtered_obs,
                     result.processed_obs,
                     0.),
                    (result.filtered_syn,
                     result.processed_syn,
                     result.cc_shift)):

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

            dtrace = result.processed_syn.copy()
            dtrace.set_ydata(
                (
                    (result.processed_syn.get_ydata() -
                     result.processed_obs.get_ydata())**2))

        target_to_result[target] = result

        dtrace.meta = dict(super_group=target.super_group)
        dtraces.append(dtrace)

        result.processed_syn.meta = dict(super_group=target.super_group)
        all_syn_trs.append(result.processed_syn)

    if not all_syn_trs:
        logger.warn('no traces to show')
        return

    aminmaxs = trace.minmax(
        all_syn_trs,
        lambda tr: tr.meta['super_group'])

    dminmaxs = trace.minmax(
        [x for x in dtraces if x is not None],
        lambda tr: tr.meta['super_group'])

    for tr in dtraces:
        if tr:
            dmin, dmax = dminmaxs[tr.meta['super_group']]
            tr.ydata /= max(abs(dmin), abs(dmax))

    cg_to_targets = gather(
        problem.targets,
        lambda t: (t.super_group, t.group, t.codes[3]),
        filter=lambda t: t in target_to_result)

    cgs = sorted(cg_to_targets.keys())

    for cg in cgs:
        targets = cg_to_targets[cg]
        nframes = len(targets)

        nx = int(math.ceil(math.sqrt(nframes)))
        ny = (nframes-1)/nx+1

        nxmax = 4
        nymax = 4

        nxx = (nx-1) / nxmax + 1
        nyy = (ny-1) / nymax + 1

        # nz = nxx * nyy

        xs = num.arange(nx) / ((max(2, nx) - 1.0) / 2.)
        ys = num.arange(ny) / ((max(2, ny) - 1.0) / 2.)

        xs -= num.mean(xs)
        ys -= num.mean(ys)

        fxs = num.tile(xs, ny)
        fys = num.repeat(ys, nx)

        data = []

        for target in targets:
            azi = source.azibazi_to(target)[0]
            dist = source.distance_to(target)
            x = dist*num.sin(num.deg2rad(azi))
            y = dist*num.cos(num.deg2rad(azi))
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

        figures = {}
        for iy in xrange(ny):
            for ix in xrange(nx):
                if (iy, ix) not in frame_to_target:
                    continue

                ixx = ix/nxmax
                iyy = iy/nymax
                if (iyy, ixx) not in figures:
                    figures[iyy, ixx] = plt.figure(figsize=(16, 9))

                fig = figures[iyy, ixx]

                target = frame_to_target[iy, ix]

                amin, amax = aminmaxs[target.super_group]
                absmax = max(abs(amin), abs(amax))

                ny_this = min(ny, nymax)
                nx_this = min(nx, nxmax)
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

                itarget = target_index[target]
                result = target_to_result[target]

                dtrace = dtraces[itarget]

                tap_color_annot = (0.35, 0.35, 0.25)
                tap_color_edge = (0.85, 0.85, 0.80)
                tap_color_fill = (0.95, 0.95, 0.90)

                plot_taper(
                    axes2, result.processed_obs.get_xdata(), result.taper,
                    fc=tap_color_fill, ec=tap_color_edge)

                obs_color = scolor('aluminium5')
                obs_color_light = light(obs_color, 0.5)

                syn_color = scolor('scarletred2')
                syn_color_light = light(syn_color, 0.5)

                misfit_color = scolor('scarletred2')
                weight_color = scolor('chocolate2')

                cc_color = scolor('aluminium5')

                if target.misfit_config.domain == 'cc_max_norm':
                    tref = (result.filtered_obs.tmin +
                            result.filtered_obs.tmax) * 0.5

                    plot_dtrace(
                        axes2, dtrace, space, -1., 1.,
                        fc=light(cc_color, 0.5),
                        ec=cc_color)

                    plot_dtrace_vline(
                        axes2, tref, space, color=tap_color_annot)

                else:
                    plot_dtrace(
                        axes2, dtrace, space, 0., 1.,
                        fc=light(misfit_color, 0.5),
                        ec=misfit_color)

                plot_trace(
                    axes, result.filtered_syn,
                    color=syn_color_light, lw=1.5)

                plot_trace(
                    axes, result.filtered_obs,
                    color=obs_color_light)

                plot_trace(
                    axes, result.processed_syn,
                    color=syn_color, lw=1.5)

                plot_trace(
                    axes, result.processed_obs,
                    color=obs_color)

                tmarks = [
                    result.processed_obs.tmin,
                    result.processed_obs.tmax]

                for tmark in tmarks:
                    axes2.plot(
                        [tmark, tmark], [-0.9, 0.1], color=tap_color_annot)

                for tmark, text, ha in [
                        (tmarks[0], str_duration(tmarks[0] - source.time),
                         'right'),
                        (tmarks[1], '+' + str_duration2(tmarks[1] - tmarks[0]),
                         'left')]:

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

                rel_w = ws[itarget] / w_max
                rel_c = gcms[itarget] / gcm_max

                sw = 0.25
                sh = 0.1
                ph = 0.01

                for (ih, rw, facecolor, edgecolor) in [
                        (0, rel_w,  light(weight_color, 0.5), weight_color),
                        (1, rel_c,  light(misfit_color, 0.5), misfit_color)]:

                    bar = patches.Rectangle(
                        (1.0-rw*sw, 1.0-(ih+1)*sh+ph), rw*sw, sh-2*ph,
                        facecolor=facecolor, edgecolor=edgecolor,
                        zorder=10,
                        transform=axes.transAxes, clip_on=False)

                    axes.add_patch(bar)

                infos = []
                infos.append('.'.join(x for x in target.codes if x))
                dist = source.distance_to(target)
                azi = source.azibazi_to(target)[0]
                infos.append(str_dist(dist))
                infos.append(u'%.0f\u00B0' % azi)
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

        for (iyy, ixx), fig in figures.iteritems():
            title = '.'.join(x for x in cg if x)
            if len(figures) > 1:
                title += ' (%i/%i, %i/%i)' % (iyy+1, nyy, ixx+1, nxx)

            fig.suptitle(title, fontsize=fontsize)
            fig.savefig('fits-%s-%i-%i.pdf' % (
                title, iyy, ixx))

    plt.show()


def draw_hudson_figure(model, plt):

    color = 'black'
    fontsize = 12.
    markersize = fontsize * 1.5
    markersize_small = markersize * 0.2
    beachballsize = markersize
    beachballsize_small = beachballsize * 0.5
    width = 7.
    figsize = (width, width / (4./3.))

    problem = model.problem
    mean_source = core.get_mean_source(problem, model.xs)
    best_source = core.get_best_source(problem, model.xs, model.misfits)

    fig = plt.figure(figsize=figsize)
    axes = fig.add_subplot(1, 1, 1)

    data = []
    for ix, x in enumerate(model.xs):
        source = problem.unpack(x)
        mt = source.pyrocko_moment_tensor()
        u, v = hudson.project(mt)

        if random.random() < 0.1:
            try:
                beachball.plot_beachball_mpl(
                    mt, axes,
                    beachball_type='dc',
                    position=(u, v),
                    size=beachballsize_small,
                    color_t=color,
                    alpha=0.5,
                    zorder=1,
                    linewidth=0.25)
            except beachball.BeachballError, e:
                logger.warn(str(e))

        else:
            data.append((u, v))

    if data:
        u, v = num.array(data).T
        axes.plot(
            u, v, 'o',
            color=color,
            ms=markersize_small,
            mec='none',
            mew=0,
            alpha=0.25,
            zorder=0)

    hudson.draw_axes(axes)

    mt = mean_source.pyrocko_moment_tensor()
    u, v = hudson.project(mt)

    beachball.plot_beachball_mpl(
        mt, axes,
        beachball_type='dc',
        position=(u, v),
        size=beachballsize,
        color_t=color,
        zorder=2,
        linewidth=0.5)

    mt = best_source.pyrocko_moment_tensor()
    u, v = hudson.project(mt)

    axes.plot(
        u, v, 's',
        markersize=markersize,
        mew=1,
        mec='black',
        color='none',
        zorder=-2)

    mt = problem.base_source.pyrocko_moment_tensor()
    u, v = hudson.project(mt)

    beachball.plot_beachball_mpl(
        mt, axes,
        beachball_type='dc',
        position=(u, v),
        size=beachballsize,
        color_t='red',
        zorder=2,
        linewidth=0.5)

    fig.savefig('hudson.pdf')


def xpop(s, k):
    try:
        s.remove(k)
        return k

    except KeyError:
        return None


plot_dispatch = {
    'bootstrap': draw_bootstrap_figure,
    'sequence': draw_sequence_figures,
    'contributions': draw_contributions_figure,
    'jointpar': draw_jointpar_figures,
    'hudson': draw_hudson_figure,
    'fits': draw_fits_figures,
    'solution': draw_solution_figure}


def plot_result(dirname, plotnames_want):
    plotnames_want = set(plotnames_want)
    plotnames_avail = set(plot_dispatch.keys())

    unavailable = plotnames_want - plotnames_avail
    if unavailable:
        raise core.GrondError(
            'unavailable plotname: %s' % ', '.join(unavailable))

    mpl_init()

    if 3 != len({'bootstrap', 'sequence', 'contributions'} - plotnames_want):
        problem, xs, misfits = core.load_problem_info_and_data(
            dirname, subset=None)

        model = GrondModel()
        model.set_problem(problem)
        model.append(xs, misfits)

        for plotname in ['bootstrap', 'sequence', 'contributions']:
            if plotname in plotnames_want:
                plot_dispatch[plotname](model, plt)

    if 4 != len({'fits', 'jointpar', 'hudson', 'solution'} - plotnames_want):
        problem, xs, misfits = core.load_problem_info_and_data(
            dirname, subset='harvest')

        model = GrondModel()
        model.set_problem(problem)
        model.append(xs, misfits)

        for plotname in ['fits']:
            if plotname in plotnames_want:
                config = guts.load(filename=op.join(dirname, 'config.yaml'))
                config.set_basepath(dirname)
                problem.set_engine(config.engine_config.get_engine())
                ds = config.get_dataset()
                plot_dispatch[plotname](ds, model, plt)

        for plotname in ['jointpar', 'hudson', 'solution']:
            if plotname in plotnames_want:
                plot_dispatch[plotname](model, plt)

    plt.show()
