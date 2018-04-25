import logging
from collections import defaultdict

import numpy as num

from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid

from pyrocko import gf, orthodrome as od, plot, model, trace
from grond import dataset

km = 1000.

logger = logging.getLogger('grond.qc')


def darken(c):
    return (c[0]*0.5, c[1]*0.5, c[2]*0.5)


def plot_color_line(axes, x, y, t, color, tmin, tmax):
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.collections import LineCollection

    cmap = LinearSegmentedColormap.from_list(
        'noname', [color, (1., 1., 1.)], 100)

    points = num.array([x, y], dtype=num.float).T.reshape(-1, 1, 2)
    segments = num.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap, norm=plt.Normalize(tmin, tmax))
    lc.set_array(t)
    axes.add_collection(lc, autolim=True)


def polarization(
        ds, store, timing, fmin, fmax, ffactor,
        time_factor_pre=2.,
        time_factor_post=2.,
        distance_min=None,
        distance_max=None,
        depth_min=None,
        depth_max=None,
        size_factor=0.05,
        nsl_to_time=None,
        output_filename=None,
        output_format=None,
        output_dpi=None):

    event = ds.get_event()
    stations = ds.get_stations()

    source = gf.Source.from_pyrocko_event(event)

    trs = []
    for station in stations:

        nsl = station.nsl()

        dist = source.distance_to(station)

        if distance_min is not None and dist < distance_min:
            continue

        if distance_max is not None and distance_max < dist:
            continue

        if depth_min is not None and station.depth < depth_min:
            continue

        if depth_max is not None and depth_max < station.depth:
            continue

        if nsl_to_time is None:
            tp = event.time + store.t(timing, source, station)

        else:
            if nsl not in nsl_to_time:
                continue

            tp = nsl_to_time[nsl]

        for component in 'ZNE':

            tmin = tp - time_factor_pre / fmin
            tmax = tp + time_factor_post / fmin

            nslc = nsl + (component,)

            freqlimits = [
                fmin / ffactor,
                fmin,
                fmax,
                fmax * ffactor]

            tfade = 1.0 / (fmin / ffactor)

            try:
                trs_projected, trs_restituted, trs_raw, _ = \
                    ds.get_waveform(
                        nslc,
                        tmin=tmin,
                        tmax=tmax,
                        tfade=tfade,
                        freqlimits=freqlimits,
                        debug=True)

                for tr in trs_projected:
                    tr.shift(-tp)

                trs.extend(trs_projected)

            except dataset.NotFound as e:
                logger.warn(str(e))
                continue

    trace.snuffle(trs, stations=stations)

    plot_polarizations(
        stations, trs,
        event=event,
        size_factor=size_factor,
        output_filename=output_filename,
        output_format=output_format,
        output_dpi=output_dpi)


def plot_polarizations(
        stations, trs,
        event=None,
        size_factor=0.05,
        fontsize=10.,
        output_filename=None,
        output_format=None,
        output_dpi=None):

    if event is None:
        slats = num.array([s.lat for s in stations], dtype=num.float)
        slons = num.array([s.lon for s in stations], dtype=num.float)
        clat, clon = od.geographic_midpoint(slats, slons)
        event = od.Loc(clat, clon)

    nsl_c_to_trs = defaultdict(dict)
    for tr in trs:
        nsl_c_to_trs[tr.nslc_id[:3]][tr.nslc_id[3]] = tr

    nsl_to_station = dict(
        (s.nsl(), s) for s in stations)

    plot.mpl_init(fontsize=fontsize)
    fig = plt.figure(figsize=plot.mpl_papersize('a4', 'landscape'))
    plot.mpl_margins(fig, w=7., h=6., units=fontsize)

    grid = ImageGrid(
        fig, 111, nrows_ncols=(2, 2),
        axes_pad=0.5,
        add_all=True,
        label_mode='L',
        aspect=True)

    axes_en = grid[0]
    axes_en.set_ylabel('Northing [km]')

    axes_dn = grid[1]
    axes_dn.locator_params(axis='x', nbins=4)
    axes_dn.set_xlabel('Depth [km]')

    axes_ed = grid[2]
    axes_ed.locator_params(axis='y', nbins=4)
    axes_ed.set_ylabel('Depth [km]')
    axes_ed.set_xlabel('Easting [km]')

    if isinstance(event, model.Event):
        axes_en.plot(0., 0., '*')
        axes_dn.plot(event.depth/km, 0., '*')
        axes_ed.plot(0., event.depth/km, '*')

    grid[3].set_axis_off()

    locations = []
    for nsl in sorted(nsl_c_to_trs.keys()):
        station = nsl_to_station[nsl]
        n, e = od.latlon_to_ne(
            event.lat, event.lon, station.lat, station.lon)

        locations.append((n, e))

    ns, es = num.array(locations, dtype=num.float).T

    n_min = num.min(ns)
    n_max = num.max(ns)
    e_min = num.min(es)
    e_max = num.max(es)

    factor = max((n_max - n_min) * size_factor, (e_max - e_min) * size_factor)

    fontsize_annot = fontsize * 0.7

    data = {}
    for insl, nsl in enumerate(sorted(nsl_c_to_trs.keys())):

        color = plot.mpl_graph_color(insl)

        try:
            tr_e = nsl_c_to_trs[nsl]['E']
            tr_n = nsl_c_to_trs[nsl]['N']
            tr_z = nsl_c_to_trs[nsl]['Z']

        except KeyError:
            continue

        station = nsl_to_station[nsl]

        n, e = od.latlon_to_ne(
            event.lat, event.lon, station.lat, station.lon)

        d = station.depth

        axes_en.annotate(
            '.'.join(x for x in nsl if x),
            xy=(e/km, n/km),
            xycoords='data',
            xytext=(fontsize_annot/3., fontsize_annot/3.),
            textcoords='offset points',
            verticalalignment='bottom',
            horizontalalignment='left',
            rotation=0.,
            size=fontsize_annot)

        axes_en.plot(e/km, n/km, '^', mfc=color, mec=darken(color))
        axes_dn.plot(d/km, n/km, '^', mfc=color, mec=darken(color))
        axes_ed.plot(e/km, d/km, '^', mfc=color, mec=darken(color))

        arr_e = tr_e.ydata
        arr_n = tr_n.ydata
        arr_z = tr_z.ydata
        arr_t = tr_z.get_xdata()

        data[nsl] = (arr_e, arr_n, arr_z, arr_t, n, e, d, color)

    amaxs = []
    amax_hors = []
    for nsl in sorted(data.keys()):
        arr_e, arr_n, arr_z, arr_t, n, e, d, color = data[nsl]
        amaxs.append(
            num.max(num.abs(num.sqrt(arr_e**2 + arr_n**2 + arr_z**2))))
        amax_hors.append(
            num.max(num.abs(num.sqrt(arr_e**2 + arr_n**2))))

    amax = num.median(amaxs)
    amax_hor = num.median(amax_hors)

    for nsl in sorted(data.keys()):
        arr_e, arr_n, arr_z, arr_t, n, e, d, color = data[nsl]
        tmin = arr_t.min()
        tmax = arr_t.max()
        plot_color_line(
            axes_en,
            (e + arr_e/amax_hor * factor)/km, (n + arr_n/amax_hor * factor)/km,
            arr_t, color, tmin, tmax)
        plot_color_line(
            axes_dn,
            (d - arr_z/amax * factor)/km, (n + arr_n/amax * factor)/km,
            arr_t, color, tmin, tmax)
        plot_color_line(
            axes_ed,
            (e + arr_e/amax * factor)/km, (d - arr_z/amax * factor)/km,
            arr_t, color, tmin, tmax)

    axes_ed.invert_yaxis()

    for axes in (axes_dn, axes_ed, axes_en):
        axes.autoscale_view(tight=True)

    if output_filename is None:
        plt.show()
    else:
        fig.savefig(output_filename, format=output_format, dpi=output_dpi)
