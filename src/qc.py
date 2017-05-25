import logging
from collections import defaultdict

import numpy as num

from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid

from pyrocko import gf, trace, orthodrome as od, plot
from grond import dataset

km = 1000.

logger = logging.getLogger('grond.qc')


def polarization(ds, store, timing, fmin, fmax, ffactor=1.5, tfactor=2.):
    event = ds.get_event()
    stations = ds.get_stations()

    nsl_to_station = dict(
        (s.nsl(), s) for s in stations)

    source = gf.Source.from_pyrocko_event(event)

    trs = []
    for station in stations:

        for component in 'ZNE':

            ttp = store.t(timing, source, station)

            tmin = event.time + ttp - tfactor / fmin
            tmax = event.time + ttp + tfactor / fmin

            nslc = station.nsl() + (component,)

            freqlimits = [
                fmin / ffactor,
                fmin,
                fmax,
                fmax * ffactor]

            tfade = 1.0 / (fmin / ffactor)

            tlong = 5.0 / fmin
            tshort = 0.2 * tlong

            try:
                trs_projected, trs_restituted, trs_raw = \
                    ds.get_waveform(
                        nslc,
                        tmin=tmin,
                        tmax=tmax,
                        tfade=tfade,
                        freqlimits=freqlimits,
                        debug=True)

                # for tr in trs_projected:
                #     tr.sta_lta_right(tshort, tlong)

                trs.extend(trs_projected)

            except dataset.NotFound, e:
                logger.warn(str(e))
                continue

    d_trs = defaultdict(dict)
    for tr in trs:
        d_trs[tr.nslc_id[:3]][tr.nslc_id[3]] = tr

    fontsize = 11.
    plot.mpl_init(fontsize=fontsize)
    fig = plt.figure(figsize=(8., 8.))
    plot.mpl_margins(fig, w=7., h=6., units=fontsize)

    grid = ImageGrid(
        fig, 111, nrows_ncols=(2, 2),
        axes_pad=0.5,
        add_all=True,
        label_mode='L',
        aspect=True)

    print dir(grid)

    axes_en = grid[0]
    axes_en.set_ylabel('Northing [km]')

    axes_dn = grid[1]
    axes_dn.locator_params(axis='x', nbins=4)
    axes_dn.set_xlabel('Depth [km]')

    axes_ed = grid[2]
    axes_ed.locator_params(axis='y', nbins=4)
    axes_ed.set_ylabel('Depth [km]')
    axes_ed.set_xlabel('Easting [km]')

    axes_en.plot(0., 0., '*')
    axes_dn.plot(event.depth/km, 0., '*')
    axes_ed.plot(0., event.depth/km, '*')

    grid[3].set_axis_off()

    factor = 5000.

    for nsl in sorted(d_trs.keys()):

        tr_e = d_trs[nsl].get('E', None)
        tr_n = d_trs[nsl].get('N', None)
        tr_z = d_trs[nsl].get('Z', None)

        station = nsl_to_station[nsl]

        n, e = od.latlon_to_ne(
            event.lat, event.lon, station.lat, station.lon)

        d = station.depth

        axes_en.plot(e/km, n/km, '^')
        axes_dn.plot(d/km, n/km, '^')
        axes_ed.plot(e/km, d/km, '^')

        arr_e = tr_e.ydata if tr_e else None
        arr_n = tr_n.ydata if tr_n else None
        arr_z = tr_z.ydata if tr_z else None

        amax = num.max(num.abs(num.sqrt(arr_e**2 + arr_n**2 + arr_z**2)))
        axes_en.plot(
            (e + arr_e/amax * factor)/km, (n + arr_n/amax * factor)/km,
            color=plot.mpl_color('skyblue2'))
        axes_dn.plot(
            (d - arr_z/amax * factor)/km, (n + arr_n/amax * factor)/km,
            color=plot.mpl_color('skyblue2'))
        axes_ed.plot(
            (e + arr_e/amax * factor)/km, (d - arr_z/amax * factor)/km,
            color=plot.mpl_color('skyblue2'))

    axes_ed.invert_yaxis()

    plt.show()
    # trace.snuffle(trs)
