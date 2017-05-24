import logging

from pyrocko import gf, trace
from grond import dataset


logger = logging.getLogger('grond.qc')



def polarization(ds, store, timing, fmin, fmax, ffactor=1.5, tfactor=10.):
    event = ds.get_event()
    stations = ds.get_stations()

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

            try:
                trs_projected, trs_restituted, trs_raw = \
                    ds.get_waveform(
                        nslc,
                        tmin=tmin,
                        tmax=tmax,
                        tfade=tfade,
                        freqlimits=freqlimits,
                        debug=True)

                trs.extend(trs_projected)

            except dataset.NotFound, e:
                logger.warn(str(e))
                continue

    trace.snuffle(trs)


