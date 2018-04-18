import logging
import math
import numpy as num

from pyrocko import gf, trace, weeding
from pyrocko.guts import (Object, String, Float, Bool, Int, StringChoice,
                          Timestamp, List)
from pyrocko.guts_array import Array

from grond.dataset import NotFound

from ..base import (MisfitConfig, MisfitTarget, MisfitResult,
                    TargetGroup, TargetAnalysisResult)

guts_prefix = 'grond'
logger = logging.getLogger('grond.targets.waveform.target')


class DomainChoice(StringChoice):
    choices = [
        'time_domain',
        'frequency_domain',
        'envelope',
        'absolute',
        'cc_max_norm']


class Trace(Object):
    pass


class WaveformMisfitConfig(MisfitConfig):
    fmin = Float.T()
    fmax = Float.T()
    ffactor = Float.T(default=1.5)
    tmin = gf.Timing.T(
        optional=True,
        help='Start of main time window used for waveform fitting.')
    tmax = gf.Timing.T(
        optional=True,
        help='End of main time window used for waveform fitting.')
    tfade = Float.T(
        optional=True,
        help='Decay time of taper prepended and appended to main time window '
             'used for waveform fitting [s].')
    pick_synthetic_traveltime = gf.Timing.T(
        optional=True,
        help='Synthetic phase arrival definition for alignment of observed '
             'and synthetic traces.')
    pick_phasename = String.T(
        optional=True,
        help='Name of picked phase for alignment of observed and synthetic '
             'traces.')
    domain = DomainChoice.T(
        default='time_domain',
        help='Type of data characteristic to be fitted.\n\nAvailable choices '
             'are: %s' % ', '.join("``'%s'``" % s
                                   for s in DomainChoice.choices))
    norm_exponent = Int.T(
        default=2,
        help='Exponent to use in norm (1: L1-norm, 2: L2-norm)')
    tautoshift_max = Float.T(
        default=0.0,
        help='If non-zero, allow synthetic and observed traces to be shifted '
             'against each other by up to +/- the given value [s].')
    autoshift_penalty_max = Float.T(
        default=0.0,
        help='If non-zero, a penalty misfit is added for non-zero shift '
             'values.\n\nThe penalty value is computed as '
             '``autoshift_penalty_max * normalization_factor * tautoshift**2 '
             '/ tautoshift_max**2``')

    ranges = {}

    def get_full_frequency_range(self):
        return self.fmin / self.ffactor, self.fmax * self.ffactor


def log_exclude(target, reason):
    logger.debug('excluding potential target %s: %s' % (
        target.string_id(), reason))


class WaveformTargetGroup(TargetGroup):
    distance_min = Float.T(optional=True)
    distance_max = Float.T(optional=True)
    distance_3d_min = Float.T(optional=True)
    distance_3d_max = Float.T(optional=True)
    depth_min = Float.T(optional=True)
    depth_max = Float.T(optional=True)
    limit = Int.T(optional=True)
    channels = List.T(String.T(), optional=True)
    misfit_config = WaveformMisfitConfig.T()

    def get_targets(self, ds, event, default_path):
        logger.debug('Selecting waveform targets...')
        origin = event
        targets = []

        for st in ds.get_stations():
            for cha in self.channels:

                nslc = st.nsl() + (cha,)

                target = WaveformMisfitTarget(
                    quantity='displacement',
                    codes=nslc,
                    lat=st.lat,
                    lon=st.lon,
                    depth=st.depth,
                    interpolation=self.interpolation,
                    store_id=self.store_id,
                    misfit_config=self.misfit_config,
                    manual_weight=self.weight,
                    normalisation_family=self.normalisation_family,
                    path=self.path or default_path)

                if ds.is_blacklisted((st.nsl() + (cha,))):
                    log_exclude(target, 'blacklisted')
                    continue

                if self.distance_min is not None and \
                   target.distance_to(origin) < self.distance_min:
                    log_exclude(target, 'distance < distance_min')
                    continue

                if self.distance_max is not None and \
                   target.distance_to(origin) > self.distance_max:
                    log_exclude(target, 'distance > distance_max')
                    continue

                if self.distance_3d_min is not None and \
                   target.distance_3d_to(origin) < self.distance_3d_min:
                    log_exclude(target, 'distance_3d < distance_3d_min')
                    continue

                if self.distance_3d_max is not None and \
                   target.distance_3d_to(origin) > self.distance_3d_max:
                    log_exclude(target, 'distance_3d > distance_3d_max')
                    continue

                if self.depth_min is not None and \
                   target.depth < self.depth_min:
                    log_exclude(target, 'depth < depth_min')
                    continue

                if self.depth_max is not None and \
                   target.depth > self.depth_max:
                    log_exclude(target, 'depth > depth_max')
                    continue

                azi, _ = target.azibazi_to(origin)
                if cha == 'R':
                    target.azimuth = azi - 180.
                    target.dip = 0.
                elif cha == 'T':
                    target.azimuth = azi - 90.
                    target.dip = 0.
                elif cha == 'Z':
                    target.azimuth = 0.
                    target.dip = -90.

                target.set_dataset(ds)
                targets.append(target)

        if self.limit:
            return weed(origin, targets, self.limit)[0]
        else:
            return targets


class TraceSpectrum(Object):
    network = String.T()
    station = String.T()
    location = String.T()
    channel = String.T()
    deltaf = Float.T(default=1.0)
    fmin = Float.T(default=0.0)
    ydata = Array.T(shape=(None,), dtype=num.complex, serialize_as='list')

    def get_ydata(self):
        return self.ydata

    def get_xdata(self):
        return self.fmin + num.arange(self.ydata.size) * self.deltaf


class WaveformMisfitResult(gf.Result, MisfitResult):
    processed_obs = Trace.T(optional=True)
    processed_syn = Trace.T(optional=True)
    filtered_obs = Trace.T(optional=True)
    filtered_syn = Trace.T(optional=True)
    spectrum_obs = TraceSpectrum.T(optional=True)
    spectrum_syn = TraceSpectrum.T(optional=True)

    taper = trace.Taper.T(optional=True)
    tobs_shift = Float.T(optional=True)
    tsyn_pick = Timestamp.T(optional=True)
    tshift = Float.T(optional=True)
    cc = Trace.T(optional=True)


class WaveformMisfitTarget(gf.Target, MisfitTarget):
    flip_norm = Bool.T(default=False)
    misfit_config = WaveformMisfitConfig.T()

    def __init__(self, **kwargs):
        gf.Target.__init__(self, **kwargs)
        MisfitTarget.__init__(self, **kwargs)

    def string_id(self):
        return '.'.join(x for x in (self.path,) + self.codes if x)

    def get_combined_weight(self, apply_balancing_weights):
        w = self.manual_weight
        if apply_balancing_weights:
            w *= self.get_balancing_weight()
        return num.array([w], dtype=num.float)

    def get_balancing_weight(self):
        if not self.analysis_result:
            raise TargetAnalysisResult.tNoResults(
                'no balancing weights available')

        return self.analysis_result.balancing_weight

    def get_taper_params(self, engine, source):
        store = engine.get_store(self.store_id)
        config = self.misfit_config
        tmin_fit = source.time + store.t(config.tmin, source, self)
        tmax_fit = source.time + store.t(config.tmax, source, self)
        tfade = 1.0/config.fmin
        if config.tfade is None:
            tfade_taper = tfade
        else:
            tfade_taper = config.tfade

        return tmin_fit, tmax_fit, tfade, tfade_taper

    def get_backazimuth_for_waveform(self):
        return backazimuth_for_waveform(self.azimuth, self.codes)

    def get_freqlimits(self):
        config = self.misfit_config

        return (
            config.fmin/config.ffactor,
            config.fmin, config.fmax,
            config.fmax*config.ffactor)

    def get_pick_shift(self, engine, source):
        config = self.misfit_config
        tobs = None
        tsyn = None
        ds = self.get_dataset()

        if config.pick_synthetic_traveltime and config.pick_phasename:
            store = engine.get_store(self.store_id)
            tsyn = source.time + store.t(
                config.pick_synthetic_traveltime, source, self)

            marker = ds.get_pick(
                source.name,
                self.codes[:3],
                config.pick_phasename)

            if marker:
                tobs = marker.tmin

        return tobs, tsyn

    def get_cutout_timespan(self, tmin, tmax, tfade):
        tinc_obs = 1.0 / self.misfit_config.fmin

        tmin_obs = (math.floor(
            (tmin - tfade) / tinc_obs) - 1.0) * tinc_obs
        tmax_obs = (math.ceil(
            (tmax + tfade) / tinc_obs) + 1.0) * tinc_obs

        return tmin_obs, tmax_obs

    def post_process(self, engine, source, tr_syn):

        tr_syn = tr_syn.pyrocko_trace()
        nslc = self.codes

        config = self.misfit_config

        tmin_fit, tmax_fit, tfade, tfade_taper = \
            self.get_taper_params(engine, source)

        ds = self.get_dataset()

        tobs, tsyn = self.get_pick_shift(engine, source)
        if None not in (tobs, tsyn):
            tobs_shift = tobs - tsyn
        else:
            tobs_shift = 0.0

        tr_syn.extend(
            tmin_fit - tfade * 2.0,
            tmax_fit + tfade * 2.0,
            fillmethod='repeat')

        freqlimits = self.get_freqlimits()

        tr_syn = tr_syn.transfer(
            freqlimits=freqlimits,
            tfade=tfade)

        tr_syn.chop(tmin_fit - 2*tfade, tmax_fit + 2*tfade)

        tmin_obs, tmax_obs = self.get_cutout_timespan(
            tmin_fit+tobs_shift, tmax_fit+tobs_shift, tfade)

        try:
            tr_obs = ds.get_waveform(
                nslc,
                tinc_cache=1.0/config.fmin,
                tmin=tmin_fit+tobs_shift-tfade,
                tmax=tmax_fit+tobs_shift+tfade,
                tfade=tfade,
                freqlimits=freqlimits,
                deltat=tr_syn.deltat,
                cache=True,
                backazimuth=self.get_backazimuth_for_waveform())

            if tobs_shift != 0.0:
                tr_obs = tr_obs.copy()
                tr_obs.shift(-tobs_shift)

            mr = misfit(
                tr_obs, tr_syn,
                taper=trace.CosTaper(
                    tmin_fit - tfade_taper,
                    tmin_fit,
                    tmax_fit,
                    tmax_fit + tfade_taper),
                domain=config.domain,
                exponent=config.norm_exponent,
                flip=self.flip_norm,
                result_mode=self._result_mode,
                tautoshift_max=config.tautoshift_max,
                autoshift_penalty_max=config.autoshift_penalty_max)

            mr.tobs_shift = float(tobs_shift)
            mr.tsyn_pick = float_or_none(tsyn)

            return mr

        except NotFound as e:
            logger.debug(str(e))
            raise gf.SeismosizerError('no waveform data, %s' % str(e))

    def prepare_modelling(self, engine, source):
        return [self]

    def finalize_modelling(
            self, engine, source, modelling_targets, modelling_results):

        return modelling_results[0]

    def get_plain_targets(self, engine, source):
        d = dict(
            (k, getattr(self, k)) for k in gf.Target.T.propnames)
        return [gf.Target(**d)]

    @classmethod
    def get_plotter_class(cls):
        from . import plot
        return plot.WaveformTargetPlotter


def misfit(
        tr_obs, tr_syn, taper, domain, exponent, tautoshift_max,
        autoshift_penalty_max, flip, result_mode='sparse'):

    '''
    Calculate misfit between observed and synthetic trace.

    :param tr_obs: observed trace as :py:class:`pyrocko.trace.Trace`
    :param tr_syn: synthetic trace as :py:class:`pyrocko.trace.Trace`
    :param taper: taper applied in timedomain as
        :py:class:`pyrocko.trace.Taper`
    :param domain: how to calculate difference, see :py:class:`DomainChoice`
    :param exponent: exponent of Lx type norms
    :param tautoshift_max: if non-zero, return lowest misfit when traces are
        allowed to shift against each other by up to +/- ``tautoshift_max``
    :param autoshift_penalty_max: if non-zero, a penalty misfit is added for
        for non-zero shift values. The penalty value is
        ``autoshift_penalty_max * normalization_factor * \
tautoshift**2 / tautoshift_max**2``
    :param flip: ``bool``, if set to ``True``, normalization factor is
        computed against *tr_syn* rather than *tr_obs*
    :param result_mode: ``'full'``, include traces and spectra or ``'sparse'``,
        include only misfit and normalization factor in result

    :returns: object of type :py:class:`WaveformMisfitResult`
    '''

    trace.assert_same_sampling_rate(tr_obs, tr_syn)
    deltat = tr_obs.deltat
    tmin, tmax = taper.time_span()

    tr_proc_obs, trspec_proc_obs = _process(tr_obs, tmin, tmax, taper, domain)
    tr_proc_syn, trspec_proc_syn = _process(tr_syn, tmin, tmax, taper, domain)

    tshift = None
    ctr = None
    deltat = tr_proc_obs.deltat
    if domain in ('time_domain', 'envelope', 'absolute'):
        a, b = tr_proc_syn.ydata, tr_proc_obs.ydata
        if flip:
            b, a = a, b

        nshift_max = max(0, min(a.size-1,
                                int(math.floor(tautoshift_max / deltat))))

        if nshift_max == 0:
            m, n = trace.Lx_norm(a, b, norm=exponent)
        else:
            mns = []
            for ishift in range(-nshift_max, nshift_max+1):
                if ishift < 0:
                    a_cut = a[-ishift:]
                    b_cut = b[:ishift]
                elif ishift == 0:
                    a_cut = a
                    b_cut = b
                elif ishift > 0:
                    a_cut = a[:-ishift]
                    b_cut = b[ishift:]

                mns.append(trace.Lx_norm(a_cut, b_cut, norm=exponent))

            ms, ns = num.array(mns).T

            iarg = num.argmin(ms)
            tshift = (iarg-nshift_max)*deltat

            m, n = ms[iarg], ns[iarg]
            m += autoshift_penalty_max * n * tshift**2 / tautoshift_max**2

    elif domain == 'cc_max_norm':

        ctr = trace.correlate(
            tr_proc_syn,
            tr_proc_obs,
            mode='same',
            normalization='normal')

        tshift, cc_max = ctr.max()
        m = 0.5 - 0.5 * cc_max
        n = 0.5

    elif domain == 'frequency_domain':
        a, b = trspec_proc_syn.ydata, trspec_proc_obs.ydata
        if flip:
            b, a = a, b

        m, n = trace.Lx_norm(num.abs(a), num.abs(b), norm=exponent)

    if result_mode == 'full':
        result = WaveformMisfitResult(
            misfits=num.array([[m, n]], dtype=num.float),
            processed_obs=tr_proc_obs,
            processed_syn=tr_proc_syn,
            filtered_obs=tr_obs.copy(),
            filtered_syn=tr_syn,
            spectrum_obs=trspec_proc_obs,
            spectrum_syn=trspec_proc_syn,
            taper=taper,
            tshift=tshift,
            cc=ctr)

    elif result_mode == 'sparse':
        result = WaveformMisfitResult(
            misfits=num.array([[m, n]], dtype=num.float))
    else:
        assert False

    return result


def _extend_extract(tr, tmin, tmax):
    deltat = tr.deltat
    itmin_frame = int(math.floor(tmin/deltat))
    itmax_frame = int(math.ceil(tmax/deltat))
    nframe = itmax_frame - itmin_frame + 1
    n = tr.data_len()
    a = num.empty(nframe, dtype=num.float)
    itmin_tr = int(round(tr.tmin / deltat))
    itmax_tr = itmin_tr + n
    icut1 = min(max(0, itmin_tr - itmin_frame), nframe)
    icut2 = min(max(0, itmax_tr - itmin_frame), nframe)
    icut1_tr = min(max(0, icut1 + itmin_frame - itmin_tr), n)
    icut2_tr = min(max(0, icut2 + itmin_frame - itmin_tr), n)
    a[:icut1] = tr.ydata[0]
    a[icut1:icut2] = tr.ydata[icut1_tr:icut2_tr]
    a[icut2:] = tr.ydata[-1]
    tr = tr.copy(data=False)
    tr.tmin = itmin_frame * deltat
    tr.set_ydata(a)
    return tr


def _process(tr, tmin, tmax, taper, domain):
    tr_proc = _extend_extract(tr, tmin, tmax)
    tr_proc.taper(taper)

    df = None
    trspec_proc = None

    if domain == 'envelope':
        tr_proc = tr_proc.envelope(inplace=False)
        tr_proc.set_ydata(num.abs(tr_proc.get_ydata()))

    elif domain == 'absolute':
        tr_proc.set_ydata(num.abs(tr_proc.get_ydata()))

    elif domain == 'frequency_domain':
        ndata = tr_proc.ydata.size
        nfft = trace.nextpow2(ndata)
        padded = num.zeros(nfft, dtype=num.float)
        padded[:ndata] = tr_proc.ydata
        spectrum = num.fft.rfft(padded)
        df = 1.0 / (tr_proc.deltat * nfft)

        trspec_proc = TraceSpectrum(
            network=tr_proc.network,
            station=tr_proc.station,
            location=tr_proc.location,
            channel=tr_proc.channel,
            deltaf=df,
            fmin=0.0,
            ydata=spectrum)

    return tr_proc, trspec_proc


def backazimuth_for_waveform(azimuth, nslc):
    if nslc[-1] == 'R':
        backazimuth = azimuth + 180.
    elif nslc[-1] == 'T':
        backazimuth = azimuth + 90.
    else:
        backazimuth = None

    return backazimuth


def float_or_none(x):
    if x is None:
        return x
    else:
        return float(x)


def weed(origin, targets, limit, neighborhood=3):
    azimuths = num.zeros(len(targets))
    dists = num.zeros(len(targets))
    for i, target in enumerate(targets):
        _, azimuths[i] = target.azibazi_to(origin)
        dists[i] = target.distance_to(origin)

    badnesses = num.ones(len(targets), dtype=float)
    deleted, meandists_kept = weeding.weed(
        azimuths, dists, badnesses,
        nwanted=limit,
        neighborhood=neighborhood)

    targets_weeded = [
        target for (delete, target) in zip(deleted, targets) if not delete]

    return targets_weeded, meandists_kept, deleted


__all__ = '''
    WaveformTargetGroup
    WaveformMisfitConfig
    WaveformMisfitTarget
    WaveformMisfitResult
'''.split()
