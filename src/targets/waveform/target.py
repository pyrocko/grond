from __future__ import print_function

import logging
import math
import numpy as num

from pyrocko import gf, trace, weeding, util
from pyrocko.guts import (Object, String, Float, Bool, Int, StringChoice,
                          Timestamp, List, Dict)
from pyrocko.guts_array import Array

from grond.dataset import NotFound
from grond.meta import GrondError, nslcs_to_patterns

from ..base import (MisfitConfig, MisfitTarget, MisfitResult, TargetGroup)
from grond.meta import has_get_plot_classes

from pyrocko import crust2x2
from string import Template

guts_prefix = 'grond'
logger = logging.getLogger('grond.targets.waveform.target')


class StoreIDSelectorError(GrondError):
    pass


class StoreIDSelector(Object):
    '''
    Base class for GF store selectors.

    GF store selectors can be implemented to select different stores, based on
    station location, source location or other characteristics.
    '''

    pass


class Crust2StoreIDSelector(StoreIDSelector):
    '''
    Store ID selector picking CRUST 2.0 model based on event location.
    '''

    template = String.T(
        help="Template for the GF store ID. For example ``'crust2_${id}'`` "
             "where ``'${id}'`` will be replaced with the corresponding CRUST "
             "2.0 profile identifier for the source location.")

    def get_store_id(self, event, st, cha):
        s = Template(self.template)
        return s.substitute(id=(
            crust2x2.get_profile(event.lat, event.lon)._ident).lower())


class StationDictStoreIDSelector(StoreIDSelector):
    '''
    Store ID selector using a manual station to store ID mapping.
    '''

    mapping = Dict.T(
        String.T(), gf.StringID.T(),
        help='Dictionary with station to store ID pairs, keys are NET.STA. '
             "Add a fallback store ID under the key ``'others'``.")

    def get_store_id(self, event, st, cha):
        try:
            store_id = self.mapping['%s.%s' % (st.network, st.station)]
        except KeyError:
            try:
                store_id = self.mapping['others']
            except KeyError:
                raise StoreIDSelectorError(
                    'No store ID found for station "%s.%s".' % (
                        st.network, st.station))

        return store_id


class DepthRangeToStoreID(Object):
    depth_min = Float.T()
    depth_max = Float.T()
    store_id = gf.StringID.T()


class StationDepthStoreIDSelector(StoreIDSelector):
    '''
    Store ID selector using a mapping from station depth range to store ID.
    '''

    depth_ranges = List.T(DepthRangeToStoreID.T())

    def get_store_id(self, event, st, cha):
        for r in self.depth_ranges:
            if r.depth_min <= st.depth < r.depth_max:
                return r.store_id

        raise StoreIDSelectorError(
            'No store ID found for station "%s.%s" at %g m depth.' % (
                st.network, st.station, st.depth))


class DomainChoice(StringChoice):
    choices = [
        'time_domain',
        'frequency_domain',
        'log_frequency_domain',
        'envelope',
        'absolute',
        'cc_max_norm']


class WaveformMisfitConfig(MisfitConfig):
    quantity = gf.QuantityType.T(default='displacement')
    fmin = Float.T(default=0.0, help='minimum frequency of bandpass filter')
    fmax = Float.T(help='maximum frequency of bandpass filter')
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
    logger.debug('Excluding potential target %s: %s' % (
        target.string_id(), reason))


class WaveformTargetGroup(TargetGroup):
    '''Handles seismogram targets or other targets of dynamic ground motion.
    '''
    distance_min = Float.T(
        optional=True,
        help='excludes targets nearer to source, along a great circle')
    distance_max = Float.T(
        optional=True,
        help='excludes targets farther from source, along a great circle')
    distance_3d_min = Float.T(
        optional=True,
        help='excludes targets nearer from source (direct distance)')
    distance_3d_max = Float.T(
        optional=True,
        help='excludes targets farther from source (direct distance)')
    depth_min = Float.T(
        optional=True,
        help='excludes targets with smaller depths')
    depth_max = Float.T(
        optional=True,
        help='excludes targets with larger depths')
    include = List.T(
        String.T(),
        optional=True,
        help='If not None, list of stations/components to include according '
             'to their STA, NET.STA, NET.STA.LOC, or NET.STA.LOC.CHA codes.')
    exclude = List.T(
        String.T(),
        help='Stations/components to be excluded according to their STA, '
             'NET.STA, NET.STA.LOC, or NET.STA.LOC.CHA codes.')
    limit = Int.T(optional=True)
    channels = List.T(
        String.T(),
        optional=True,
        help="set channels to include, e.g. ['Z', 'T']")
    misfit_config = WaveformMisfitConfig.T()
    store_id_selector = StoreIDSelector.T(
        optional=True,
        help='select GF store based on event-station geometry.')

    def get_targets(self, ds, event, default_path='none'):
        logger.debug('Selecting waveform targets...')
        origin = event
        targets = []

        stations = ds.get_stations()
        if len(stations) == 0:
            logger.warning(
                'No stations found to create waveform target group.')

        for st in ds.get_stations():
            logger.debug('Selecting waveforms for station %s.%s.%s' % st.nsl())
            for cha in self.channels:
                nslc = st.nsl() + (cha,)

                logger.debug('Selecting waveforms for %s.%s.%s.%s' % nslc)

                if self.store_id_selector:
                    store_id = self.store_id_selector.get_store_id(
                        event, st, cha)
                else:
                    store_id = self.store_id

                logger.debug('Selecting waveforms for %s.%s.%s.%s' % nslc)

                target = WaveformMisfitTarget(
                    quantity='displacement',
                    codes=nslc,
                    lat=st.lat,
                    lon=st.lon,
                    north_shift=st.north_shift,
                    east_shift=st.east_shift,
                    depth=st.depth,
                    interpolation=self.interpolation,
                    store_id=store_id,
                    misfit_config=self.misfit_config,
                    manual_weight=self.weight,
                    normalisation_family=self.normalisation_family,
                    path=self.path or default_path)

                if ds.is_blacklisted(nslc):
                    log_exclude(target, 'excluded by dataset')
                    continue

                if util.match_nslc(
                        nslcs_to_patterns(self.exclude), nslc):
                    log_exclude(target, 'excluded by target group')
                    continue

                if self.include is not None and not util.match_nslc(
                        nslcs_to_patterns(self.include), nslc):
                    log_exclude(target, 'excluded by target group')
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


class WaveformPiggybackSubtarget(Object):
    piggy_id = Int.T()

    _next_piggy_id = 0

    @classmethod
    def new_piggy_id(cls):
        piggy_id = WaveformPiggybackSubtarget._next_piggy_id
        WaveformPiggybackSubtarget._next_piggy_id += 1
        return piggy_id

    def __init__(self, piggy_id=None, **kwargs):
        if piggy_id is None:
            piggy_id = self.new_piggy_id()

        Object.__init__(self, piggy_id=piggy_id, **kwargs)

    def evaluate(
            self, tr_proc_obs, trspec_proc_obs, tr_proc_syn, trspec_proc_syn):

        raise NotImplementedError()


class WaveformPiggybackSubresult(Object):
    piggy_id = Int.T()


class WaveformMisfitResult(gf.Result, MisfitResult):
    '''Carries the observations for a target and corresponding synthetics.

    A number of different waveform  or phase representations are possible.
    '''
    processed_obs = trace.Trace.T(optional=True)
    processed_syn = trace.Trace.T(optional=True)
    filtered_obs = trace.Trace.T(optional=True)
    filtered_syn = trace.Trace.T(optional=True)
    spectrum_obs = TraceSpectrum.T(optional=True)
    spectrum_syn = TraceSpectrum.T(optional=True)

    taper = trace.Taper.T(optional=True)
    tobs_shift = Float.T(optional=True)
    tsyn_pick = Timestamp.T(optional=True)
    tshift = Float.T(optional=True)
    cc = trace.Trace.T(optional=True)

    piggyback_subresults = List.T(WaveformPiggybackSubresult.T())


@has_get_plot_classes
class WaveformMisfitTarget(gf.Target, MisfitTarget):
    flip_norm = Bool.T(default=False)
    misfit_config = WaveformMisfitConfig.T()

    can_bootstrap_weights = True

    def __init__(self, **kwargs):
        gf.Target.__init__(self, **kwargs)
        MisfitTarget.__init__(self, **kwargs)
        self._piggyback_subtargets = []

    def string_id(self):
        return '.'.join(x for x in (self.path,) + self.codes)

    @classmethod
    def get_plot_classes(cls):
        from . import plot
        plots = super(WaveformMisfitTarget, cls).get_plot_classes()
        plots.extend(plot.get_plot_classes())
        return plots

    def get_combined_weight(self):
        if self._combined_weight is None:
            w = self.manual_weight
            for analyser in self.analyser_results.values():
                w *= analyser.weight
            self._combined_weight = num.array([w], dtype=num.float)
        return self._combined_weight

    def get_taper_params(self, engine, source):
        store = engine.get_store(self.store_id)
        config = self.misfit_config
        tmin_fit = source.time + store.t(config.tmin, source, self)
        tmax_fit = source.time + store.t(config.tmax, source, self)
        if config.fmin > 0.0:
            tfade = 1.0/config.fmin
        else:
            tfade = 1.0/config.fmax

        if config.tfade is None:
            tfade_taper = tfade
        else:
            tfade_taper = config.tfade

        return tmin_fit, tmax_fit, tfade, tfade_taper

    def get_backazimuth_for_waveform(self):
        return backazimuth_for_waveform(self.azimuth, self.codes)

    @property
    def backazimuth(self):
        return self.azimuth - 180.

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

        if self.misfit_config.fmin > 0:
            tinc_obs = 1.0 / self.misfit_config.fmin
        else:
            tinc_obs = 10.0 / self.misfit_config.fmax

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

        if config.quantity == 'displacement':
            syn_resp = None
        elif config.quantity == 'velocity':
            syn_resp = trace.DifferentiationResponse(1)
        elif config.quantity == 'acceleration':
            syn_resp = trace.DifferentiationResponse(2)
        else:
            GrondError('Unsupported quantity: %s' % config.quantity)

        tr_syn = tr_syn.transfer(
            freqlimits=freqlimits,
            tfade=tfade,
            transfer_function=syn_resp)

        tr_syn.chop(tmin_fit - 2*tfade, tmax_fit + 2*tfade)

        tmin_obs, tmax_obs = self.get_cutout_timespan(
            tmin_fit+tobs_shift, tmax_fit+tobs_shift, tfade)

        try:
            tr_obs = ds.get_waveform(
                nslc,
                quantity=config.quantity,
                tinc_cache=1.0/(config.fmin or 0.1*config.fmax),
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
                autoshift_penalty_max=config.autoshift_penalty_max,
                subtargets=self._piggyback_subtargets)

            self._piggyback_subtargets = []

            mr.tobs_shift = float(tobs_shift)
            mr.tsyn_pick = float_or_none(tsyn)

            return mr

        except NotFound as e:
            logger.debug(str(e))
            raise gf.SeismosizerError('No waveform data: %s' % str(e))

    def get_plain_targets(self, engine, source):
        d = dict(
            (k, getattr(self, k)) for k in gf.Target.T.propnames)
        return [gf.Target(**d)]

    def add_piggyback_subtarget(self, subtarget):
        self._piggyback_subtargets.append(subtarget)


def misfit(
        tr_obs, tr_syn, taper, domain, exponent, tautoshift_max,
        autoshift_penalty_max, flip, result_mode='sparse', subtargets=[]):

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

    piggyback_results = []
    for subtarget in subtargets:
        piggyback_results.append(
            subtarget.evaluate(
                tr_proc_obs, trspec_proc_obs, tr_proc_syn, trspec_proc_syn))

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

    elif domain == 'log_frequency_domain':
        a, b = trspec_proc_syn.ydata, trspec_proc_obs.ydata
        if flip:
            b, a = a, b

        a = num.abs(a)
        b = num.abs(b)

        eps = (num.mean(a) + num.mean(b)) * 1e-7
        if eps == 0.0:
            eps = 1e-7

        a = num.log(a + eps)
        b = num.log(b + eps)

        m, n = trace.Lx_norm(a, b, norm=exponent)

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

    result.piggyback_subresults = piggyback_results

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

    elif domain in ('frequency_domain', 'log_frequency_domain'):
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
    StoreIDSelectorError
    StoreIDSelector
    Crust2StoreIDSelector
    StationDictStoreIDSelector
    DepthRangeToStoreID
    StationDepthStoreIDSelector
    WaveformTargetGroup
    WaveformMisfitConfig
    WaveformMisfitTarget
    WaveformMisfitResult
    WaveformPiggybackSubtarget
    WaveformPiggybackSubresult
'''.split()
