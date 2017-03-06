import logging
import math
import numpy as num

from pyrocko import gf, trace, weeding
from pyrocko.guts import (Object, String, Float, Bool, Int, List, Dict,
                          StringChoice, Timestamp)
from pyrocko.guts_array import Array

from .dataset import NotFound

guts_prefix = 'grond'
logger = logging.getLogger('grond.target')


def float_or_none(x):
    if x is None:
        return x
    else:
        return float(x)


class DomainChoice(StringChoice):
    choices = [
        'time_domain',
        'frequency_domain',
        'envelope',
        'absolute',
        'cc_max_norm']


class TargetAnalysisResult(Object):
    balancing_weight = Float.T()


class NoAnalysisResults(Exception):
    pass


class Trace(Object):
    pass


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


class InnerMisfitConfig(Object):
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

    def get_full_frequency_range(self):
        return self.fmin / self.ffactor, self.fmax * self.ffactor


class MisfitResult(gf.Result):
    misfit_value = Float.T()
    misfit_norm = Float.T()
    processed_obs = Trace.T(optional=True)
    processed_syn = Trace.T(optional=True)
    filtered_obs = Trace.T(optional=True)
    filtered_syn = Trace.T(optional=True)
    spectrum_obs = TraceSpectrum.T(optional=True)
    spectrum_syn = TraceSpectrum.T(optional=True)
    taper = trace.Taper.T(optional=True)
    tobs_shift = Float.T(optional=True)
    tsyn_pick = Timestamp.T(optional=True)
    cc_shift = Float.T(optional=True)
    cc = Trace.T(optional=True)


class MisfitTarget(gf.Target):
    misfit_config = InnerMisfitConfig.T()
    flip_norm = Bool.T(default=False)
    manual_weight = Float.T(default=1.0)
    analysis_result = TargetAnalysisResult.T(optional=True)
    super_group = gf.StringID.T()
    group = gf.StringID.T()

    def __init__(self, **kwargs):
        gf.Target.__init__(self, **kwargs)
        self._ds = None
        self._result_mode = 'sparse'

    def string_id(self):
        return '.'.join(x for x in (
            self.super_group, self.group) + self.codes if x)

    @property
    def id(self):
        return self.codes

    def get_plain_target(self):
        d = dict(
            (k, getattr(self, k)) for k in gf.Target.T.propnames)
        return gf.Target(**d)

    def get_dataset(self):
        return self._ds

    def set_dataset(self, ds):
        self._ds = ds

    def set_result_mode(self, result_mode):
        self._result_mode = result_mode

    def get_combined_weight(self, apply_balancing_weights):
        w = self.manual_weight
        if apply_balancing_weights:
            w *= self.get_balancing_weight()

        return w

    def get_balancing_weight(self):
        if not self.analysis_result:
            raise NoAnalysisResults('no balancing weights available')

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
        nslc = self.codes
        if nslc[-1] == 'R':
            backazimuth = self.azimuth + 180.
        elif nslc[-1] == 'T':
            backazimuth = self.azimuth + 90.
        else:
            backazimuth = None

        return backazimuth

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
                tmin=tmin_obs,
                tmax=tmax_obs,
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
                exponent=2,
                flip=self.flip_norm,
                result_mode=self._result_mode,
                tautoshift_max=config.tautoshift_max,
                autoshift_penalty_max=config.autoshift_penalty_max)

            mr.tobs_shift = float(tobs_shift)
            mr.tsyn_pick = float_or_none(tsyn)

            return mr

        except NotFound, e:
            logger.debug(str(e))
            raise gf.SeismosizerError('no waveform data, %s' % str(e))


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

    :returns: object of type :py:class:`MisfitResult`
    '''

    trace.assert_same_sampling_rate(tr_obs, tr_syn)
    tmin, tmax = taper.time_span()

    tr_proc_obs, trspec_proc_obs = _process(tr_obs, tmin, tmax, taper, domain)
    tr_proc_syn, trspec_proc_syn = _process(tr_syn, tmin, tmax, taper, domain)

    cc_shift = None
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
            for ishift in xrange(-nshift_max, nshift_max+1):
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

        cc_shift, cc_max = ctr.max()
        m = 0.5 - 0.5 * cc_max
        n = 0.5

    elif domain == 'frequency_domain':
        a, b = trspec_proc_syn.ydata, trspec_proc_obs.ydata
        if flip:
            b, a = a, b

        m, n = trace.Lx_norm(num.abs(a), num.abs(b), norm=exponent)

    if result_mode == 'full':
        result = MisfitResult(
            misfit_value=m,
            misfit_norm=n,
            processed_obs=tr_proc_obs,
            processed_syn=tr_proc_syn,
            filtered_obs=tr_obs.copy(),
            filtered_syn=tr_syn,
            spectrum_obs=trspec_proc_obs,
            spectrum_syn=trspec_proc_syn,
            taper=taper,
            cc_shift=cc_shift,
            cc=ctr)

    elif result_mode == 'sparse':
        result = MisfitResult(
            misfit_value=m,
            misfit_norm=n)
    else:
        assert False

    return result


def _extend_extract(tr, tmin, tmax):
    deltat = tr.deltat
    itmin_frame = int(math.floor(tmin/deltat))
    itmax_frame = int(math.ceil(tmax/deltat))
    nframe = itmax_frame - itmin_frame
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
    tr.tmin = tmin
    tr.set_ydata(a)
    return tr


def _process(tr, tmin, tmax, taper, domain):
    tr_proc = _extend_extract(tr, tmin, tmax)
    tr_proc.taper(taper)

    df = None
    trspec_proc = None

    if domain == 'envelope':
        tr_proc = tr_proc.envelope(inplace=False)

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


class InnerSatelliteMisfitConfig(Object):
    use_weight_focal = Bool.T(default=False)
    ranges = Dict.T(String.T(), gf.Range.T(),
                             default={'waterlevel': '-0.5 .. 0.5',
                                      'ramp_north': '-1e-4 .. 1e-4',
                                      'ramp_east': '-1e-4 .. 1e-4'})



class MisfitSatelliteTarget(gf.SatelliteTarget):
    scene_id = String.T()
    super_group = gf.StringID.T()
    inner_misfit_config = InnerSatelliteMisfitConfig.T()
    manual_weight = Float.T(default=1.0)
    group = gf.StringID.T()

    target_parameters = [
        ...]

    ranges = Dict.T(String.T(), gf.Range.T())

    def __init__(self, *args, **kwargs):
        gf.SatelliteTarget.__init__(self, *args, **kwargs)
        self.waterlevel = 0.
        self.ramp_north = 0.
        self.ramp_east = 0.
        self._ds = None
        self._distance_cache = None

    def set_dataset(self, ds):
        self._ds = ds

    @property
    def id(self):
        return self.scene_id

    def get_dataset(self):
        return self._ds

    def set_result_mode(self, result_mode):
        pass

    def string_id(self):
        return '.'.join([self.super_group, self.group, self.scene_id])

    def set_scene_levels(self, waterlevel, ramp_north, ramp_east):
        self.waterlevel = waterlevel
        self.ramp_north = ramp_north
        self.ramp_east = ramp_east

    def post_process(self, engine, source, statics):
        scene = self._ds.get_kite_scene(self.scene_id)
        quadtree = scene.quadtree

        stat_obs = scene.quadtree.leaf_medians
        stat_syn = statics['displacement.los']

        stat_level = num.zeros_like(stat_obs)
        stat_level.fill(self.waterlevel)
        stat_level += (quadtree.leaf_center_distance[:, 0] * self.ramp_east)
        stat_level += (quadtree.leaf_center_distance[:, 1] * self.ramp_north)

        res = num.abs(stat_obs - (stat_syn + stat_level))

        misfit_value = num.sqrt(
            num.sum((res * scene.covariance.weight_vector)**2))

        result = MisfitResult(
            misfit_value=misfit_value,
            misfit_norm=2)
        return result

    def get_combined_weight(self, apply_balancing_weights=False):
        return 1.


class GrondTarget(MisfitSatelliteTarget, MisfitTarget):
    def __init__(self, args):
        super(MisfitSatelliteTarget, self).__init__()


class TargetConfig(Object):

    super_group = gf.StringID.T(default='', optional=True)
    group = gf.StringID.T(optional=True)
    distance_min = Float.T(optional=True)
    distance_max = Float.T(optional=True)
    distance_3d_min = Float.T(optional=True)
    distance_3d_max = Float.T(optional=True)
    depth_min = Float.T(optional=True)
    depth_max = Float.T(optional=True)
    limit = Int.T(optional=True)
    channels = List.T(String.T(), optional=True)
    inner_misfit_config = InnerMisfitConfig.T(optional=True)
    inner_satellite_misfit_config = InnerSatelliteMisfitConfig.T(
        optional=True)
    interpolation = gf.InterpolationMethod.T()
    store_id = gf.StringID.T()
    weight = Float.T(default=1.0)

    def get_targets(self, ds, event, default_group):

        origin = event

        targets = []

        for scene in ds.get_kite_scenes():
            qt = scene.quadtree

            lats = num.empty(qt.nleafs)
            lons = num.empty(qt.nleafs)
            lats.fill(qt.frame.llLat)
            lons.fill(qt.frame.llLon)

            east_shifts = qt.leaf_focal_points[:, 0]
            north_shifts = qt.leaf_focal_points[:, 1]

            sat_target = MisfitSatelliteTarget(
                quantity='displacement',
                scene_id=scene.meta.scene_id,
                lats=lats,
                lons=lons,
                east_shifts=east_shifts,
                north_shifts=north_shifts,
                theta=qt.leaf_thetas,
                phi=qt.leaf_phis,
                tsnapshot=None,
                interpolation=self.interpolation,
                store_id=self.store_id,
                super_group=self.super_group,
                group=self.group or default_group,
                inner_misfit_config=self.inner_satellite_misfit_config)

            sat_target.set_dataset(ds)
            targets.append(sat_target)

        for st in ds.get_stations():
            for cha in self.channels:
                if ds.is_blacklisted((st.nsl() + (cha,))):
                    continue

                target = MisfitTarget(
                    quantity='displacement',
                    codes=st.nsl() + (cha,),
                    lat=st.lat,
                    lon=st.lon,
                    depth=st.depth,
                    interpolation=self.interpolation,
                    store_id=self.store_id,
                    misfit_config=self.inner_misfit_config,
                    manual_weight=self.weight,
                    super_group=self.super_group,
                    group=self.group or default_group)

                if self.distance_min is not None and \
                        target.distance_to(origin) < self.distance_min:
                    continue

                if self.distance_max is not None and \
                        target.distance_to(origin) > self.distance_max:
                    continue

                if self.distance_3d_min is not None and \
                        target.distance_3d_to(origin) < self.distance_3d_min:
                    continue

                if self.distance_3d_max is not None and \
                        target.distance_3d_to(origin) > self.distance_3d_max:
                    continue

                if self.depth_min is not None and \
                        target.depth < self.depth_min:
                    continue

                if self.depth_max is not None and \
                        target.depth > self.depth_max:
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
