
from pyrocko.client import catalog

import logging
import numpy as num
from pyrocko.guts import Int, Bool, Float, String, StringChoice
from pyrocko.gf.meta import OutOfBounds
from ..base import Analyser, AnalyserConfig, AnalyserResult
from grond.dataset import NotFound

logger = logging.getLogger('grond.analysers.NoiseAnalyser')


guts_prefix = 'grond'


def get_phase_arrival_time(engine, source, target, wavename):
    """
    Get arrival time from Green's Function store for respective
    :class:`pyrocko.gf.seismosizer.Target`,
    :class:`pyrocko.gf.meta.Location` pair.

    Parameters
    ----------
    engine : :class:`pyrocko.gf.seismosizer.LocalEngine`
    source : :class:`pyrocko.gf.meta.Location`
        can be therefore :class:`pyrocko.gf.seismosizer.Source` or
        :class:`pyrocko.model.Event`
    target : :class:`pyrocko.gf.seismosizer.Target`
    wavename : string
        of the tabulated phase_def that determines the phase arrival

    Returns
    -------
    scalar, float of the arrival time of the wave
    """
    store = engine.get_store(target.store_id)
    dist = target.distance_to(source)
    depth = source.depth
    return store.t(wavename, (depth, dist)) + source.time


def seismic_noise_variance(traces, engine, source, targets,
                           nwindows, pre_event_noise_duration,
                           check_events, phase_def):
    """
    Calculate variance of noise in a given time before P-Phase onset.

    Optionally check the gCMT earthquake catalogue for M>5 events interfering.

    Parameters
    ----------
    data_traces : list
        of :class:`pyrocko.trace.Trace` containing observed data
    engine : :class:`pyrocko.gf.seismosizer.LocalEngine`
        processing object for synthetics calculation
    source : :class:`pyrocko.gf.Source`
        reference source
    targets : list
        of :class:`pyrocko.gf.seismosizer.Targets`
    nwindows : integer
        number of windows in which the noise trace is split. If not 1, the
        variance is calculated for each window separately and a mean
        variance is returned. Else, the variance is calculated on the
        entire pre-event noise window.
    pre_event_noise_duration : Time before the first arrival to include in the
        noise analysis
    phase_def : :class:'pyrocko.gf.Timing'
    arrivals : list
        of :class'pyrocko.gf.Timing' arrivals of waveforms
        at station

    Returns
    -------
    :class:`numpy.ndarray`
    """

    var_ds = []
    global_cmt_catalog = catalog.GlobalCMT()
    var_ds = []
    ev_ws = []
    for tr, target in zip(traces, targets):
        stat_w = 1.

        if tr is None:
            var_ds.append(num.nan)
            ev_ws.append(num.nan)
        else:

            arrival_time = get_phase_arrival_time(
                engine=engine, source=source,
                target=target, wavename=phase_def)
            if check_events:
                events = global_cmt_catalog.get_events(
                    time_range=(
                        arrival_time-pre_event_noise_duration-50.*60.,
                        arrival_time),
                    magmin=5.,)
                for ev in events:
                    try:
                        arrival_time_pre = get_phase_arrival_time(
                            engine=engine,
                            source=ev,
                            target=target,
                            wavename=phase_def)

                        if arrival_time_pre > arrival_time \
                                - pre_event_noise_duration \
                                and arrival_time_pre < arrival_time:

                            stat_w = 0.
                            logger.info(
                                'Noise analyser found event "%s" phase onset '
                                'of "%s" for target "%s".' % (
                                    ev.name, phase_def, target.name))

                        if arrival_time_pre > arrival_time-30.*60.\
                                and arrival_time_pre < arrival_time - \
                                pre_event_noise_duration:
                            stat_w *= 1.
                            logger.info(
                                'Noise analyser found event "%s" possibly '
                                'contaminating the noise.' % ev.name)

                            # this should be magnitude dependent
                    except Exception:
                        pass
            ev_ws.append(stat_w)

            if nwindows == 1:
                vtrace_var = num.nanvar(tr.ydata)
                var_ds.append(vtrace_var)
            else:
                win = arrival_time - (arrival_time -
                                      pre_event_noise_duration)
                win_len = win/nwindows
                v_traces_w = []
                for i in range(0, nwindows):
                    vtrace_w = tr.chop(
                        tmin=win+win_len*i,
                        tmax=arrival_time+win_len*i+1,
                        inplace=False)
                    v_traces_w.append(vtrace_w.ydata)
                v_traces_w = num.nanmean(v_traces_w)
                var_ds.append(v_traces_w)
    var_ds = num.array(var_ds, dtype=num.float)
    ev_ws = num.array(ev_ws, dtype=num.float)
    return var_ds, ev_ws


class NoiseAnalyser(Analyser):
    '''
    From the pre-event station noise variance-based trace weights are formed.

    By default, the trace weights are the inverse of the noise variance. The
    correlation of the noise is neglected. Optionally, using a the gCMT global
    earthquake catalogue, the station data are checked for theoretical phase
    arrivals of M>5 earthquakes. In case of a very probable contamination the
    trace weights are set to zero. In case global earthquake phase arrivals are
    within a 30 min time window before the start of the set pre-event noise
    window, only a warning is thrown.

    It is further possible to disregard data with a noise level exceeding the
    median by a given ``cutoff`` factor. These weights are set to 0. This can
    be done exclusively (``mode='weeding'``) such that noise weights are either
    1 or 0, or in combination with weighting below the median-times-cutoff
    noise level (``mode='weighting'``).
    '''

    def __init__(self, nwindows, pre_event_noise_duration,
                 check_events, phase_def, statistic, mode, cutoff,
                 cutoff_exception_on_high_snr):

        Analyser.__init__(self)
        self.nwindows = nwindows
        self.pre_event_noise_duration = pre_event_noise_duration
        self.check_events = check_events
        self.phase_def = phase_def
        self.statistic = statistic
        self.mode = mode
        self.cutoff = cutoff
        self.cutoff_exception_on_high_snr = cutoff_exception_on_high_snr

    def analyse(self, problem, ds):

        tdur = self.pre_event_noise_duration

        if tdur == 0:
            return

        if not problem.has_waveforms:
            return

        engine = problem.get_engine()
        source = problem.base_source

        paths = sorted(set(t.path for t in problem.waveform_targets))

        for path in paths:
            targets = [t for t in problem.waveform_targets if t.path == path]

            deltats = set()
            for target in targets:  # deltat diff check?
                store = engine.get_store(target.store_id)
                deltats.add(float(store.config.deltat))

            if len(deltats) > 1:
                logger.warn(
                    'Differing sampling rates in stores used. Using highest.')

            deltat = min(deltats)

            data = []
            for target in targets:
                try:
                    freqlimits = list(target.get_freqlimits())
                    freqlimits = tuple(freqlimits)

                    source = problem.base_source

                    arrival_time = get_phase_arrival_time(
                        engine=engine,
                        source=source,
                        target=target,
                        wavename=self.phase_def)

                    tmin_fit, tmax_fit, tfade, tfade_taper = \
                        target.get_taper_params(engine, source)

                    data.append([
                        tmin_fit,
                        tmax_fit,
                        tfade_taper,
                        ds.get_waveform(
                            target.codes,
                            tmin=arrival_time-tdur-tfade,
                            tmax=tmax_fit+tfade,
                            tfade=tfade,
                            freqlimits=freqlimits,
                            deltat=deltat,
                            backazimuth=target.get_backazimuth_for_waveform(),
                            tinc_cache=1./freqlimits[0],
                            debug=False)])

                except (NotFound, OutOfBounds) as e:
                    logger.debug(str(e))
                    data.append([None, None, None, None])

            traces_noise = []
            traces_signal = []
            for tmin_fit, tmax_fit, tfade_taper, tr in data:
                if tr:
                    traces_noise.append(
                        tr.chop(tr.tmin, tr.tmin + tdur, inplace=False))
                    traces_signal.append(
                        tr.chop(
                            tmin_fit-tfade_taper,
                            tmax_fit+tfade_taper,
                            inplace=False))
                else:
                    traces_noise.append(None)
                    traces_signal.append(None)

            var_ds, ev_ws = seismic_noise_variance(
                traces_noise, engine, source, targets,
                self.nwindows, tdur,
                self.check_events, self.phase_def)

            amp_maxs = num.array([
                (tr.absmax()[1] if tr else num.nan) for tr in traces_signal])

            if self.statistic == 'var':
                noise = var_ds
            elif self.statistic == 'std':
                noise = num.sqrt(var_ds)
            else:
                assert False, 'invalid statistic argument'

            ok = num.isfinite(noise)

            if num.sum(ok) == 0:
                norm_noise = 0.0
            else:
                norm_noise = num.median(noise[ok])

            if norm_noise == 0.0:
                logger.info(
                    'Noise Analyser returned a weight of 0 for all stations.')

            assert num.all(noise[ok] >= 0.0)

            ce_factor = self.cutoff_exception_on_high_snr
            high_snr = num.zeros(ok.size, dtype=num.bool)
            if ce_factor is not None:
                high_snr[ok] = amp_maxs[ok] > ce_factor * num.sqrt(var_ds)[ok]

            weights = num.zeros(noise.size)
            if self.mode == 'weighting':
                weights[ok] = norm_noise / noise[ok]
            elif self.mode == 'weeding':
                weights[ok] = 1.0
            else:
                assert False, 'invalid mode argument'

            if self.cutoff is not None:
                weights[ok] = num.where(
                    num.logical_or(
                        noise[ok] <= norm_noise * self.cutoff,
                        high_snr[ok]),
                    weights[ok], 0.0)

            if self.check_events:
                weights = weights*ev_ws

            for itarget, target in enumerate(targets):
                logger.info((
                    'Noise analysis for target "%s":\n'
                    '  var: %g\n'
                    '  std: %g\n'
                    '  max/std: %g\n'
                    '  %s/median(%s): %g\n'
                    '  contamination_weight: %g\n'
                    '  weight: %g') % (
                        target.string_id(),
                        var_ds[itarget],
                        num.sqrt(var_ds[itarget]),
                        amp_maxs[itarget] / num.sqrt(var_ds[itarget]),
                        self.statistic, self.statistic,
                        noise[itarget] / norm_noise,
                        ev_ws[itarget],
                        weights[itarget]))

            for weight, target in zip(weights, targets):
                target.analyser_results['noise'] = \
                    NoiseAnalyserResult(weight=float(weight))


class NoiseAnalyserResult(AnalyserResult):
    weight = Float.T(
        help='The inverse of the pre-event data variance or standard '
             'deviation. If traces were checked for other event phase '
             'arrivals, the weight can be zero for contaminated traces.')


class NoiseAnalyserConfig(AnalyserConfig):
    """Configuration parameters for the pre-event noise analysis."""

    nwindows = Int.T(
        default=1,
        help='number of windows for trace splitting')

    pre_event_noise_duration = Float.T(
        default=0.,
        help='Total length of noise trace in the analysis')

    phase_def = String.T(
        default='P',
        help='Onset of phase_def used for upper limit of window')

    check_events = Bool.T(
        default=False,
        help='check the GlobalCMT for M>5 earthquakes'
             ' that produce phase arrivals'
             ' contaminating and affecting the noise analysis')

    statistic = StringChoice.T(
        choices=('var', 'std'),
        default='var',
        help='Set weight to inverse of noise variance (var) or standard '
             'deviation (std).')

    mode = StringChoice.T(
        choices=('weighting', 'weeding'),
        default='weighting',
        help='Generate weights based on inverse of noise measure (weighting), '
             'or discrete on/off style in combination with cutoff value '
             '(weeding).')

    cutoff = Float.T(
        optional=True,
        help='Set weight to zero, when noise level exceeds median by the '
             'given cutoff factor.')

    cutoff_exception_on_high_snr = Float.T(
        optional=True,
        help='Exclude from cutoff when max amplitude exceeds standard '
             'deviation times this factor.')

    def get_analyser(self):
        return NoiseAnalyser(
            nwindows=self.nwindows,
            pre_event_noise_duration=self.pre_event_noise_duration,
            check_events=self.check_events, phase_def=self.phase_def,
            statistic=self.statistic, mode=self.mode, cutoff=self.cutoff,
            cutoff_exception_on_high_snr=self.cutoff_exception_on_high_snr)


__all__ = '''
    NoiseAnalyser
    NoiseAnalyserConfig
'''.split()
