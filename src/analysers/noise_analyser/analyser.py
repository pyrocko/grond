
from pyrocko.client import catalog

import logging
import numpy as num
from pyrocko.guts import Int, Bool, Float, String
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


def seismic_noise_variance(traces, engine, event, targets,
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
    event : :class:`pyrocko.meta.Event`
        reference event from catalog
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
                engine=engine, source=event,
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
                                'Noise analyser found event %s phase onset of '
                                '%s for %s' % (
                                    ev.name, phase_def, target.name))

                        if arrival_time_pre > arrival_time-30.*60.\
                                and arrival_time_pre < arrival_time - \
                                pre_event_noise_duration:
                            stat_w *= 1.
                            logger.info(
                                'Noise analyser found event %s possibly '
                                'contaminating the noise' % ev.name)

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
    """From the pre-event station noise variance-based trace weights are formed.

    The trace weights are the inverse of the noise variance. The correlation
    of the noise is neglected.
    Optionally, using a the gCMT global earthquake catalogue, the station data
    are checked for theoretical phase arrivals of M>5 earthquakes. In case of
    a very probable contamination the trace weights are set to zero. In case
    global earthquake phase arrivals are within a 30 min time window before
    the start of the set pre-event noise window, only a warning is thrown.
    """
    def __init__(self, nwindows, pre_event_noise_duration,
                 check_events, phase_def):
        Analyser.__init__(self)
        self.nwindows = nwindows
        self.pre_event_noise_duration = pre_event_noise_duration
        self.check_events = check_events
        self.phase_def = phase_def

    def analyse(self, problem, ds):

        if self.pre_event_noise_duration == 0:
            return

        if not problem.has_waveforms:
            return

        traces = []
        engine = problem.get_engine()
        event = ds.get_event()
        for target in problem.waveform_targets:  # deltat diff check?
                store = engine.get_store(target.store_id)
                deltat = store.config.deltat
        for target in problem.waveform_targets:
                try:
                    freqlimits = list(target.get_freqlimits())
                    freqlimits = tuple(freqlimits)
                    tfade = 1./freqlimits[0]
                    arrival_time = get_phase_arrival_time(
                        engine=engine, source=event,
                        target=target, wavename=self.phase_def)
                    traces.append(
                        ds.get_waveform(
                            target.codes,
                            tmin=arrival_time-self.pre_event_noise_duration,
                            tmax=arrival_time-1./freqlimits[0],
                            tfade=tfade,
                            freqlimits=freqlimits,
                            deltat=deltat,
                            backazimuth=target.
                            get_backazimuth_for_waveform(),
                            tinc_cache=1./freqlimits[0],
                            debug=False))

                except (NotFound, OutOfBounds) as e:
                    logger.debug(str(e))
                    traces.append(None)

        var_ds, ev_ws = seismic_noise_variance(
            traces, engine, event, problem.waveform_targets,
            self.nwindows, self.pre_event_noise_duration,
            self.check_events, self.phase_def)
        norm_noise = num.nanmedian(var_ds)
        if norm_noise == 0:
                logger.info(
                    'Noise Analyser returned a weight of 0 for all stations')

        weights0 = norm_noise/var_ds
        weights = num.where(var_ds==num.nan, 0., weights0)
        if self.check_events:
            weights = weights*ev_ws

        for weight, target in zip(weights, problem.waveform_targets):
            target.analyser_results['noise'] = \
                NoiseAnalyserResult(weight=float(weight))


class NoiseAnalyserResult(AnalyserResult):
    weight = Float.T(
        help='The inverse of the pre-event data variance. '
             'If traces were checked for other event phase arrivals,'
             ' the weight can be zero for contaminated traces.')


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
        help='check the gCMT global catalogue for M>5 earthquakes'
             ' that produce phase arrivals'
             ' contaminating and affecting the noise analysis')

    def get_analyser(self):
        return NoiseAnalyser(
            nwindows=self.nwindows,
            pre_event_noise_duration=self.pre_event_noise_duration,
            check_events=self.check_events, phase_def=self.phase_def)


__all__ = '''
    NoiseAnalyser
    NoiseAnalyserConfig
'''.split()
