
from pyrocko.client import catalog

import logging
import numpy as num
from pyrocko.guts import Int, Bool, Float

from ..base import Analyser, AnalyserConfig, AnalyserResult

logger = logging.getLogger('grond.analysers.target_balancer')


guts_prefix = 'grond'


def get_phase_arrival_time(engine, source, target, wavename):
    """
    Get arrival time from Greens Function store for respective
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
        of the tabulated phase that determines the phase arrival

    Returns
    -------
    scalar, float of the arrival time of the wave
    """
    store = engine.get_store(target.store_id)
    dist = target.distance_to(source)
    depth = source.depth
    return store.t(wavename, (depth, dist)) + source.time


def seismic_noise_variance(data_traces, engine, event, targets,
                           nsamples, pre_ev_noise_len, check_evs):
    '''
    Calculate variance of noise (half an hour) before P-Phase onset, and check
    for other events interfering

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
    nsamples : integer
        if given, windowed analysis of noise, else
        variance is calculated on the entire pre-event
        noise

    Returns
    -------
    :class:`numpy.ndarray`
    '''

    wavename = 'P'   # hardcode here, want always pre P time
    var_ds = []
    global_cmt_catalog = catalog.GlobalCMT()
    var_ds = []
    ev_ws = []
    for tr, target in zip(data_traces, targets):
        stat_w = 1.
        arrival_time = get_phase_arrival_time(
            engine=engine, source=event,
            target=target, wavename=wavename)
        if check_evs:
            events = global_cmt_catalog.get_events(
                time_range=(
                    arrival_time-pre_ev_noise_len-40.*60.,
                    arrival_time-60.),
                magmin=6.,)
            for ev in events:
                try:
                    arrival_time_pre = get_phase_arrival_time(
                        engine=engine,
                        source=ev,
                        target=target,
                        wavename=wavename)

                    if arrival_time_pre > arrival_time-pre_ev_noise_len and \
                            arrival_time_pre < arrival_time:
                        stat_w = 0.

                    ev_ws.append(stat_w)

                except Exception:
                    pass

        if nsamples == 1:
            vtrace = tr.chop(
                tmin=arrival_time-pre_ev_noise_len,
                tmax=arrival_time-60.,
                inplace=False)

            vtrace_var = num.var(vtrace.ydata)
            var_ds.append(vtrace_var)
        else:
            win = arrival_time - 60. - tr.tmin
            win_len = win/nsamples
            v_traces_w = []
            for i in range(0, nsamples):
                vtrace_w = tr.chop(
                    tmin=tr.tmin+win_len*i,
                    tmax=tr.tmin+win_len*i+1,
                    inplace=False)
                v_traces_w.append(vtrace_w.ydata)
            v_traces_w = num.mean(v_traces_w)
            var_ds.append(v_traces_w)

    return var_ds, ev_ws


class NoiseAnalyser(Analyser):

    def __init__(self, nsamples, pre_ev_noise_len, check_evs):
        Analyser.__init__(self)
        self.nsamples = nsamples
        self.pre_ev_noise_len = pre_ev_noise_len
        self.check_evs = check_evs

    def analyse(self, problem, ds):

        if self.pre_ev_noise_len == 0:
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
                tmin_fit, tmax_fit, tfade, tfade_taper = \
                    target.get_taper_params(engine, problem.base_source)

                freqlimits = list(target.get_freqlimits())
                freqlimits = tuple(freqlimits)

                trs_projected, trs_restituted, trs_raw, tr_return = \
                    ds.get_waveform(
                        target.codes,
                        tmin=tmin_fit-self.pre_ev_noise_len,
                        tmax=tmin_fit,
                        tfade=tfade,
                        freqlimits=freqlimits,
                        deltat=deltat,
                        backazimuth=target.
                        get_backazimuth_for_waveform(),
                        tinc_cache=1./freqlimits[0],
                        debug=True)
                traces.append(tr_return)

        wproblem = problem.copy()

        var_ds, ev_ws = seismic_noise_variance(
            traces, engine, event, problem.waveform_targets,
            self.nsamples, self.pre_ev_noise_len, self.check_evs)

        norm_noise = num.median(var_ds)
        weights = norm_noise/var_ds
        if self.check_evs:
            weights = weights*ev_ws

        families, nfamilies = wproblem.get_family_mask()

        for ifamily in range(nfamilies):
            weights[families == ifamily] /= (
                num.nansum(weights[families == ifamily]) /
                num.nansum(num.isfinite(weights[families == ifamily])))

        for weight, target in zip(weights, problem.waveform_targets):
            target.analyser_results['noise'] = \
                NoiseAnalyserResult(weight=float(weight))


class NoiseAnalyserResult(AnalyserResult):
    weight = Float.T()


class NoiseAnalyserConfig(AnalyserConfig):
    nsamples = Int.T(default=1)
    pre_ev_noise_len = Float.T(default=0.)
    check_evs = Bool.T(default=False)

    def get_analyser(self):
        return NoiseAnalyser(
            nsamples=self.nsamples,
            pre_ev_noise_len=self.pre_ev_noise_len,
            check_evs=self.check_evs)


__all__ = '''
    NoiseAnalyser
    NoiseAnalyserConfig
'''.split()
