import numpy as num
from pyrocko import gf, trace
from pyrocko.guts import Object, Float, StringChoice, List, String
from pyrocko.gui import marker
from grond.meta import store_t


guts_prefix = 'grond'


response_wwssn_lp = trace.PoleZeroResponse(
    poles=[
        -0.2513 + 0.3351J,
        -0.2513 - 0.3351J,
        -0.0628 + 0.0304J,
        -0.0628 - 0.0304J],
    zeros=[0., 0., 0.],
    constant=383.6)


response_wwssn_lp_2 = trace.PoleZeroResponse(
    poles=[
        -0.40180 + 0.08559J,
        -0.40180 - 0.08559J,
        -0.04841,
        -0.08816],
    zeros=[0., 0., 0.])

response_wa = trace.PoleZeroResponse(
    poles=[
        -5.49779 - 5.60886J,
        -5.49779 + 5.60886J],
    zeros=[0., 0.])


class NamedResponse(StringChoice):
    choices = [
        'wood-anderson',
        'wwssn-lp']

    map = {
        'wood-anderson': response_wa,
        'wwssn-lp': response_wwssn_lp}


class BruneResponse(trace.FrequencyResponse):

    duration = Float.T()

    def evaluate(self, freqs):
        return 1.0 / (1.0 + (freqs*self.duration)**2)

# from pyrocko.plot import response as pr
# pr.plot([response_wwssn_lp, response_wwssn_lp_2, response_wa])


class FeatureMethod(StringChoice):
    choices = [
        'peak_component',
        'peak_to_peak_component',
        'peak_absolute_vector',
        'spectral_average']


class WaveformQuantity(StringChoice):
    choices = [
        'displacement',
        'velocity',
        'acceleration']


class FeatureMeasurementFailed(Exception):
    pass


class FeatureMeasure(Object):
    name = String.T()
    timing_tmin = gf.Timing.T(default='vel:8')
    timing_tmax = gf.Timing.T(default='vel:2')
    fmin = Float.T(optional=True)
    fmax = Float.T(optional=True)
    response = trace.FrequencyResponse.T(optional=True)
    named_response = NamedResponse.T(optional=True)
    channels = List.T(String.T())
    quantity = WaveformQuantity.T(default='displacement')
    method = FeatureMethod.T(default='peak_component')

    def get_nmodelling_targets(self):
        return len(self.channels)

    def get_modelling_targets(
            self, codes, lat, lon, depth, store_id, backazimuth):

        mtargets = []
        for channel in self.channels:
            target = gf.Target(
                quantity='displacement',
                codes=codes + (channel, ),
                lat=lat,
                lon=lon,
                depth=depth,
                store_id=store_id)

            if channel == 'R':
                target.azimuth = backazimuth - 180.
                target.dip = 0.
            elif channel == 'T':
                target.azimuth = backazimuth - 90.
                target.dip = 0.
            elif channel == 'Z':
                target.azimuth = 0.
                target.dip = -90.

            mtargets.append(target)

        return mtargets

    def evaluate(
            self, engine, source, targets,
            dataset=None,
            trs=None,
            extra_responses=[],
            debug=False):

        from ..waveform import target as base

        trs_processed = []
        trs_orig = []
        for itarget, target in enumerate(targets):
            if target.codes[-1] not in self.channels:
                continue

            store = engine.get_store(target.store_id)

            tmin = source.time + store_t(
                store, self.timing_tmin, source, target)
            tmax = source.time + store_t(
                store, self.timing_tmax, source, target)

            if self.fmin is not None and self.fmax is not None:
                freqlimits = [
                    self.fmin/2.,
                    self.fmin,
                    self.fmax,
                    self.fmax*2.]
                tfade = 1./self.fmin

            else:
                freqlimits = None
                tfade = 0.0

            if dataset is not None:
                bazi = base.backazimuth_for_waveform(
                    target.azimuth, target.codes)

                tr = dataset.get_waveform(
                    target.codes,
                    tinc_cache=1.0/self.fmin,
                    quantity=self.quantity,
                    tmin=tmin,
                    tmax=tmax,
                    freqlimits=freqlimits,
                    tfade=tfade,
                    deltat=store.config.deltat,
                    cache=True,
                    backazimuth=bazi)

            else:
                tr = trs[itarget]

                tr.extend(
                    tmin - tfade,
                    tmax + tfade,
                    fillmethod='repeat')

                tr = tr.transfer(
                    freqlimits=freqlimits,
                    tfade=tfade)

            trs_orig.append(tr)

            tr = tr.copy()

            responses = []
            responses.extend(extra_responses)

            ndiff = \
                WaveformQuantity.choices.index(self.quantity) - \
                WaveformQuantity.choices.index(target.quantity)

            if ndiff > 0:
                responses.append(trace.DifferentiationResponse(ndiff))

            if ndiff < 0:
                responses.append(trace.IntegrationResponse(-ndiff))

            if self.response:
                responses.append(self.response)

            if self.named_response:
                responses.append(
                    NamedResponse.map[self.named_response])

            if responses:
                trans = trace.MultiplyResponse(responses)
                try:
                    tr = tr.transfer(transfer_function=trans)

                except trace.TraceTooShort:
                    raise FeatureMeasurementFailed(
                        'transfer: trace too short')

            if tmin is None or tmax is None:
                raise FeatureMeasurementFailed(
                    'timing determination failed (phase unavailable?)')

            tr.chop(tmin, tmax)

            tr.set_location(tr.location + '-' + self.name + '-proc')
            trs_processed.append(tr)

        markers = []
        marker_candidates = []
        if self.method in ['peak_component', 'peak_to_peak_component']:
            component_amp_maxs = []
            for tr in trs_processed:
                y = tr.get_ydata()
                if self.method == 'peak_component':
                    yabs = num.abs(y)
                    i_at_amax = num.argmax(yabs)
                    amax = yabs[i_at_amax]
                    if debug:
                        t_at_amax = tr.tmin + i_at_amax * tr.deltat
                        mark = marker.Marker(
                            [tr.nslc_id],
                            t_at_amax,
                            t_at_amax,
                            0)

                        marker_candidates.append(mark)

                    component_amp_maxs.append(amax)
                else:
                    i_at_amax = num.argmax(y)
                    i_at_amin = num.argmin(y)
                    amax = y[i_at_amax]
                    amin = y[i_at_amin]
                    if debug:
                        t_at_amax = tr.tmin + i_at_amax * tr.deltat
                        t_at_amin = tr.tmin + i_at_amin * tr.deltat
                        ts = sorted([t_at_amax, t_at_amin])
                        mark = marker.Marker(
                            [tr.nslc_id],
                            ts[0],
                            ts[1],
                            0)

                        marker_candidates.append(mark)

                    component_amp_maxs.append(amax - amin)

            i_at_amax = num.argmax(component_amp_maxs)
            if debug:
                markers.append(marker_candidates[i_at_amax])
            amp_max = component_amp_maxs[i_at_amax]

        elif self.method == 'peak_absolute_vector':
            trsum = None
            for tr in trs_processed:
                tr.set_ydata(tr.get_ydata()**2)
                if trsum is None:
                    trsum = tr
                else:
                    trsum.add(tr)

            trsum.set_ydata(num.sqrt(tr.get_ydata))
            trsum.set_codes(channel='SUM')

            yabs = trsum.get_ydata()

            i_at_amax = num.argmax(yabs)
            amax = yabs[i_at_amax]
            t_at_amax = tr.tmin + i_at_amax * tr.deltat
            amp_max = amax

            if debug:
                markers.append(marker.Marker(
                    [trsum.nslc_id],
                    t_at_amax,
                    t_at_amax,
                    0))

            trs_processed.append(trsum)

        elif self.method == 'spectral_average':
            component_amp_maxs = []
            for tr in trs_processed:
                freqs, values = tr.spectrum()
                component_amp_maxs.append(num.mean(num.abs(values[
                    num.logical_and(self.fmin <= freqs, freqs <= self.fmax)])))

            amp_max = num.mean(component_amp_maxs)

        if debug:
            trs_out = []
            for tr in trs_orig:
                tr_out = tr.copy()
                tr_out.set_location(tr_out.location + '-' + self.name)
                trs_out.append(tr_out)

            return amp_max, (trs_out + trs_processed, markers)

        return amp_max, None
