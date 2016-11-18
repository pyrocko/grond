
import copy
import logging
from collections import defaultdict
import numpy as num
from pyrocko import util, pile, model, config, trace, snuffling, gui_util
from pyrocko.fdsn import enhanced_sacpz, station as fs
from pyrocko.guts import Object, Tuple, String, Float, dump_all, load_all

guts_prefix = 'grond'

logger = logging.getLogger('grond.dataset')


class InvalidObject(Exception):
    pass


class NotFound(Exception):
    def __init__(self, reason, codes=None):
        self.reason = reason
        self.codes = codes

    def __str__(self):
        s = self.reason
        if self.codes:
            s += ' (%s)' % '.'.join(self.codes)

        return s


class DatasetError(Exception):
    pass


class StationCorrection(Object):
    codes = Tuple.T(4, String.T())
    delay = Float.T()
    factor = Float.T()


def load_station_corrections(filename):
    scs = load_all(filename=filename)
    for sc in scs:
        assert isinstance(sc, StationCorrection)

    return scs


def dump_station_corrections(station_corrections, filename):
    return dump_all(station_corrections, filename=filename)


class Dataset(object):

    def __init__(self, event_name=None):
        self.events = []
        self.pile = pile.Pile()
        self.stations = {}
        self.responses = defaultdict(list)
        self.responses_stationxml = []
        self.clippings = {}
        self.blacklist = set()
        self.whitelist_nslc = None
        self.whitelist_nsl = None
        self.station_corrections = {}
        self.station_factors = {}
        self.pick_markers = []
        self.apply_correction_delays = True
        self.apply_correction_factors = True
        self.clip_handling = 'by_nsl'
        self.synthetic_test = None
        self._picks = None
        self._cache = {}
        self._event_name = event_name

    def empty_cache(self):
        self._cache = {}

    def set_synthetic_test(self, synthetic_test):
        self.synthetic_test = synthetic_test

    def add_stations(
            self,
            stations=None,
            pyrocko_stations_filename=None,
            stationxml_filenames=None):

        if stations is not None:
            for station in stations:
                self.stations[station.nsl()] = station

        if pyrocko_stations_filename is not None:
            logger.debug(
                'Loading stations from file %s' %
                pyrocko_stations_filename)

            for station in model.load_stations(pyrocko_stations_filename):
                self.stations[station.nsl()] = station

        if stationxml_filenames is not None and len(stationxml_filenames) > 0:

            for stationxml_filename in stationxml_filenames:
                logger.debug(
                    'Loading stations from StationXML file %s' %
                    stationxml_filename)

                sx = fs.load_xml(filename=stationxml_filename)
                for station in sx.get_pyrocko_stations():
                    self.stations[station.nsl()] = station

    def add_events(self, events=None, filename=None):
        if events is not None:
            self.events.extend(events)

        if filename is not None:
            logger.debug('Loading events from file %s' % filename)
            self.events.extend(model.load_events(filename))

    def add_waveforms(self, paths, regex=None, fileformat='detect',
                      show_progress=True):
        logger.debug('Loading waveform data from %s' % paths)
        cachedirname = config.config().cache_dir
        fns = util.select_files(paths, regex=regex,
                                show_progress=show_progress)
        cache = pile.get_cache(cachedirname)
        self.pile.load_files(sorted(fns), cache=cache,
                             fileformat=fileformat,
                             show_progress=show_progress)

    def add_responses(self, sacpz_dirname=None, stationxml_filenames=None):
        if sacpz_dirname:
            logger.debug('Loading SAC PZ responses from %s' % sacpz_dirname)
            for x in enhanced_sacpz.iload_dirname(sacpz_dirname):
                self.responses[x.codes].append(x)

        if stationxml_filenames:
            for stationxml_filename in stationxml_filenames:
                logger.debug(
                    'Loading StationXML responses from %s' %
                    stationxml_filename)

                self.responses_stationxml.append(
                    fs.load_xml(filename=stationxml_filename))

    def add_clippings(self, markers_filename):
        markers = snuffling.load_markers(markers_filename)
        clippings = {}
        for marker in markers:
            nslc = marker.one_nslc()
            nsl = nslc[:3]
            if nsl not in clippings:
                clippings[nsl] = []

            if nslc not in clippings:
                clippings[nslc] = []

            clippings[nsl].append(marker.tmin)
            clippings[nslc].append(marker.tmin)

        for k, times in clippings.iteritems():
            atimes = num.array(times, dtype=num.float)
            if k not in self.clippings:
                self.clippings[k] = atimes
            else:
                self.clippings[k] = num.concatenate(self.clippings, atimes)

    def add_blacklist(self, blacklist=[], filenames=None):
        logger.debug('Loading blacklisted stations')
        if filenames:
            blacklist = list(blacklist)
            for filename in filenames:
                with open(filename, 'r') as f:
                    blacklist.extend(s.strip() for s in f.splitlines())

        for x in blacklist:
            if isinstance(x, basestring):
                x = tuple(x.split('.'))
            self.blacklist.add(x)

    def add_whitelist(self, whitelist=[], filenames=None):
        logger.debug('Loading whitelisted stations')
        if filenames:
            whitelist = list(whitelist)
            for filename in filenames:
                with open(filename, 'r') as f:
                    whitelist.extend(s.strip() for s in f.splitlines())

        if self.whitelist_nslc is None:
            self.whitelist_nslc = set()
            self.whitelist_nsl = set()
            self.whitelist_nsl_xx = set()

        for x in whitelist:
            if isinstance(x, basestring):
                x = tuple(x.split('.'))
            assert len(x) in (3, 4)
            if len(x) == 4:
                self.whitelist_nslc.add(x)
                self.whitelist_nsl_xx.add(x[:3])
            if len(x) == 3:
                self.whitelist_nsl.add(x)

    def add_station_corrections(self, filename):
        self.station_corrections.update(
            (sc.codes, sc) for sc in load_station_corrections(filename))

    def add_picks(self, filename):
        self.pick_markers.extend(
            gui_util.load_markers(filename))

        self._picks = None

    def is_blacklisted(self, obj):
        try:
            nslc = self.get_nslc(obj)
            if nslc in self.blacklist:
                return True

        except InvalidObject:
            pass

        nsl = self.get_nsl(obj)
        return (
            nsl in self.blacklist or
            nsl[1:2] in self.blacklist or
            nsl[:2] in self.blacklist)

    def is_whitelisted(self, obj):
        if self.whitelist_nslc is None:
            return True

        nsl = self.get_nsl(obj)
        try:
            nslc = self.get_nslc(obj)
            if nslc in self.whitelist_nslc:
                return True

            return nsl in self.whitelist_nsl

        except InvalidObject:
            return nsl in self.whitelist_nsl_xx or nsl in self.whitelist_nsl

    def has_clipping(self, nsl_or_nslc, tmin, tmax):
        if nsl_or_nslc not in self.clippings:
            return False

        atimes = self.clippings[nsl_or_nslc]
        return num.any(num.logical_and(tmin < atimes, atimes <= tmax))

    def get_nsl(self, obj):
        if isinstance(obj, trace.Trace):
            net, sta, loc, _ = obj.nslc_id
        elif isinstance(obj, model.Station):
            net, sta, loc = obj.nsl()
        elif isinstance(obj, tuple) and len(obj) in (3, 4):
            net, sta, loc = obj[:3]
        else:
            raise InvalidObject(
                'cannot get nsl code from given object of type %s' % type(obj))

        return net, sta, loc

    def get_nslc(self, obj):
        if isinstance(obj, trace.Trace):
            return obj.nslc_id
        elif isinstance(obj, tuple) and len(obj) == 4:
            return obj
        else:
            raise InvalidObject(
                'cannot get nslc code from given object %s' % type(obj))

    def get_tmin_tmax(self, obj):
        if isinstance(obj, trace.Trace):
            return obj.tmin, obj.tmax
        else:
            raise InvalidObject(
                'cannot get tmin and tmax from given object of type %s' %
                type(obj))

    def get_station(self, obj):
        if self.is_blacklisted(obj):
            raise NotFound('station is blacklisted', self.get_nsl(obj))

        if not self.is_whitelisted(obj):
            raise NotFound('station is not on whitelist', self.get_nsl(obj))

        if isinstance(obj, model.Station):
            return obj

        net, sta, loc = self.get_nsl(obj)

        keys = [(net, sta, loc), (net, sta, ''), ('', sta, '')]
        for k in keys:
            if k in self.stations:
                return self.stations[k]

        raise NotFound('no station information', keys)

    def get_stations(self):
        return [self.stations[k] for k in sorted(self.stations)
                if not self.is_blacklisted(self.stations[k])
                and self.is_whitelisted(self.stations[k])]

    def get_response(self, obj):
        if (self.responses is None or len(self.responses) == 0) \
                and (self.responses_stationxml is None
                     or len(self.responses_stationxml) == 0):

            raise NotFound('no response information available')

        if self.is_blacklisted(obj):
            raise NotFound('response is blacklisted', self.get_nslc(obj))

        if not self.is_whitelisted(obj):
            raise NotFound('response is not on whitelist', self.get_nslc(obj))

        net, sta, loc, cha = self.get_nslc(obj)
        tmin, tmax = self.get_tmin_tmax(obj)

        keys_x = [
            (net, sta, loc, cha), (net, sta, '', cha), ('', sta, '', cha)]

        keys = []
        for k in keys_x:
            if k not in keys:
                keys.append(k)

        candidates = []
        for k in keys:
            if k in self.responses:
                for x in self.responses[k]:
                    if x.tmin < tmin and (x.tmax is None or tmax < x.tmax):
                        candidates.append(x.response)

        for sx in self.responses_stationxml:
            try:
                candidates.append(
                    sx.get_pyrocko_response(
                        (net, sta, loc, cha),
                        timespan=(tmin, tmax),
                        fake_input_units='M'))

            except fs.NoResponseInformation, fs.MultipleResponseInformation:
                pass

        if len(candidates) == 1:
            return candidates[0]

        elif len(candidates) == 0:
            raise NotFound('no response found', (net, sta, loc, cha))
        else:
            raise NotFound('multiple responses found', (net, sta, loc, cha))

    def get_waveforms_raw(self, obj, tmin=None, tmax=None, tpad=0.):
        net, sta, loc = self.get_nsl(obj)

        trs = self.pile.all(
            tmin=tmin, tmax=tmax, tpad=tpad,
            trace_selector=lambda tr: tr.nslc_id[:3] == (net, sta, loc),
            want_incomplete=False)

        return trs

    def get_waveform_raw(
            self, obj,
            tmin=None,
            tmax=None,
            tpad=0.,
            toffset_noise_extract=0.):

        net, sta, loc, cha = self.get_nslc(obj)

        if self.is_blacklisted((net, sta, loc, cha)):
            raise NotFound(
                'waveform is blacklisted', (net, sta, loc, cha))

        if not self.is_whitelisted((net, sta, loc, cha)):
            raise NotFound(
                'waveform is not on whitelist', (net, sta, loc, cha))

        if self.clip_handling == 'by_nsl':
            if self.has_clipping((net, sta, loc), tmin, tmax):
                raise NotFound(
                    'waveform clipped', (net, sta, loc))

        elif self.clip_handling == 'by_nslc':
            if self.has_clipping((net, sta, loc, cha), tmin, tmax):
                raise NotFound(
                    'waveform clipped', (net, sta, loc, cha))

        trs = self.pile.all(
            tmin=tmin+toffset_noise_extract,
            tmax=tmax+toffset_noise_extract,
            tpad=tpad,
            trace_selector=lambda tr: tr.nslc_id == (net, sta, loc, cha),
            want_incomplete=False)

        if toffset_noise_extract != 0.0:
            for tr in trs:
                tr.shift(-toffset_noise_extract)

        if len(trs) == 1:
            return trs[0]

        else:
            raise NotFound(
                'waveform missing or incomplete', (net, sta, loc, cha))

    def get_waveform_restituted(
            self,
            obj, quantity='displacement',
            tmin=None, tmax=None, tpad=0.,
            tfade=0., freqlimits=None, deltat=None,
            toffset_noise_extract=0.):

        assert quantity == 'displacement'  # others not yet implemented

        tr = self.get_waveform_raw(
            obj, tmin=tmin, tmax=tmax, tpad=tpad+tfade,
            toffset_noise_extract=toffset_noise_extract)

        if deltat is not None:
            tr.downsample_to(deltat, snap=True, allow_upsample_max=5)
            tr.deltat = deltat

        resp = self.get_response(tr)
        return tr.transfer(tfade=tfade, freqlimits=freqlimits,
                           transfer_function=resp, invert=True)

    def _get_projections(
            self, station, backazimuth, source, target, tmin, tmax):

        # fill in missing channel information (happens when station file
        # does not contain any channel information)
        if not station.get_channels():
            station = copy.deepcopy(station)

            nsl = station.nsl()
            trs = self.pile.all(
                tmin=tmin, tmax=tmax,
                trace_selector=lambda tr: tr.nslc_id[:3] == nsl,
                load_data=False)

            channels = list(set(tr.channel for tr in trs))
            station.set_channels_by_name(*channels)

        projections = []
        projections.extend(station.guess_projections_to_enu(
            out_channels=('E', 'N', 'Z')))

        if source is not None and target is not None:
            backazimuth = source.azibazi_to(target)[1]

        if backazimuth is not None:
            projections.extend(station.guess_projections_to_rtu(
                out_channels=('R', 'T', 'Z'),
                backazimuth=backazimuth))

        if not projections:
            raise NotFound(
                'cannot determine projection of data components',
                station.nsl())

        return projections

    def get_waveform(
            self,
            obj, quantity='displacement',
            tmin=None, tmax=None, tpad=0.,
            tfade=0., freqlimits=None, deltat=None, cache=None,
            backazimuth=None,
            source=None,
            target=None):

        assert quantity == 'displacement'  # others not yet implemented

        if cache is True:
            cache = self._cache

        _, _, _, channel = self.get_nslc(obj)
        station = self.get_station(self.get_nsl(obj))

        nslc = station.nsl() + (channel,)

        if self.is_blacklisted(nslc):
            raise NotFound(
                'waveform is blacklisted', nslc)

        if not self.is_whitelisted(nslc):
            raise NotFound(
                'waveform is not on whitelist', nslc)

        if tmin is not None:
            tmin = float(tmin)

        if tmax is not None:
            tmax = float(tmax)

        if cache is not None and (nslc, tmin, tmax) in cache:
            obj = cache[nslc, tmin, tmax]
            if isinstance(obj, Exception):
                raise obj
            else:
                return obj

        syn_test = self.synthetic_test
        toffset_noise_extract = 0.0
        if syn_test:
            if not syn_test.respect_data_availability:
                if syn_test.add_real_noise:
                    raise DatasetError(
                        'respect_data_availability=False and '
                        'add_real_noise=True cannot be combined.')

                tr = syn_test.get_waveform(
                    nslc, tmin, tmax,
                    tfade=tfade,
                    freqlimits=freqlimits)

                if cache is not None:
                    cache[tr.nslc_id, tmin, tmax] = tr

                return tr

            if syn_test.add_real_noise:
                toffset_noise_extract = syn_test.toffset_real_noise

        abs_delays = []
        for ocha in 'ENZRT':
            sc = self.station_corrections.get(station.nsl() + (channel,), None)
            if sc:
                abs_delays.append(abs(sc.delay))

        if abs_delays:
            abs_delay_max = max(abs_delays)
        else:
            abs_delay_max = 0.0

        projections = self._get_projections(
            station, backazimuth, source, target, tmin, tmax)

        try:
            trs_projected = []
            for matrix, in_channels, out_channels in projections:
                deps = trace.project_dependencies(
                    matrix, in_channels, out_channels)

                trs = []
                if channel in deps:
                    for cha in deps[channel]:
                        trs.append(self.get_waveform_restituted(
                            station.nsl() + (cha,),
                            tmin=tmin, tmax=tmax, tpad=tpad+abs_delay_max,
                            toffset_noise_extract=toffset_noise_extract,
                            tfade=tfade, freqlimits=freqlimits, deltat=deltat))

                    trs_projected.extend(
                        trace.project(trs, matrix, in_channels, out_channels))

            for tr in trs_projected:
                sc = self.station_corrections.get(tr.nslc_id, None)
                if sc:
                    if self.apply_correction_factors:
                        tr.ydata /= sc.factor

                    if self.apply_correction_delays:
                        tr.shift(-sc.delay)

                if tmin is not None and tmax is not None:
                    tr.chop(tmin, tmax)

            if syn_test:
                trs_projected_synthetic = []
                for tr in trs_projected:
                    tr_syn = syn_test.get_waveform(
                        tr.nslc_id, tmin, tmax,
                        tfade=tfade, freqlimits=freqlimits)

                    if tr_syn:
                        if syn_test.add_real_noise:
                            tr_syn = tr_syn.copy()
                            tr_syn.add(tr)

                        trs_projected_synthetic.append(tr_syn)

                trs_projected = trs_projected_synthetic

            if cache is not None:
                for tr in trs_projected:
                    cache[tr.nslc_id, tmin, tmax] = tr

            for tr in trs_projected:
                if tr.channel == channel:
                    return tr

            raise NotFound(
                'waveform not available', station.nsl() + (channel,))

        except NotFound, e:
            cache[nslc, tmin, tmax] = e
            raise

    def get_events(self, magmin=None, event_names=None):
        evs = []
        for ev in self.events:
            if ((magmin is None or ev.magnitude >= magmin) and
                    (event_names is None or ev.name in event_names)):
                evs.append(ev)

        return evs

    def get_event_by_time(self, t, magmin=None):
        evs = self.get_events(magmin=magmin)
        ev_x = None
        for ev in evs:
            if ev_x is None or abs(ev.time - t) < abs(ev_x.time - t):
                ev_x = ev

        if not ev_x:
            raise NotFound(
                'no event information matching criteria (t=%s, magmin=%s)' %
                (t, magmin))

        return ev_x

    def get_event(self):
        if self._event_name is None:
            raise NotFound('no main event selected in dataset')

        for ev in self.events:
            if ev.name == self._event_name:
                return ev

        raise NotFound('no such event: %s' % self._event_name)

    def get_picks(self):
        if self._picks is None:
            hash_to_name = {}
            names = set()
            for marker in self.pick_markers:
                if isinstance(marker, gui_util.EventMarker):
                    name = marker.get_event().name
                    if name in names:
                        raise DatasetError(
                            'duplicate event name "%s" in picks' % name)

                    names.add(name)
                    hash_to_name[marker.get_event_hash()] = name

            picks = {}
            for marker in self.pick_markers:
                if isinstance(marker, gui_util.PhaseMarker):
                    ehash = marker.get_event_hash()

                    nsl = marker.one_nslc()[:3]
                    phasename = marker.get_phasename()

                    if ehash is None or ehash not in hash_to_name:
                        raise DatasetError(
                            'unassociated pick %s.%s.%s, %s' %
                            (nsl + (phasename, )))

                    eventname = hash_to_name[ehash]

                    if (nsl, phasename, eventname) in picks:
                        raise DatasetError(
                            'duplicate pick %s.%s.%s, %s' %
                            (nsl + (phasename, )))

                    picks[nsl, phasename, eventname] = marker

            self._picks = picks

        return self._picks

    def get_pick(self, eventname, obj, phasename):
        nsl = self.get_nsl(obj)
        return self.get_picks().get((nsl, phasename, eventname), None)


__all__ = '''
    DatasetError
    InvalidObject
    NotFound
    StationCorrection
    Dataset
    load_station_corrections
    dump_station_corrections
'''.split()
