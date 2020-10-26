import glob
import copy
import os.path as op
import logging
import math
import numpy as num

from collections import defaultdict
from pyrocko import util, pile, model, config, trace, \
    marker as pmarker
from pyrocko.io.io_common import FileLoadError
from pyrocko.fdsn import enhanced_sacpz, station as fs
from pyrocko.guts import (Object, Tuple, String, Float, List, Bool, dump_all,
                          load_all)

from pyrocko import gf

from .meta import Path, HasPaths, expand_template, GrondError

from .synthetic_tests import SyntheticTest

guts_prefix = 'grond'
logger = logging.getLogger('grond.dataset')


def quote_paths(paths):
    return ', '.join('"%s"' % path for path in paths)


class InvalidObject(Exception):
    pass


class NotFound(Exception):
    def __init__(self, reason, codes=None, time_range=None):
        self.reason = reason
        self.time_range = time_range
        self.codes = codes

    def __str__(self):
        s = self.reason
        if self.codes:
            s += ' (%s)' % '.'.join(self.codes)

        if self.time_range:
            s += ' (%s - %s)' % (
                util.time_to_str(self.time_range[0]),
                util.time_to_str(self.time_range[1]))

        return s


class DatasetError(GrondError):
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
        self.whitelist_nsl_xx = None
        self.whitelist = None
        self.station_corrections = {}
        self.station_factors = {}
        self.pick_markers = []
        self.apply_correction_delays = True
        self.apply_correction_factors = True
        self.apply_displaced_sampling_workaround = False
        self.extend_incomplete = False
        self.clip_handling = 'by_nsl'
        self.kite_scenes = []
        self.gnss_campaigns = []
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
                'Loading stations from file "%s"...' %
                pyrocko_stations_filename)

            for station in model.load_stations(pyrocko_stations_filename):
                self.stations[station.nsl()] = station

        if stationxml_filenames is not None and len(stationxml_filenames) > 0:

            for stationxml_filename in stationxml_filenames:
                if not op.exists(stationxml_filename):
                    continue

                logger.debug(
                    'Loading stations from StationXML file "%s"...' %
                    stationxml_filename)

                sx = fs.load_xml(filename=stationxml_filename)
                ev = self.get_event()
                stations = sx.get_pyrocko_stations(time=ev.time)
                if len(stations) == 0:
                    logger.warning(
                        'No stations found for time %s in file "%s".' % (
                            util.time_to_str(ev.time), stationxml_filename))

                for station in stations:
                    logger.debug('Adding station: %s.%s.%s' % station.nsl())
                    channels = station.get_channels()
                    if len(channels) == 1 and channels[0].name.endswith('Z'):
                        logger.warning(
                            'Station "%s" has vertical component'
                            ' information only, adding mocked channels.'
                            % station.nsl_string())
                        station.add_channel(model.Channel('N'))
                        station.add_channel(model.Channel('E'))

                    self.stations[station.nsl()] = station

    def add_events(self, events=None, filename=None):
        if events is not None:
            self.events.extend(events)

        if filename is not None:
            logger.debug('Loading events from file "%s"...' % filename)
            try:
                events = model.load_events(filename)
                self.events.extend(events)
                logger.info(
                    'Loading events from %s: %i events found.' %
                    (filename, len(events)))
            except Exception as e:
                logger.warning('Could not load events from %s!', filename)
                raise e

    def add_waveforms(self, paths, regex=None, fileformat='detect',
                      show_progress=False):
        cachedirname = config.config().cache_dir

        logger.debug('Selecting waveform files %s...' % quote_paths(paths))
        fns = util.select_files(paths, regex=regex,
                                show_progress=show_progress)
        cache = pile.get_cache(cachedirname)
        logger.debug('Scanning waveform files %s...' % quote_paths(paths))
        self.pile.load_files(sorted(fns), cache=cache,
                             fileformat=fileformat,
                             show_progress=show_progress)

    def add_responses(self, sacpz_dirname=None, stationxml_filenames=None):
        if sacpz_dirname:
            logger.debug(
                'Loading SAC PZ responses from "%s"...' % sacpz_dirname)
            for x in enhanced_sacpz.iload_dirname(sacpz_dirname):
                self.responses[x.codes].append(x)

        if stationxml_filenames:
            for stationxml_filename in stationxml_filenames:
                if not op.exists(stationxml_filename):
                    continue

                logger.debug(
                    'Loading StationXML responses from "%s"...' %
                    stationxml_filename)

                self.responses_stationxml.append(
                    fs.load_xml(filename=stationxml_filename))

    def add_clippings(self, markers_filename):
        markers = pmarker.load_markers(markers_filename)
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

        for k, times in clippings.items():
            atimes = num.array(times, dtype=num.float)
            if k not in self.clippings:
                self.clippings[k] = atimes
            else:
                self.clippings[k] = num.concatenate(self.clippings, atimes)

    def add_blacklist(self, blacklist=[], filenames=None):
        logger.debug('Loading blacklisted stations...')
        if filenames:
            blacklist = list(blacklist)
            for filename in filenames:
                if op.exists(filename):
                    with open(filename, 'r') as f:
                        blacklist.extend(
                            s.strip() for s in f.read().splitlines())
                else:
                    logger.warning('No such blacklist file: %s' % filename)

        for x in blacklist:
            if isinstance(x, str):
                x = tuple(x.split('.'))
            self.blacklist.add(x)

    def add_whitelist(self, whitelist=[], filenames=None):
        logger.debug('Loading whitelisted stations...')
        if filenames:
            whitelist = list(whitelist)
            for filename in filenames:
                with open(filename, 'r') as f:
                    whitelist.extend(s.strip() for s in f.read().splitlines())

        if self.whitelist_nslc is None:
            self.whitelist_nslc = set()
            self.whitelist = set()
            self.whitelist_nsl_xx = set()

        for x in whitelist:
            if isinstance(x, str):
                x = tuple(x.split('.'))
            if len(x) == 4:
                self.whitelist_nslc.add(x)
                self.whitelist_nsl_xx.add(x[:3])
            else:
                self.whitelist.add(x)

    def add_station_corrections(self, filename):
        self.station_corrections.update(
            (sc.codes, sc) for sc in load_station_corrections(filename))

    def add_picks(self, filename):
        self.pick_markers.extend(
            pmarker.load_markers(filename))

        self._picks = None

    def add_gnss_campaigns(self, paths):
        paths = util.select_files(
            paths,
            regex=r'\.yml|\.yaml',
            show_progress=False)

        for path in paths:
            self.add_gnss_campaign(filename=path)

    def add_gnss_campaign(self, filename):
        try:
            from pyrocko.model import gnss  # noqa
        except ImportError:
            raise ImportError('Module pyrocko.model.gnss not found,'
                              ' please upgrade pyrocko!')
        logger.debug('Loading GNSS campaign from "%s"...' % filename)

        campaign = load_all(filename=filename)
        self.gnss_campaigns.append(campaign[0])

    def add_kite_scenes(self, paths):
        logger.info('Loading kite InSAR scenes...')
        paths = util.select_files(
            paths,
            regex=r'\.npz',
            show_progress=False)

        for path in paths:
            self.add_kite_scene(filename=path)

        if not self.kite_scenes:
            logger.warning('Could not find any kite scenes at %s.' %
                           quote_paths(self.kite_scene_paths))

    def add_kite_scene(self, filename):
        try:
            from kite import Scene
        except ImportError:
            raise ImportError('Module kite could not be imported,'
                              ' please install from https://pyrocko.org.')
        logger.debug('Loading kite scene from "%s"...' % filename)

        scene = Scene.load(filename)
        scene._log.setLevel(logger.level)

        try:
            self.get_kite_scene(scene.meta.scene_id)
        except NotFound:
            self.kite_scenes.append(scene)
        else:
            raise AttributeError('Kite scene_id not unique for "%s".'
                                 % filename)

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

        if (
                nsl in self.whitelist or
                nsl[1:2] in self.whitelist or
                nsl[:2] in self.whitelist):

            return True

        try:
            nslc = self.get_nslc(obj)
            if nslc in self.whitelist_nslc:
                return True

        except InvalidObject:
            return nsl in self.whitelist_nsl_xx

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
        elif isinstance(obj, gf.Target):
            net, sta, loc, _ = obj.codes
        elif isinstance(obj, tuple) and len(obj) in (3, 4):
            net, sta, loc = obj[:3]
        else:
            raise InvalidObject(
                'Cannot get nsl code from given object of type "%s".'
                % type(obj))

        return net, sta, loc

    def get_nslc(self, obj):
        if isinstance(obj, trace.Trace):
            return obj.nslc_id
        elif isinstance(obj, gf.Target):
            return obj.codes
        elif isinstance(obj, tuple) and len(obj) == 4:
            return obj
        else:
            raise InvalidObject(
                'Cannot get nslc code from given object "%s"' % type(obj))

    def get_tmin_tmax(self, obj):
        if isinstance(obj, trace.Trace):
            return obj.tmin, obj.tmax
        else:
            raise InvalidObject(
                'Cannot get tmin and tmax from given object of type "%s"' %
                type(obj))

    def get_station(self, obj):
        if self.is_blacklisted(obj):
            raise NotFound('Station is blacklisted:', self.get_nsl(obj))

        if not self.is_whitelisted(obj):
            raise NotFound('Station is not on whitelist:', self.get_nsl(obj))

        if isinstance(obj, model.Station):
            return obj

        net, sta, loc = self.get_nsl(obj)

        keys = [(net, sta, loc), (net, sta, ''), ('', sta, '')]
        for k in keys:
            if k in self.stations:
                return self.stations[k]

        raise NotFound('No station information:', keys)

    def get_stations(self):
        return [self.stations[k] for k in sorted(self.stations)
                if not self.is_blacklisted(self.stations[k])
                and self.is_whitelisted(self.stations[k])]

    def get_kite_scenes(self):
        return self.kite_scenes

    def get_kite_scene(self, scene_id=None):
        if scene_id is None:
            if len(self.kite_scenes) == 0:
                raise AttributeError('No kite displacements defined.')
            return self.kite_scenes[0]
        else:
            for scene in self.kite_scenes:
                if scene.meta.scene_id == scene_id:
                    return scene
        raise NotFound('No kite scene with id "%s" defined.' % scene_id)

    def get_gnss_campaigns(self):
        return self.gnss_campaigns

    def get_gnss_campaign(self, name):
        for camp in self.gnss_campaigns:
            if camp.name == name:
                return camp
        raise NotFound('GNSS campaign %s not found!' % name)

    def get_response(self, obj, quantity='displacement'):
        if (self.responses is None or len(self.responses) == 0) \
                and (self.responses_stationxml is None
                     or len(self.responses_stationxml) == 0):

            raise NotFound('No response information available.')

        quantity_to_unit = {
            'displacement': 'M',
            'velocity': 'M/S',
            'acceleration': 'M/S**2'}

        if self.is_blacklisted(obj):
            raise NotFound('Response is blacklisted:', self.get_nslc(obj))

        if not self.is_whitelisted(obj):
            raise NotFound('Response is not on whitelist:', self.get_nslc(obj))

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
                        if quantity == 'displacement':
                            candidates.append(x.response)
                        elif quantity == 'velocity':
                            candidates.append(trace.MultiplyResponse([
                                x.response,
                                trace.DifferentiationResponse()]))
                        elif quantity == 'acceleration':
                            candidates.append(trace.MultiplyResponse([
                                x.response,
                                trace.DifferentiationResponse(2)]))
                        else:
                            assert False

        for sx in self.responses_stationxml:
            try:
                candidates.append(
                    sx.get_pyrocko_response(
                        (net, sta, loc, cha),
                        timespan=(tmin, tmax),
                        fake_input_units=quantity_to_unit[quantity]))

            except (fs.NoResponseInformation, fs.MultipleResponseInformation):
                pass

        if len(candidates) == 1:
            return candidates[0]

        elif len(candidates) == 0:
            raise NotFound('No response found:', (net, sta, loc, cha))
        else:
            raise NotFound('Multiple responses found:', (net, sta, loc, cha))

    def get_waveform_raw(
            self, obj,
            tmin,
            tmax,
            tpad=0.,
            toffset_noise_extract=0.,
            want_incomplete=False,
            extend_incomplete=False):

        net, sta, loc, cha = self.get_nslc(obj)

        if self.is_blacklisted((net, sta, loc, cha)):
            raise NotFound(
                'Waveform is blacklisted:', (net, sta, loc, cha))

        if not self.is_whitelisted((net, sta, loc, cha)):
            raise NotFound(
                'Waveform is not on whitelist:', (net, sta, loc, cha))

        if self.clip_handling == 'by_nsl':
            if self.has_clipping((net, sta, loc), tmin, tmax):
                raise NotFound(
                    'Waveform clipped:', (net, sta, loc))

        elif self.clip_handling == 'by_nslc':
            if self.has_clipping((net, sta, loc, cha), tmin, tmax):
                raise NotFound(
                    'Waveform clipped:', (net, sta, loc, cha))

        trs = self.pile.all(
            tmin=tmin+toffset_noise_extract,
            tmax=tmax+toffset_noise_extract,
            tpad=tpad,
            trace_selector=lambda tr: tr.nslc_id == (net, sta, loc, cha),
            want_incomplete=want_incomplete or extend_incomplete)

        if self.apply_displaced_sampling_workaround:
            for tr in trs:
                tr.snap()

        if toffset_noise_extract != 0.0:
            for tr in trs:
                tr.shift(-toffset_noise_extract)

        if extend_incomplete and len(trs) == 1:
            trs[0].extend(
                tmin + toffset_noise_extract - tpad,
                tmax + toffset_noise_extract + tpad,
                fillmethod='repeat')

        if not want_incomplete and len(trs) != 1:
            if len(trs) == 0:
                message = 'Waveform missing or incomplete.'
            else:
                message = 'Waveform has gaps.'

            raise NotFound(
                message,
                codes=(net, sta, loc, cha),
                time_range=(
                    tmin + toffset_noise_extract - tpad,
                    tmax + toffset_noise_extract + tpad))

        return trs

    def get_waveform_restituted(
            self,
            obj, quantity='displacement',
            tmin=None, tmax=None, tpad=0.,
            tfade=0., freqlimits=None, deltat=None,
            toffset_noise_extract=0.,
            want_incomplete=False,
            extend_incomplete=False):

        trs_raw = self.get_waveform_raw(
            obj, tmin=tmin, tmax=tmax, tpad=tpad+tfade,
            toffset_noise_extract=toffset_noise_extract,
            want_incomplete=want_incomplete,
            extend_incomplete=extend_incomplete)

        trs_restituted = []
        for tr in trs_raw:
            if deltat is not None:
                tr.downsample_to(deltat, snap=True, allow_upsample_max=5)
                tr.deltat = deltat

            resp = self.get_response(tr, quantity=quantity)
            try:
                trs_restituted.append(
                    tr.transfer(
                        tfade=tfade, freqlimits=freqlimits,
                        transfer_function=resp, invert=True, demean=True))

            except trace.InfiniteResponse:
                raise NotFound(
                    'Instrument response deconvolution failed '
                    '(divide by zero)', tr.nslc_id)

        return trs_restituted, trs_raw

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
                'Cannot determine projection of data components:',
                station.nsl())

        return projections

    def _get_waveform(
            self,
            obj, quantity='displacement',
            tmin=None, tmax=None, tpad=0.,
            tfade=0., freqlimits=None, deltat=None, cache=None,
            backazimuth=None,
            source=None,
            target=None,
            debug=False):

        assert not debug or (debug and cache is None)

        if cache is True:
            cache = self._cache

        _, _, _, channel = self.get_nslc(obj)
        station = self.get_station(self.get_nsl(obj))

        nslc = station.nsl() + (channel,)

        if self.is_blacklisted(nslc):
            raise NotFound(
                'Waveform is blacklisted:', nslc)

        if not self.is_whitelisted(nslc):
            raise NotFound(
                'Waveform is not on whitelist:', nslc)

        assert tmin is not None
        assert tmax is not None

        tmin = float(tmin)
        tmax = float(tmax)

        nslc = tuple(nslc)

        cache_k = nslc + (
            tmin, tmax, tuple(freqlimits), tfade, deltat, tpad, quantity)
        if cache is not None and (nslc + cache_k) in cache:
            obj = cache[nslc + cache_k]
            if isinstance(obj, Exception):
                raise obj
            elif obj is None:
                raise NotFound('Waveform not found!', nslc)
            else:
                return obj

        syn_test = self.synthetic_test
        toffset_noise_extract = 0.0
        if syn_test:
            if not syn_test.respect_data_availability:
                if syn_test.real_noise_scale != 0.0:
                    raise DatasetError(
                        'respect_data_availability=False and '
                        'addition of real noise cannot be combined.')

                tr = syn_test.get_waveform(
                    nslc, tmin, tmax,
                    tfade=tfade,
                    freqlimits=freqlimits)

                if cache is not None:
                    cache[tr.nslc_id + cache_k] = tr

                if debug:
                    return [], [], []
                else:
                    return tr

            if syn_test.real_noise_scale != 0.0:
                toffset_noise_extract = syn_test.real_noise_toffset

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
            trs_restituted = []
            trs_raw = []
            exceptions = []
            for matrix, in_channels, out_channels in projections:
                deps = trace.project_dependencies(
                    matrix, in_channels, out_channels)

                try:
                    trs_restituted_group = []
                    trs_raw_group = []
                    if channel in deps:
                        for cha in deps[channel]:
                            trs_restituted_this, trs_raw_this = \
                                self.get_waveform_restituted(
                                    station.nsl() + (cha,),
                                    quantity=quantity,
                                    tmin=tmin, tmax=tmax,
                                    tpad=tpad+abs_delay_max,
                                    toffset_noise_extract=toffset_noise_extract,  # noqa
                                    tfade=tfade,
                                    freqlimits=freqlimits,
                                    deltat=deltat,
                                    want_incomplete=debug,
                                    extend_incomplete=self.extend_incomplete)

                            trs_restituted_group.extend(trs_restituted_this)
                            trs_raw_group.extend(trs_raw_this)

                        trs_projected.extend(
                            trace.project(
                                trs_restituted_group, matrix,
                                in_channels, out_channels))

                        trs_restituted.extend(trs_restituted_group)
                        trs_raw.extend(trs_raw_group)

                except NotFound as e:
                    exceptions.append((in_channels, out_channels, e))

            if not trs_projected:
                err = []
                for (in_channels, out_channels, e) in exceptions:
                    sin = ', '.join(c.name for c in in_channels)
                    sout = ', '.join(c.name for c in out_channels)
                    err.append('(%s) -> (%s): %s' % (sin, sout, e))

                raise NotFound('\n'.join(err))

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
                    if tr.channel == channel:
                        tr_syn = syn_test.get_waveform(
                            tr.nslc_id, tmin, tmax,
                            tfade=tfade, freqlimits=freqlimits)

                        if tr_syn:
                            if syn_test.real_noise_scale != 0.0:
                                tr_syn = tr_syn.copy()
                                tr_noise = tr.copy()
                                tr_noise.set_ydata(
                                    tr_noise.get_ydata()
                                    * syn_test.real_noise_scale)

                                tr_syn.add(tr_noise)

                            trs_projected_synthetic.append(tr_syn)

                trs_projected = trs_projected_synthetic

            if cache is not None:
                for tr in trs_projected:
                    cache[tr.nslc_id + cache_k] = tr

            tr_return = None
            for tr in trs_projected:
                if tr.channel == channel:
                    tr_return = tr

            if debug:
                return trs_projected, trs_restituted, trs_raw, tr_return

            elif tr_return:
                return tr_return

            else:
                raise NotFound(
                    'waveform not available', station.nsl() + (channel,))

        except NotFound:
            if cache is not None:
                cache[nslc + cache_k] = None
            raise

    def get_waveform(self, obj, tinc_cache=None, **kwargs):
        tmin = kwargs['tmin']
        tmax = kwargs['tmax']
        if tinc_cache is not None:
            tmin_r = (math.floor(tmin / tinc_cache) - 1.0) * tinc_cache
            tmax_r = (math.ceil(tmax / tinc_cache) + 1.0) * tinc_cache
        else:
            tmin_r = tmin
            tmax_r = tmax

        kwargs['tmin'] = tmin_r
        kwargs['tmax'] = tmax_r

        if kwargs.get('debug', None):
            return self._get_waveform(obj, **kwargs)
        else:
            tr = self._get_waveform(obj, **kwargs)
            return tr.chop(tmin, tmax, inplace=False)

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
                'No event information matching criteria (t=%s, magmin=%s).' %
                (t, magmin))

        return ev_x

    def get_event(self):
        if self._event_name is None:
            raise NotFound('No main event selected in dataset!')

        for ev in self.events:
            if ev.name == self._event_name:
                return ev

        raise NotFound('No such event: %s' % self._event_name)

    def get_picks(self):
        if self._picks is None:
            hash_to_name = {}
            names = set()
            for marker in self.pick_markers:
                if isinstance(marker, pmarker.EventMarker):
                    name = marker.get_event().name
                    if name in names:
                        raise DatasetError(
                            'Duplicate event name "%s" in picks.' % name)

                    names.add(name)
                    hash_to_name[marker.get_event_hash()] = name

            for ev in self.events:
                hash_to_name[ev.get_hash()] = ev.name

            picks = {}
            for marker in self.pick_markers:
                if isinstance(marker, pmarker.PhaseMarker):
                    ehash = marker.get_event_hash()

                    nsl = marker.one_nslc()[:3]
                    phasename = marker.get_phasename()

                    if ehash is None or ehash not in hash_to_name:
                        raise DatasetError(
                            'Unassociated pick: %s.%s.%s, %s' %
                            (nsl + (phasename, )))

                    eventname = hash_to_name[ehash]

                    if (nsl, phasename, eventname) in picks:
                        raise DatasetError(
                            'Duplicate pick: %s.%s.%s, %s' %
                            (nsl + (phasename, )))

                    picks[nsl, phasename, eventname] = marker

            self._picks = picks

        return self._picks

    def get_pick(self, eventname, obj, phasename):
        nsl = self.get_nsl(obj)
        return self.get_picks().get((nsl, phasename, eventname), None)


class DatasetConfig(HasPaths):
    ''' Configuration for a Grond `Dataset`  object. '''

    stations_path = Path.T(
        optional=True,
        help='List of files with station coordinates in Pyrocko format.')
    stations_stationxml_paths = List.T(
        Path.T(),
        optional=True,
        help='List of files with station coordinates in StationXML format.')
    events_path = Path.T(
        optional=True,
        help='File with hypocenter information and possibly'
             ' reference solution')
    waveform_paths = List.T(
        Path.T(),
        optional=True,
        help='List of directories with raw waveform data')
    clippings_path = Path.T(
        optional=True)
    responses_sacpz_path = Path.T(
        optional=True,
        help='List of SACPZ response files for restitution of'
             ' the raw waveform data.')
    responses_stationxml_paths = List.T(
        Path.T(),
        optional=True,
        help='List of StationXML response files for restitution of'
             ' the raw waveform data.')
    station_corrections_path = Path.T(
        optional=True,
        help='File containing station correction informations.')
    apply_correction_factors = Bool.T(
        optional=True,
        default=True,
        help='Apply correction factors from station corrections.')
    apply_correction_delays = Bool.T(
        optional=True,
        default=True,
        help='Apply correction delays from station corrections.')
    apply_displaced_sampling_workaround = Bool.T(
        optional=True,
        default=False,
        help='Work around displaced sampling issues.')
    extend_incomplete = Bool.T(
        default=False,
        help='Extend incomplete seismic traces.')
    picks_paths = List.T(
        Path.T())
    blacklist_paths = List.T(
        Path.T(),
        help='List of text files with blacklisted stations.')
    blacklist = List.T(
        String.T(),
        help='Stations/components to be excluded according to their STA, '
             'NET.STA, NET.STA.LOC, or NET.STA.LOC.CHA codes.')
    whitelist_paths = List.T(
        Path.T(),
        help='List of text files with whitelisted stations.')
    whitelist = List.T(
        String.T(),
        optional=True,
        help='If not None, list of stations/components to include according '
             'to their STA, NET.STA, NET.STA.LOC, or NET.STA.LOC.CHA codes. '
             'Note: ''when whitelisting on channel level, both, the raw and '
             'the processed channel codes have to be listed.')
    synthetic_test = SyntheticTest.T(
        optional=True)

    kite_scene_paths = List.T(
        Path.T(),
        optional=True,
        help='List of directories for the InSAR scenes.')

    gnss_campaign_paths = List.T(
        Path.T(),
        optional=True,
        help='List of directories for the GNSS campaign data.')

    def __init__(self, *args, **kwargs):
        HasPaths.__init__(self, *args, **kwargs)
        self._ds = {}

    def get_event_names(self):
        logger.info('Loading events ...')

        def extra(path):
            return expand_template(path, dict(
                event_name='*'))

        def fp(path):
            return self.expand_path(path, extra=extra)

        def check_events(events, fn):
            for ev in events:
                if not ev.name:
                    logger.warning('Event in %s has no name!', fn)
                    return
                if not ev.lat or not ev.lon:
                    logger.warning('Event %s has inconsistent coordinates!',
                                   ev.name)
                if not ev.depth:
                    logger.warning('Event %s has no depth!', ev.name)
                if not ev.time:
                    logger.warning('Event %s has no time!', ev.name)

        events = []
        events_path = fp(self.events_path)
        fns = glob.glob(events_path)
        if not fns:
            raise DatasetError('No event files matching "%s".' % events_path)

        for fn in fns:
            logger.debug('Loading from file %s' % fn)
            ev = model.load_events(filename=fn)
            check_events(ev, fn)

            events.extend(ev)

        event_names = [ev_.name for ev_ in events]
        event_names.sort()
        return event_names

    def get_dataset(self, event_name):
        if event_name not in self._ds:
            def extra(path):
                return expand_template(path, dict(
                    event_name=event_name))

            def fp(path):
                p = self.expand_path(path, extra=extra)
                if p is None:
                    return None

                if isinstance(p, list):
                    for path in p:
                        if not op.exists(path):
                            logger.warn('Path %s does not exist.' % path)
                else:
                    if not op.exists(p):
                        logger.warn('Path %s does not exist.' % p)

                return p

            ds = Dataset(event_name)
            try:
                ds.add_events(filename=fp(self.events_path))

                ds.add_stations(
                    pyrocko_stations_filename=fp(self.stations_path),
                    stationxml_filenames=fp(self.stations_stationxml_paths))

                if self.waveform_paths:
                    ds.add_waveforms(paths=fp(self.waveform_paths))

                if self.kite_scene_paths:
                    ds.add_kite_scenes(paths=fp(self.kite_scene_paths))

                if self.gnss_campaign_paths:
                    ds.add_gnss_campaigns(paths=fp(self.gnss_campaign_paths))

                if self.clippings_path:
                    ds.add_clippings(markers_filename=fp(self.clippings_path))

                if self.responses_sacpz_path:
                    ds.add_responses(
                        sacpz_dirname=fp(self.responses_sacpz_path))

                if self.responses_stationxml_paths:
                    ds.add_responses(
                        stationxml_filenames=fp(
                            self.responses_stationxml_paths))

                if self.station_corrections_path:
                    ds.add_station_corrections(
                        filename=fp(self.station_corrections_path))

                ds.apply_correction_factors = self.apply_correction_factors
                ds.apply_correction_delays = self.apply_correction_delays
                ds.apply_displaced_sampling_workaround = \
                    self.apply_displaced_sampling_workaround
                ds.extend_incomplete = self.extend_incomplete

                for picks_path in self.picks_paths:
                    ds.add_picks(
                        filename=fp(picks_path))

                ds.add_blacklist(self.blacklist)
                ds.add_blacklist(filenames=fp(self.blacklist_paths))
                if self.whitelist:
                    ds.add_whitelist(self.whitelist)
                if self.whitelist_paths:
                    ds.add_whitelist(filenames=fp(self.whitelist_paths))

                ds.set_synthetic_test(copy.deepcopy(self.synthetic_test))
                self._ds[event_name] = ds
            except (FileLoadError, OSError) as e:
                raise DatasetError(str(e))

        return self._ds[event_name]


__all__ = '''
    Dataset
    DatasetConfig
    DatasetError
    InvalidObject
    NotFound
    StationCorrection
    load_station_corrections
    dump_station_corrections
'''.split()
