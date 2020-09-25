from __future__ import print_function

import sys
import logging
import time
import copy
import shutil
import glob
import math
import os
import numpy as num
from contextlib import contextmanager

from pyrocko.guts import Object, String, Float, List
from pyrocko import gf, trace, guts, util, weeding
from pyrocko import parimap, model, marker as pmarker

from .dataset import NotFound, InvalidObject
from .problems.base import Problem, load_problem_info_and_data, \
    load_problem_data, ProblemDataNotAvailable

from .optimisers.base import BadProblem
from .targets.waveform.target import WaveformMisfitResult
from .targets.base import dump_misfit_result_collection, \
    MisfitResultCollection, MisfitResult, MisfitResultError
from .meta import expand_template, GrondError, selected
from .environment import Environment
from .monitor import GrondMonitor

logger = logging.getLogger('grond.core')
guts_prefix = 'grond'
op = os.path


class RingBuffer(num.ndarray):
    def __new__(cls, *args, **kwargs):
        cls = num.ndarray.__new__(cls, *args, **kwargs)
        cls.fill(0.)
        return cls

    def __init__(self, *args, **kwargs):
        self.pos = 0

    def put(self, value):
        self[self.pos] = value
        self.pos += 1
        self.pos %= self.size


def mahalanobis_distance(xs, mx, cov):
    imask = num.diag(cov) != 0.
    icov = num.linalg.inv(cov[imask, :][:, imask])
    temp = xs[:, imask] - mx[imask]
    return num.sqrt(num.sum(temp * num.dot(icov, temp.T).T, axis=1))


@contextmanager
def lock_rundir(rundir):
    statefn = op.join(rundir, '.running')
    if op.exists(statefn):
        raise EnvironmentError('file %s already exists!' % statefn)
    try:
        with open(statefn, 'w') as f:
            f.write('')
        yield True
    finally:
        os.remove(statefn)


class DirectoryAlreadyExists(GrondError):
    pass


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


def sarr(a):
    return ' '.join('%15g' % x for x in a)


def forward(env, show='filtered'):
    payload = []
    if env.have_rundir():
        env.setup_modelling()
        history = env.get_history(subset='harvest')
        xbest = history.get_best_model()
        problem = env.get_problem()
        ds = env.get_dataset()
        payload.append((ds, problem, xbest))

    else:
        for event_name in env.get_selected_event_names():
            env.set_current_event_name(event_name)
            env.setup_modelling()
            problem = env.get_problem()
            ds = env.get_dataset()
            xref = problem.preconstrain(problem.get_reference_model())
            payload.append((ds, problem, xref))

    all_trs = []
    events = []
    stations = {}
    for (ds, problem, x) in payload:
        results = problem.evaluate(x)

        event = problem.get_source(x).pyrocko_event()
        events.append(event)

        for result in results:
            if isinstance(result, WaveformMisfitResult):
                if show == 'filtered':
                    result.filtered_obs.set_codes(location='ob')
                    result.filtered_syn.set_codes(location='sy')
                    all_trs.append(result.filtered_obs)
                    all_trs.append(result.filtered_syn)
                elif show == 'processed':
                    result.processed_obs.set_codes(location='ob')
                    result.processed_syn.set_codes(location='sy')
                    all_trs.append(result.processed_obs)
                    all_trs.append(result.processed_syn)
                else:
                    raise ValueError('Invalid argument for show: %s' % show)

        for station in ds.get_stations():
            stations[station.nsl()] = station

    markers = []
    for ev in events:
        markers.append(pmarker.EventMarker(ev))

    trace.snuffle(all_trs, markers=markers, stations=list(stations.values()))


def harvest(
        rundir, problem=None, nbest=10, force=False, weed=0,
        export_fits=[]):

    env = Environment([rundir])
    optimiser = env.get_optimiser()
    nchains = env.get_optimiser().nchains

    if problem is None:
        problem, xs, misfits, bootstrap_misfits, _ = \
            load_problem_info_and_data(rundir, nchains=nchains)
    else:
        xs, misfits, bootstrap_misfits, _ = \
            load_problem_data(rundir, problem, nchains=nchains)

    logger.info('Harvesting problem "%s"...' % problem.name)

    dumpdir = op.join(rundir, 'harvest')
    if op.exists(dumpdir):
        if force:
            shutil.rmtree(dumpdir)
        else:
            raise DirectoryAlreadyExists(
                'Harvest directory already exists: %s' % dumpdir)

    util.ensuredir(dumpdir)

    ibests_list = []
    ibests = []
    gms = bootstrap_misfits[:, 0]
    isort = num.argsort(gms)

    ibests_list.append(isort[:nbest])

    if weed != 3:
        for ibootstrap in range(optimiser.nbootstrap):
            bms = bootstrap_misfits[:, ibootstrap]
            isort = num.argsort(bms)
            ibests_list.append(isort[:nbest])
            ibests.append(isort[0])

        if weed:
            mean_gm_best = num.median(gms[ibests])
            std_gm_best = num.std(gms[ibests])
            ibad = set()

            for ibootstrap, ibest in enumerate(ibests):
                if gms[ibest] > mean_gm_best + std_gm_best:
                    ibad.add(ibootstrap)

            ibests_list = [
                ibests_ for (ibootstrap, ibests_) in enumerate(ibests_list)
                if ibootstrap not in ibad]

    ibests = num.concatenate(ibests_list)

    if weed == 2:
        ibests = ibests[gms[ibests] < mean_gm_best]

    for i in ibests:
        problem.dump_problem_data(dumpdir, xs[i], misfits[i, :, :])

    if export_fits:
        env.setup_modelling()
        problem = env.get_problem()
        history = env.get_history(subset='harvest')

        for what in export_fits:
            if what == 'best':
                models = [history.get_best_model()]
            elif what == 'mean':
                models = [history.get_mean_model()]
            elif what == 'ensemble':
                models = history.models
            else:
                raise GrondError(
                    'Invalid option for harvest\'s export_fits argument: %s'
                    % what)

            results = []
            for x in models:
                results.append([
                    (result if isinstance(result, MisfitResult)
                     else MisfitResultError(message=str(result))) for
                    result in problem.evaluate(x)])

            emr = MisfitResultCollection(results=results)

            dump_misfit_result_collection(
                emr,
                op.join(dumpdir, 'fits-%s.yaml' % what))

    logger.info('Done harvesting problem "%s".' % problem.name)


def cluster(rundir, clustering, metric):
    env = Environment([rundir])
    history = env.get_history(subset='harvest')
    problem = history.problem
    models = history.models

    events = [problem.get_source(model).pyrocko_event() for model in models]

    from grond.clustering import metrics

    if metric not in metrics.metrics:
        raise GrondError('Unknown metric: %s' % metric)

    mat = metrics.compute_similarity_matrix(events, metric)

    clusters = clustering.perform(mat)

    labels = num.sort(num.unique(clusters))
    bins = num.concatenate((labels, [labels[-1]+1]))
    ns = num.histogram(clusters, bins)[0]

    history.set_attribute('cluster', clusters)

    for i in range(labels.size):
        if labels[i] == -1:
            logging.info(
                'Number of unclustered events: %5i' % ns[i])
        else:
            logging.info(
                'Number of events in cluster %i: %5i' % (labels[i], ns[i]))


def get_event_names(config):
    return config.get_event_names()


def check_problem(problem):
    if len(problem.targets) == 0:
        raise GrondError('No targets available')


def check(
        config,
        event_names=None,
        target_string_ids=None,
        show_waveforms=False,
        n_random_synthetics=10,
        stations_used_path=None):

    markers = []
    stations_used = {}
    erroneous = []
    for ievent, event_name in enumerate(event_names):
        ds = config.get_dataset(event_name)
        event = ds.get_event()
        trs_all = []
        try:
            problem = config.get_problem(event)

            _, nfamilies = problem.get_family_mask()
            logger.info('Problem: %s' % problem.name)
            logger.info('Number of target families: %i' % nfamilies)
            logger.info('Number of targets (total): %i' % len(problem.targets))

            if target_string_ids:
                problem.targets = [
                    target for target in problem.targets
                    if util.match_nslc(target_string_ids, target.string_id())]

            logger.info(
                'Number of targets (selected): %i' % len(problem.targets))

            check_problem(problem)

            results_list = []
            sources = []
            if n_random_synthetics == 0:
                x = problem.preconstrain(problem.get_reference_model())
                sources.append(problem.base_source)
                results = problem.evaluate(x)
                results_list.append(results)

            else:
                for i in range(n_random_synthetics):
                    x = problem.get_random_model()
                    sources.append(problem.get_source(x))
                    results = problem.evaluate(x)
                    results_list.append(results)

            if show_waveforms:
                engine = config.engine_config.get_engine()
                times = []
                tdata = []
                for target in problem.targets:
                    tobs_shift_group = []
                    tcuts = []
                    for source in sources:
                        tmin_fit, tmax_fit, tfade, tfade_taper = \
                            target.get_taper_params(engine, source)

                        times.extend((tmin_fit-tfade*2., tmax_fit+tfade*2.))

                        tobs, tsyn = target.get_pick_shift(engine, source)
                        if None not in (tobs, tsyn):
                            tobs_shift = tobs - tsyn
                        else:
                            tobs_shift = 0.0

                        tcuts.append(target.get_cutout_timespan(
                            tmin_fit+tobs_shift, tmax_fit+tobs_shift, tfade))

                        tobs_shift_group.append(tobs_shift)

                    tcuts = num.array(tcuts, dtype=num.float)

                    tdata.append((
                        tfade,
                        num.mean(tobs_shift_group),
                        (num.min(tcuts[:, 0]), num.max(tcuts[:, 1]))))

                tmin = min(times)
                tmax = max(times)

                tmax += (tmax-tmin)*2

                for (tfade, tobs_shift, tcut), target in zip(
                        tdata, problem.targets):

                    store = engine.get_store(target.store_id)

                    deltat = store.config.deltat

                    freqlimits = list(target.get_freqlimits())
                    freqlimits[2] = 0.45/deltat
                    freqlimits[3] = 0.5/deltat
                    freqlimits = tuple(freqlimits)

                    try:
                        trs_projected, trs_restituted, trs_raw, _ = \
                            ds.get_waveform(
                                target.codes,
                                tmin=tmin+tobs_shift,
                                tmax=tmax+tobs_shift,
                                tfade=tfade,
                                freqlimits=freqlimits,
                                deltat=deltat,
                                backazimuth=target.
                                get_backazimuth_for_waveform(),
                                debug=True)

                    except NotFound as e:
                        logger.warn(str(e))
                        continue

                    trs_projected = copy.deepcopy(trs_projected)
                    trs_restituted = copy.deepcopy(trs_restituted)
                    trs_raw = copy.deepcopy(trs_raw)

                    for trx in trs_projected + trs_restituted + trs_raw:
                        trx.shift(-tobs_shift)
                        trx.set_codes(
                            network='',
                            station=target.string_id(),
                            location='')

                    for trx in trs_projected:
                        trx.set_codes(location=trx.location + '2_proj')

                    for trx in trs_restituted:
                        trx.set_codes(location=trx.location + '1_rest')

                    for trx in trs_raw:
                        trx.set_codes(location=trx.location + '0_raw')

                    trs_all.extend(trs_projected)
                    trs_all.extend(trs_restituted)
                    trs_all.extend(trs_raw)

                    for source in sources:
                        tmin_fit, tmax_fit, tfade, tfade_taper = \
                            target.get_taper_params(engine, source)

                        markers.append(pmarker.Marker(
                            nslc_ids=[('', target.string_id(), '*_proj', '*')],
                            tmin=tmin_fit, tmax=tmax_fit))

                    markers.append(pmarker.Marker(
                        nslc_ids=[('', target.string_id(), '*_raw', '*')],
                        tmin=tcut[0]-tobs_shift, tmax=tcut[1]-tobs_shift,
                        kind=1))

            else:
                for itarget, target in enumerate(problem.targets):

                    nok = 0
                    for results in results_list:
                        result = results[itarget]
                        if not isinstance(result, gf.SeismosizerError):
                            nok += 1

                    if nok == 0:
                        sok = 'not used'
                    elif nok == len(results_list):
                        sok = 'ok'
                        try:
                            s = ds.get_station(target)
                            stations_used[s.nsl()] = s
                        except (NotFound, InvalidObject):
                            pass
                    else:
                        sok = 'not used (%i/%i ok)' % (nok, len(results_list))

                    logger.info('%-40s %s' % (
                        (target.string_id() + ':', sok)))

        except GrondError as e:
            logger.error('Event %i, "%s": %s' % (
                ievent,
                event.name or util.time_to_str(event.time),
                str(e)))

            erroneous.append(event)

        if show_waveforms:
            trace.snuffle(trs_all, stations=ds.get_stations(), markers=markers)

        if stations_used_path:
            stations = list(stations_used.values())
            stations.sort(key=lambda s: s.nsl())
            model.dump_stations(stations, stations_used_path)

    if erroneous:
        raise GrondError(
            'Check failed for events: %s'
            % ', '.join(ev.name for ev in erroneous))


g_state = {}


def go(environment,
       force=False, preserve=False,
       nparallel=1, status='state', nthreads=0):

    g_data = (environment, force, preserve,
              status, nparallel, nthreads)
    g_state[id(g_data)] = g_data

    nevents = environment.nevents_selected
    for x in parimap.parimap(
            process_event,
            range(environment.nevents_selected),
            [id(g_data)] * nevents,
            nprocs=nparallel):

        pass


def process_event(ievent, g_data_id):

    environment, force, preserve, status, nparallel, nthreads = \
        g_state[g_data_id]

    config = environment.get_config()
    event_name = environment.get_selected_event_names()[ievent]
    nevents = environment.nevents_selected
    tstart = time.time()

    ds = config.get_dataset(event_name)
    event = ds.get_event()
    problem = config.get_problem(event)

    synt = ds.synthetic_test
    if synt:
        problem.base_source = problem.get_source(synt.get_x())

    check_problem(problem)

    rundir = expand_template(
        config.rundir_template,
        dict(problem_name=problem.name))
    environment.set_rundir_path(rundir)

    if op.exists(rundir):
        if preserve:
            nold_rundirs = len(glob.glob(rundir + '*'))
            shutil.move(rundir, rundir+'-old-%d' % (nold_rundirs))
        elif force:
            shutil.rmtree(rundir)
        else:
            logger.warn('Skipping problem "%s": rundir already exists: %s' %
                        (problem.name, rundir))
            return

    util.ensuredir(rundir)

    logger.info(
        'Starting event %i / %i' % (ievent+1, nevents))

    logger.info('Rundir: %s' % rundir)

    logger.info('Analysing problem "%s".' % problem.name)

    for analyser_conf in config.analyser_configs:
        analyser = analyser_conf.get_analyser()
        analyser.analyse(problem, ds)

    basepath = config.get_basepath()
    config.change_basepath(rundir)
    guts.dump(config, filename=op.join(rundir, 'config.yaml'))
    config.change_basepath(basepath)

    optimiser = config.optimiser_config.get_optimiser()
    optimiser.set_nthreads(nthreads)

    optimiser.init_bootstraps(problem)
    problem.dump_problem_info(rundir)

    monitor = None

    xs_inject = None
    synt = ds.synthetic_test
    if synt and synt.inject_solution:
        xs_inject = synt.get_x()[num.newaxis, :]

    try:
        if xs_inject is not None:
            from .optimisers import highscore
            if not isinstance(optimiser, highscore.HighScoreOptimiser):
                raise GrondError(
                    'Optimiser does not support injections.')

            optimiser.sampler_phases[0:0] = [
                highscore.InjectionSamplerPhase(xs_inject=xs_inject)]

        with lock_rundir(rundir):
            if status == 'state':
                monitor = GrondMonitor.watch(rundir)
            optimiser.optimise(
                problem,
                rundir=rundir)

        harvest(rundir, problem, force=True)

    except BadProblem as e:
        logger.error(str(e))

    except GrondError as e:
        logger.error(str(e))

    finally:
        if monitor:
            monitor.terminate()

    tstop = time.time()
    logger.info(
        'Stop %i / %i (%g min)' % (ievent+1, nevents, (tstop - tstart)/60.))

    logger.info(
        'Done with problem "%s", rundir is "%s".' % (problem.name, rundir))


class ParameterStats(Object):
    name = String.T()
    mean = Float.T()
    std = Float.T()
    best = Float.T()
    minimum = Float.T()
    percentile5 = Float.T()
    percentile16 = Float.T()
    median = Float.T()
    percentile84 = Float.T()
    percentile95 = Float.T()
    maximum = Float.T()

    def __init__(self, *args, **kwargs):
        kwargs.update(zip(self.T.propnames, args))
        Object.__init__(self, **kwargs)

    def get_values_dict(self):
        return dict(
            (self.name+'.' + k, getattr(self, k))
            for k in self.T.propnames
            if k != 'name')


class ResultStats(Object):
    problem = Problem.T()
    parameter_stats_list = List.T(ParameterStats.T())

    def get_values_dict(self):
        d = {}
        for ps in self.parameter_stats_list:
            d.update(ps.get_values_dict())
        return d


def make_stats(problem, models, gms, pnames=None):
    ibest = num.argmin(gms)
    rs = ResultStats(problem=problem)
    if pnames is None:
        pnames = problem.parameter_names

    for pname in pnames:
        iparam = problem.name_to_index(pname)
        vs = problem.extract(models, iparam)
        mi, p5, p16, median, p84, p95, ma = map(float, num.percentile(
            vs, [0., 5., 16., 50., 84., 95., 100.]))

        mean = float(num.mean(vs))
        std = float(num.std(vs))
        best = float(vs[ibest])
        s = ParameterStats(
            pname, mean, std, best, mi, p5, p16, median, p84, p95, ma)

        rs.parameter_stats_list.append(s)

    return rs


def try_add_location_uncertainty(data, types):
    vs = [data.get(k, None) for k in (
        'north_shift.std', 'east_shift.std', 'depth.std')]

    if None not in vs:
        data['location_uncertainty'] = math.sqrt(sum(v**2 for v in vs))
        types['location_uncertainty'] = float


def format_stats(rs, fmt):
    pname_to_pindex = dict(
        (p.name, i) for (i, p) in enumerate(rs.parameter_stats_list))

    values = []
    headers = []
    for x in fmt:
        pname, qname = x.split('.')
        pindex = pname_to_pindex[pname]
        values.append(getattr(rs.parameter_stats_list[pindex], qname))
        headers.append(x)

    return ' '.join('%16.7g' % v for v in values)


def export(
        what, rundirs, type=None, pnames=None, filename=None, selection=None):

    if pnames is not None:
        pnames_clean = [pname.split('.')[0] for pname in pnames]
        shortform = all(len(pname.split('.')) == 2 for pname in pnames)
    else:
        pnames_clean = None
        shortform = False

    if what == 'stats' and type is not None:
        raise GrondError('Invalid argument combination: what=%s, type=%s' % (
            repr(what), repr(type)))

    if what != 'stats' and shortform:
        raise GrondError('Invalid argument combination: what=%s, pnames=%s' % (
            repr(what), repr(pnames)))

    if what != 'stats' and type != 'vector' and pnames is not None:
        raise GrondError(
            'Invalid argument combination: what=%s, type=%s, pnames=%s' % (
                repr(what), repr(type), repr(pnames)))

    if filename is None:
        out = sys.stdout
    else:
        out = open(filename, 'w')

    if type is None:
        type = 'event'

    if shortform:
        print('#', ' '.join(['%16s' % x for x in pnames]), file=out)

    def dump(x, gm, indices):
        if type == 'vector':
            print(' ', ' '.join(
                '%16.7g' % problem.extract(x, i) for i in indices),
                '%16.7g' % gm, file=out)

        elif type == 'source':
            source = problem.get_source(x)
            guts.dump(source, stream=out)

        elif type == 'event':
            ev = problem.get_source(x).pyrocko_event()
            model.dump_events([ev], stream=out)

        elif type == 'event-yaml':
            ev = problem.get_source(x).pyrocko_event()
            guts.dump_all([ev], stream=out)

        else:
            raise GrondError('Invalid argument: type=%s' % repr(type))

    header = None
    for rundir in rundirs:
        env = Environment(rundir)
        info = env.get_run_info()

        try:
            history = env.get_history(subset='harvest')
        except ProblemDataNotAvailable as e:
            logger.error(
                'Harvest not available (Did the run succeed?): %s' % str(e))
            continue

        problem = history.problem
        models = history.models
        misfits = history.get_primary_chain_misfits()

        if selection:
            rs = make_stats(
                problem, models,
                history.get_primary_chain_misfits())

            data = dict(tags=info.tags)
            types = dict(tags=(list, str))

            for k, v in rs.get_values_dict().items():
                data[k] = v
                types[k] = float

            try_add_location_uncertainty(data, types)

            if not selected(selection, data=data, types=types):
                continue

        else:
            rs = None

        if type == 'vector':
            pnames_take = pnames_clean or \
                problem.parameter_names[:problem.nparameters]

            indices = num.array(
                [problem.name_to_index(pname) for pname in pnames_take])

            if type == 'vector' and what in ('best', 'mean', 'ensemble'):
                extra = ['global_misfit']
            else:
                extra = []

            new_header = '# ' + ' '.join(
                '%16s' % x for x in pnames_take + extra)

            if type == 'vector' and header != new_header:
                print(new_header, file=out)

            header = new_header
        else:
            indices = None

        if what == 'best':
            x_best = history.get_best_model()
            gm_best = history.get_best_misfit()
            dump(x_best, gm_best, indices)

        elif what == 'mean':
            x_mean = history.get_mean_model()
            gm_mean = history.get_mean_misfit(chain=0)
            dump(x_mean, gm_mean, indices)

        elif what == 'ensemble':
            isort = num.argsort(misfits)
            for i in isort:
                dump(models[i], misfits[i], indices)

        elif what == 'stats':
            if not rs:
                rs = make_stats(problem, models,
                                history.get_primary_chain_misfits(),
                                pnames_clean)

            if shortform:
                print(' ', format_stats(rs, pnames), file=out)
            else:
                print(rs, file=out)

        else:
            raise GrondError('Invalid argument: what=%s' % repr(what))

    if out is not sys.stdout:
        out.close()


__all__ = '''
    DirectoryAlreadyExists
    forward
    harvest
    cluster
    go
    get_event_names
    check
    export
'''.split()
