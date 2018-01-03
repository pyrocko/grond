from __future__ import print_function

import os
import sys
import logging
import time
import copy
import shutil
import os.path as op
import numpy as num

from pyrocko.guts import load, Object, String, Float, Bool, List
from pyrocko import orthodrome as od, gf, trace, guts, util, weeding
from pyrocko import parimap, model, marker as pmarker

from .dataset import DatasetConfig, NotFound
from .problems.base import ProblemConfig, Problem, \
    load_problem_info_and_data, load_problem_data

from .optimizers.base import OptimizerConfig, BadProblem
from .targets.base import TargetGroup
from .analysers.base import AnalyserConfig
from .meta import Path, HasPaths, expand_template, GrondError, Forbidden

logger = logging.getLogger('grond.core')
guts_prefix = 'grond'


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


class DirectoryAlreadyExists(Exception):
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


class EngineConfig(HasPaths):
    gf_stores_from_pyrocko_config = Bool.T(default=True)
    gf_store_superdirs = List.T(Path.T())
    gf_store_dirs = List.T(Path.T())

    def __init__(self, *args, **kwargs):
        HasPaths.__init__(self, *args, **kwargs)
        self._engine = None

    def get_engine(self):
        if self._engine is None:
            fp = self.expand_path
            self._engine = gf.LocalEngine(
                use_config=self.gf_stores_from_pyrocko_config,
                store_superdirs=fp(self.gf_store_superdirs),
                store_dirs=fp(self.gf_store_dirs))

        return self._engine


class Config(HasPaths):
    rundir_template = Path.T()
    dataset_config = DatasetConfig.T()
    target_groups = List.T(TargetGroup.T())
    problem_config = ProblemConfig.T()
    analyser_config = AnalyserConfig.T(default=AnalyserConfig.D())
    optimizer_config = OptimizerConfig.T()
    engine_config = EngineConfig.T(default=EngineConfig.D())

    def __init__(self, *args, **kwargs):
        HasPaths.__init__(self, *args, **kwargs)

    def get_event_names(self):
        return self.dataset_config.get_event_names()

    def get_dataset(self, event_name):
        return self.dataset_config.get_dataset(event_name)

    def get_targets(self, event):
        ds = self.get_dataset(event.name)

        targets = []
        for igroup, target_group in enumerate(self.target_groups):
            targets.extend(target_group.get_targets(
                ds, event, 'target.%i' % igroup))

        return targets

    def setup_modelling_environment(self, problem):
        problem.set_engine(self.engine_config.get_engine())
        ds = self.get_dataset(problem.base_source.name)
        synt = ds.synthetic_test
        if synt:
            synt.set_problem(problem)
            problem.base_source = problem.get_source(synt.get_x())

    def get_problem(self, event):
        targets = self.get_targets(event)
        problem = self.problem_config.get_problem(event, targets)
        self.setup_modelling_environment(problem)
        return problem


def sarr(a):
    return ' '.join('%15g' % x for x in a)


def get_mean_x(xs):
    return num.mean(xs, axis=0)


def get_mean_x_and_gm(problem, xs, misfits):
    gms = problem.combine_misfits(misfits)
    return num.mean(xs, axis=0), num.mean(gms)


def get_best_x(problem, xs, misfits):
    gms = problem.combine_misfits(misfits)
    ibest = num.argmin(gms)
    return xs[ibest, :]


def get_best_x_and_gm(problem, xs, misfits):
    gms = problem.combine_misfits(misfits)
    ibest = num.argmin(gms)
    return xs[ibest, :], gms[ibest]


def get_mean_source(problem, xs):
    x_mean = get_mean_x(xs)
    source = problem.get_source(x_mean)
    return source


def get_best_source(problem, xs, misfits):
    x_best = get_best_x(problem, xs, misfits)
    source = problem.get_source(x_best)
    return source


def mean_latlondist(lats, lons):
    if len(lats) == 0:
        return 0., 0., 1000.
    else:
        ns, es = od.latlon_to_ne_numpy(lats[0], lons[0], lats, lons)
        n, e = num.mean(ns), num.mean(es)
        dists = num.sqrt((ns-n)**2 + (es-e)**2)
        lat, lon = od.ne_to_latlon(lats[0], lons[0], n, e)
        return float(lat), float(lon), float(num.max(dists))


def stations_mean_latlondist(stations):
    lats = num.array([s.lat for s in stations])
    lons = num.array([s.lon for s in stations])
    return mean_latlondist(lats, lons)


def read_config(path):
    config = load(filename=path)
    if not isinstance(config, Config):
        raise GrondError('invalid Grond configuration in file "%s"' % path)

    config.set_basepath(op.dirname(path) or '.')
    return config


def write_config(config, path):
    basepath = config.get_basepath()
    dirname = op.dirname(path) or '.'
    config.change_basepath(dirname)
    guts.dump(config, filename=path)
    config.change_basepath(basepath)


def forward(rundir_or_config_path, event_names):

    if not event_names:
        return

    if os.path.isdir(rundir_or_config_path):
        rundir = rundir_or_config_path
        config = guts.load(
            filename=op.join(rundir, 'config.yaml'))

        config.set_basepath(rundir)
        problem, xs, misfits = load_problem_info_and_data(
            rundir, subset='harvest')

        gms = problem.combine_misfits(misfits)
        ibest = num.argmin(gms)
        xbest = xs[ibest, :]

        ds = config.get_dataset(problem.base_source.name)
        problem.set_engine(config.engine_config.get_engine())

        for target in problem.targets:
            target.set_dataset(ds)

        payload = [(problem, xbest)]

    else:
        config = read_config(rundir_or_config_path)

        payload = []
        for event_name in event_names:
            ds = config.get_dataset(event_name)
            event = ds.get_event()
            problem = config.get_problem(event)
            xref = problem.preconstrain(
                problem.pack(problem.base_source))
            payload.append((problem, xref))

    all_trs = []
    events = []
    for (problem, x) in payload:
        ds.empty_cache()
        _, results = problem.evaluate(x, result_mode='full')

        event = problem.get_source(x).pyrocko_event()
        events.append(event)

        for result in results:
            if not isinstance(result, gf.SeismosizerError):
                result.filtered_obs.set_codes(location='ob')
                result.filtered_syn.set_codes(location='sy')
                all_trs.append(result.filtered_obs)
                all_trs.append(result.filtered_syn)

    markers = []
    for ev in events:
        markers.append(pmarker.EventMarker(ev))

    trace.snuffle(all_trs, markers=markers, stations=ds.get_stations())


def harvest(rundir, problem=None, nbest=10, force=False, weed=0):

    if problem is None:
        problem, xs, misfits = load_problem_info_and_data(rundir)
    else:
        xs, misfits = load_problem_data(rundir, problem)

    optimizer_fn = op.join(rundir, 'optimizer.yaml')
    optimizer = guts.load(filename=optimizer_fn)

    dumpdir = op.join(rundir, 'harvest')
    if op.exists(dumpdir):
        if force:
            shutil.rmtree(dumpdir)
        else:
            raise DirectoryAlreadyExists(dumpdir)

    util.ensuredir(dumpdir)

    ibests_list = []
    ibests = []
    gms = problem.combine_misfits(misfits)
    isort = num.argsort(gms)

    ibests_list.append(isort[:nbest])

    if weed != 3:
        for ibootstrap in range(optimizer.nbootstrap):
            bms = optimizer.bootstrap_misfits(problem, misfits, ibootstrap)
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


def get_event_names(config):
    return config.get_event_names()


def check_problem(problem):
    if len(problem.targets) == 0:
        raise GrondError('no targets available')


def check(
        config,
        event_names=None,
        target_string_ids=None,
        show_plot=False,
        show_waveforms=False,
        n_random_synthetics=10):

    if show_plot:
        from matplotlib import pyplot as plt
        from grond.plot import colors

    markers = []
    for ievent, event_name in enumerate(event_names):
        ds = config.get_dataset(event_name)
        event = ds.get_event()
        trs_all = []
        try:
            problem = config.get_problem(event)

            _, nfamilies = problem.get_family_mask()
            logger.info('number of target families: %i' % nfamilies)
            logger.info('number of targets (total): %i' % len(problem.targets))

            if target_string_ids:
                problem.targets = [
                    target for target in problem.targets
                    if util.match_nslc(target_string_ids, target.string_id())]

            logger.info(
                'number of targets (selected): %i' % len(problem.targets))

            check_problem(problem)

            xbounds = num.array(
                problem.get_parameter_bounds(), dtype=num.float)

            results_list = []

            sources = []
            if n_random_synthetics == 0:
                x = problem.pack(problem.base_source)
                sources.append(problem.base_source)
                _, results = problem.evaluate(x, result_mode='full')
                results_list.append(results)

            else:
                for i in range(n_random_synthetics):
                    while True:
                        x = problem.random_uniform(xbounds)
                        try:
                            x = problem.preconstrain(x)
                            break

                        except Forbidden:
                            pass

                    sources.append(problem.get_source(x))
                    _, results = problem.evaluate(x, result_mode='full')
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
                        trs_projected, trs_restituted, trs_raw = \
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
                            nslc_ids=[('', target.string_id(), '*', '*')],
                            tmin=tmin_fit, tmax=tmax_fit))

                    markers.append(pmarker.Marker(
                        nslc_ids=[('', target.string_id(), '*', '*')],
                        tmin=tcut[0]-tobs_shift, tmax=tcut[1]-tobs_shift,
                        kind=1))

            if show_plot:
                for itarget, target in enumerate(problem.targets):
                    yabsmaxs = []
                    for results in results_list:
                        result = results[itarget]
                        if not isinstance(result, gf.SeismosizerError):
                            yabsmaxs.append(
                                num.max(num.abs(
                                    result.filtered_obs.get_ydata())))

                    if yabsmaxs:
                        yabsmax = max(yabsmaxs) or 1.0
                    else:
                        yabsmax = None

                    fig = None
                    ii = 0
                    for results in results_list:
                        result = results[itarget]
                        if not isinstance(result, gf.SeismosizerError):
                            if fig is None:
                                fig = plt.figure()
                                axes = fig.add_subplot(1, 1, 1)
                                axes.set_ylim(0., 4.)
                                axes.set_title('%s' % target.string_id())

                            xdata = result.filtered_obs.get_xdata()
                            ydata = result.filtered_obs.get_ydata() / yabsmax
                            axes.plot(xdata, ydata*0.5 + 3.5, color='black')

                            color = colors[ii % len(colors)]

                            xdata = result.filtered_syn.get_xdata()
                            ydata = result.filtered_syn.get_ydata()
                            ydata = ydata / (num.max(num.abs(ydata)) or 1.0)

                            axes.plot(xdata, ydata*0.5 + 2.5, color=color)

                            xdata = result.processed_syn.get_xdata()
                            ydata = result.processed_syn.get_ydata()
                            ydata = ydata / (num.max(num.abs(ydata)) or 1.0)

                            axes.plot(xdata, ydata*0.5 + 1.5, color=color)
                            if result.tsyn_pick:
                                axes.axvline(
                                    result.tsyn_pick,
                                    color=(0.7, 0.7, 0.7),
                                    zorder=2)

                            t = result.processed_syn.get_xdata()
                            taper = result.taper

                            y = num.ones(t.size) * 0.9
                            taper(y, t[0], t[1] - t[0])
                            y2 = num.concatenate((y, -y[::-1]))
                            t2 = num.concatenate((t, t[::-1]))
                            axes.plot(t2, y2 * 0.5 + 0.5, color='gray')
                            ii += 1
                        else:
                            logger.info(str(result))

                    if fig:
                        plt.show()

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
                    else:
                        sok = 'not used (%i/%i ok)' % (nok, len(results_list))

                    logger.info('%-40s %s' % (
                        (target.string_id() + ':', sok)))

        except GrondError as e:
            logger.error('event %i, %s: %s' % (
                ievent,
                event.name or util.time_to_str(event.time),
                str(e)))

        if show_waveforms:
            trace.snuffle(trs_all, stations=ds.get_stations(), markers=markers)


g_state = {}


def go(config, event_names=None, force=False, nparallel=1, status=('state',)):

    status = tuple(status)
    g_data = (config, force, status, nparallel, event_names)
    g_state[id(g_data)] = g_data

    nevents = len(event_names)

    for x in parimap.parimap(
            process_event,
            range(nevents),
            [id(g_data)] * nevents,
            nprocs=nparallel):

        pass


def process_event(ievent, g_data_id):

    config, force, status, nparallel, event_names = g_state[g_data_id]

    event_name = event_names[ievent]

    ds = config.get_dataset(event_name)

    nevents = len(event_names)

    tstart = time.time()

    event = ds.get_event()

    problem = config.get_problem(event)

    synt = ds.synthetic_test
    if synt:
        problem.base_source = problem.get_source(synt.get_x())

    check_problem(problem)

    rundir = expand_template(
        config.rundir_template,
        dict(problem_name=problem.name))

    if op.exists(rundir):
        if force:
            shutil.rmtree(rundir)
        else:
            logger.warn('skipping problem %s: rundir already exists: %s' %
                        (problem.name, rundir))
            return

    util.ensuredir(rundir)

    logger.info(
        'start %i / %i' % (ievent+1, nevents))

    analyser = config.analyser_config.get_analyser()
    analyser.analyse(problem)

    basepath = config.get_basepath()
    config.change_basepath(rundir)
    guts.dump(config, filename=op.join(rundir, 'config.yaml'))
    config.change_basepath(basepath)

    problem.dump_problem_info(rundir)

    xs_inject = None
    synt = ds.synthetic_test
    if synt and synt.inject_solution:
        xs_inject = synt.get_x()[num.newaxis, :]

    # from matplotlib import pyplot as plt
    # from grond import plot
    # splot = plot.SolverPlot(
    #     plt, 'time', 'magnitude',
    #     show=False,
    #     update_every=10,
    #     movie_filename='grond_opt_time_magnitude.mp4')

    def startThreads():
        from .listeners import terminal
        term = terminal.TerminalListener(rundir)
        term.start()
        return term

    term = startThreads()

    try:
        optimizer = config.optimizer_config.get_optimizer()
        if xs_inject is not None:
            from .optimizers import highscore
            if not isinstance(optimizer, highscore.HighScoreOptimizer):
                raise GrondError(
                    'optimizer does not support injections')

            optimizer.sampler_phases[0:0] = [
                highscore.InjectionSamplerPhase(xs_inject=xs_inject)]

        optimizer.optimize(
            problem,
            rundir=rundir)

        harvest(rundir, problem, force=True)

    except BadProblem as e:
        logger.error(str(e))
    finally:
        term.join()

    tstop = time.time()
    logger.info(
        'stop %i / %i (%g min)' % (ievent, nevents, (tstop - tstart)/60.))

    logger.info(
        'done with problem %s, rundir is %s' % (problem.name, rundir))


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


class ResultStats(Object):
    problem = Problem.T()
    parameter_stats_list = List.T(ParameterStats.T())


def make_stats(problem, xs, misfits, pnames=None):
    gms = problem.combine_misfits(misfits)
    ibest = num.argmin(gms)
    rs = ResultStats(problem=problem)
    if pnames is None:
        pnames = problem.parameter_names

    for pname in pnames:
        iparam = problem.name_to_index(pname)
        vs = problem.extract(xs, iparam)
        mi, p5, p16, median, p84, p95, ma = map(float, num.percentile(
            vs, [0., 5., 16., 50., 84., 95., 100.]))

        mean = float(num.mean(vs))
        std = float(num.std(vs))
        best = float(vs[ibest])
        s = ParameterStats(
            pname, mean, std, best, mi, p5, p16, median, p84, p95, ma)

        rs.parameter_stats_list.append(s)

    return rs


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


def export(what, rundirs, type=None, pnames=None, filename=None):
    if pnames is not None:
        pnames_clean = [pname.split('.')[0] for pname in pnames]
        shortform = all(len(pname.split('.')) == 2 for pname in pnames)
    else:
        pnames_clean = None
        shortform = False

    if what == 'stats' and type is not None:
        raise GrondError('invalid argument combination: what=%s, type=%s' % (
            repr(what), repr(type)))

    if what != 'stats' and shortform:
        raise GrondError('invalid argument combination: what=%s, pnames=%s' % (
            repr(what), repr(pnames)))

    if what != 'stats' and type != 'vector' and pnames is not None:
        raise GrondError(
            'invalid argument combination: what=%s, type=%s, pnames=%s' % (
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

        else:
            raise GrondError('invalid argument: type=%s' % repr(type))

    header = None
    for rundir in rundirs:
        problem, xs, misfits = load_problem_info_and_data(
            rundir, subset='harvest')

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
            x_best, gm_best = get_best_x_and_gm(problem, xs, misfits)
            dump(x_best, gm_best, indices)

        elif what == 'mean':
            x_mean, gm_mean = get_mean_x_and_gm(problem, xs, misfits)
            dump(x_mean, gm_mean, indices)

        elif what == 'ensemble':
            gms = problem.combine_misfits(misfits)
            isort = num.argsort(gms)
            for i in isort:
                dump(xs[i], gms[i], indices)

        elif what == 'stats':
            rs = make_stats(problem, xs, misfits, pnames_clean)
            if shortform:
                print(' ', format_stats(rs, pnames), file=out)
            else:
                print(rs, file=out)

        else:
            raise GrondError('invalid argument: what=%s' % repr(what))

    if out is not sys.stdout:
        out.close()


__all__ = '''
    EngineConfig
    Config
    read_config
    write_config
    forward
    harvest
    go
    get_event_names
    check
    export
'''.split()
