#!/usr/bin/env python

import math
import sys
import shutil
import os
import os.path as op
import logging
from optparse import OptionParser

from pyrocko import util, marker
from pyrocko.gf import Range

import grond

logger = logging.getLogger('grond.main')
km = 1e3


def d2u(d):
    if isinstance(d, dict):
        return dict((k.replace('-', '_'), v) for (k, v) in d.items())
    else:
        return d.replace('-', '_')


subcommand_descriptions = {
    'scenario': 'create an scenario project',
    'init': 'create project structure or print example configuration',
    'events': 'print available event names for given configuration',
    'check': 'check data and configuration',
    'go': 'run Grond optimization',
    'forward': 'run forward modelling',
    'harvest': 'manually run harvesting',
    'plot': 'plot optimization result',
    'movie': 'visualize optimizer evolution',
    'export': 'export results',
    'report': 'create result report',
    'report-index': 'create report index',
    'qc-polarization': 'check sensor orientations with polarization analysis',
}

subcommand_usages = {
    'scenario': 'scenario [options] <project_dir>',
    'init': 'init [options] <project_dir>',
    'events': 'events <configfile>',
    'check': 'check <configfile> <eventnames> ... [options]',
    'go': 'go <configfile> <eventnames> ... [options]',
    'forward': (
        'forward <rundir> [options]',
        'forward <configfile> <eventnames> ... [options]'),
    'harvest': 'harvest <rundir> [options]',
    'plot': 'plot <plotnames> <rundir> [options]',
    'movie': 'movie <rundir> <xpar> <ypar> <filetemplate> [options]',
    'export': 'export (best|mean|ensemble|stats) <rundirs> ... [options]',
    'report': (
        'report <rundir> ... [options]',
        'report <configfile> <eventnames> ...'),
    'report-index': 'report <reportdir> [options]',
    'qc-polarization': 'qc-polarization <configfile> <eventname> '
                       '<targetconfigname> [options]',
}

subcommands = subcommand_descriptions.keys()

program_name = 'grond'

usage_tdata = d2u(subcommand_descriptions)
usage_tdata['program_name'] = program_name

usage = '''%(program_name)s <subcommand> [options] [--] <arguments> ...

Subcommands:

    scenario         %(scenario)s
    init            %(init)s
    events          %(events)s
    check           %(check)s
    go              %(go)s
    forward         %(forward)s
    harvest         %(harvest)s
    plot            %(plot)s
    movie           %(movie)s
    export          %(export)s
    report          %(report)s
    report-index    %(report_index)s
    qc-polarization %(qc_polarization)s

To get further help and a list of available options for any subcommand run:

    %(program_name)s <subcommand> --help

''' % usage_tdata


def main(args=None):
    if not args:
        args = sys.argv

    args = list(sys.argv)
    if len(args) < 2:
        sys.exit('Usage: %s' % usage)

    args.pop(0)
    command = args.pop(0)

    if command in subcommands:
        globals()['command_' + d2u(command)](args)

    elif command in ('--help', '-h', 'help'):
        if command == 'help' and args:
            acommand = args[0]
            if acommand in subcommands:
                globals()['command_' + acommand](['--help'])

        sys.exit('Usage: %s' % usage)

    else:
        die('no such subcommand: %s' % command)


def add_common_options(parser):
    parser.add_option(
        '--loglevel',
        action='store',
        dest='loglevel',
        type='choice',
        choices=('critical', 'error', 'warning', 'info', 'debug'),
        default='info',
        help='set logger level to '
             '"critical", "error", "warning", "info", or "debug". '
             'Default is "%default".')


def process_common_options(options):
    util.setup_logging(program_name, options.loglevel)


def cl_parse(command, args, setup=None, details=None):
    usage = subcommand_usages[command]
    descr = subcommand_descriptions[command]

    if isinstance(usage, str):
        usage = [usage]

    susage = '%s %s' % (program_name, usage[0])
    for s in usage[1:]:
        susage += '\n%s%s %s' % (' '*7, program_name, s)

    description = descr[0].upper() + descr[1:] + '.'

    if details:
        description = description + '\n\n%s' % details

    parser = OptionParser(usage=susage, description=description)

    if setup:
        setup(parser)

    add_common_options(parser)
    (options, args) = parser.parse_args(args)
    process_common_options(options)
    return parser, options, args


def die(message, err=''):
    if err:
        sys.exit('%s: error: %s \n %s' % (program_name, message, err))
    else:
        sys.exit('%s: error: %s' % (program_name, message))


def help_and_die(parser, message):
    parser.print_help(sys.stderr)
    sys.stderr.write('\n')
    die(message)


def command_scenario(args):

    STORE_STATIC = 'ak135_static'
    STORE_WAVEFORMS = 'global_2s'

    def setup(parser):
        parser.add_option(
            '--waveforms', dest='waveforms', action='store_true',
            help='Add waveform configuration. '
                 '(default)')
        parser.add_option(
            '--insar', dest='insar', action='store_true',
            help='Add InSAR displacement scenes using kite containers. '
                 '(see https://pyrocko.org)')
        parser.add_option(
            '--nstations', dest='nstations', type=int, default=20,
            help='Number of seismic stations to create (default: %default)')
        parser.add_option(
            '--nevents', dest='nevents', type=int, default=1,
            help='Number of events to create (default: %default)')
        parser.add_option(
            '--lat', dest='lat', type=float, default=29.58,
            help='Center latitude of the scenario (default: %default)')
        parser.add_option(
            '--lon', dest='lon', type=float, default=81.72,
            help='Center latitude of the scenario (default: %default)')
        parser.add_option(
            '--radius', dest='radius', type=float, default=200.,
            help='Radius of the the scenario in [km] (default: %default)')
        parser.add_option(
            '--source', dest='source', type=str, default='dc',
            help='Source to generate \'dc\' (double couple)'
                 ' or\'rectangular\' (rectangular finite fault)'
                 ' (default: \'%default\')')
        parser.add_option(
            '--gf_waveforms', dest='store_waveforms', type=str,
            default=STORE_WAVEFORMS,
            help='Green\'s function store for waveform modelling, '
                 '(default: %default)')
        parser.add_option(
            '--gf_static', dest='store_statics', type=str,
            default=STORE_STATIC,
            help='Green\'s function store for static modelling, '
                 '(default: %default)')
        parser.add_option(
            '--force', dest='force', action='store_true',
            help='Overwrite existing project folder.')

    parser, options, args = cl_parse('init', args, setup)

    if len(args) == 1:
        project_dir = op.join(op.curdir, args[0])
    else:
        parser.print_help()
        sys.exit(1)

    from grond import scenario as grond_scenario
    scenario = grond_scenario.GrondScenario(
        project_dir,
        center_lat=options.lat, center_lon=options.lon,
        radius=options.radius*km)

    if not options.waveforms and not options.insar:
        options.waveforms = True

    if options.waveforms:
        obs = grond_scenario.WaveformObservation(
            nstations=options.nstations,
            store_id=options.store_waveforms)
        scenario.add_observation(obs)

    if options.insar:
        obs = grond_scenario.InSARObservation(
            store_id=options.store_statics)
        scenario.add_observation(obs)

    if options.source == 'dc':
        problem = grond_scenario.DCSourceProblem(
            nevents=options.nevents)
    elif options.source == 'rectangular':
        problem = grond_scenario.RectangularSourceProblem(
            nevents=options.nevents)
    scenario.set_problem(problem)

    scenario.build(force=options.force, interactive=True)


def command_init(args):

    def setup(parser):
        parser.add_option(
            '--waveforms', dest='waveforms', action='store_true',
            help='Add waveform configuration. '
                 '(default)')
        parser.add_option(
            '--insar', dest='insar', action='store_true',
            help='Add InSAR displacement scenes using kite containers. '
                 '(https://pyrocko.org')
        parser.add_option(
            '--force', dest='force', action='store_true',
            help='Overwrite existing project folder.')

    parser, options, args = cl_parse('init', args, setup)

    project_dir = None
    if len(args) == 1:
        project_dir = op.join(op.curdir, args[0])
        if op.exists(project_dir) and not options.force:
            raise EnvironmentError(
                'Directory %s already exists! Use --force to overwrite'
                % args[0])
        elif op.exists(project_dir) and options.force:
            logger.info('Overwriting directory %s.' % project_dir)
            shutil.rmtree(project_dir)

    if not options.insar and not options.waveforms:
        options.waveforms = True

    config_types = []

    sub_dirs = ['gf_store']
    empty_files = []

    dataset_config = grond.DatasetConfig(events_path='events.txt')
    target_groups = []

    if options.waveforms:
        config_types.append('seismic waveforms')

        sub_dirs += ['data']
        dataset_config.waveform_paths = ['data']

        empty_files += ['stations.xml']
        dataset_config.stations_path = 'stations.txt'

        target_groups.append(
            grond.WaveformTargetGroup(
                normalisation_family='time_domain',
                path='all',
                distance_min=10*km,
                distance_max=1000*km,
                channels=['Z', 'R', 'T'],
                interpolation='multilinear',
                store_id='gf_store',
                misfit_config=grond.WaveformMisfitConfig(
                    fmin=0.01,
                    fmax=0.1)))

        pi2 = math.pi/2
        problem_config = grond.CMTProblemConfig(
            name_template='cmt_%(event_name)s',
            distance_min=2.*km,
            mt_type='deviatoric',
            ranges=dict(
                time=Range(0, 10.0, relative='add'),
                north_shift=Range(-16*km, 16*km),
                east_shift=Range(-16*km, 16*km),
                depth=Range(1*km, 11*km),
                magnitude=Range(4.0, 6.0),
                rmnn=Range(-pi2, pi2),
                rmee=Range(-pi2, pi2),
                rmdd=Range(-pi2, pi2),
                rmne=Range(-1.0, 1.0),
                rmnd=Range(-1.0, 1.0),
                rmed=Range(-1.0, 1.0),
                duration=Range(1.0, 15.0))
            )

    if options.insar:
        config_types.append('InSAR')

        sub_dirs += ['scenes', 'gnss']

        dataset_config.kite_scene_paths = ['scenes']

        target_groups.append(
            grond.SatelliteTargetGroup(
                normalisation_family='insar_target',
                path='all',
                interpolation='multilinear',
                store_id='gf_store',
                kite_scenes=['*all'],
                misfit_config=grond.SatelliteMisfitConfig(
                    use_weight_focal=False,
                    optimize_orbital_ramp=True,
                    ranges={
                        'offset': '-0.5 .. 0.5',
                        'ramp_north': '-1e-4 .. 1e-4',
                        'ramp_east': '-1e-4 .. 1e-4'
                        }
                    )))

        problem_config = grond.RectangularProblemConfig(
            name_template='rect_source',
            ranges=dict(
                north_shift=Range(-20*km, 20*km),
                east_shift=Range(-20*km, 20*km),
                depth=Range(0*km, 10*km),
                length=Range(20*km, 40*km),
                width=Range(5*km, 12*km),
                dip=Range(20, 70),
                strike=Range(0, 180),
                rake=Range(0, 90),
                slip=Range(1, 3))
            )

    engine_config = grond.EngineConfig(
        gf_store_superdirs=['.'])

    optimizer_config = grond.HighScoreOptimizerConfig()

    config = grond.Config(
        rundir_template=op.join('rundir', '${problem_name}.grun'),
        dataset_config=dataset_config,
        target_groups=target_groups,
        problem_config=problem_config,
        optimizer_config=optimizer_config,
        engine_config=engine_config)

    events = '''name = 2011-myanmar
time = 2011-03-24 13:55:12.010
latitude = 20.687
longitude = 99.822
magnitude = 6.9
moment = 1.9228e+19
depth = 8000
region = Myanmar
--------------------------------------------'''

    if project_dir is not None:
        logger.info('Creating empty project for %s in folder %s'
                    % (' and '.join(config_types), args[0]))

        def p(fn):
            return op.join(project_dir, fn)

        os.mkdir(op.abspath(project_dir))
        for d in sub_dirs:
            os.mkdir(p(d))

        with open(p('config.yml'), 'w') as cf:
            cf.write(str(config))
        with open(p('events.txt'), 'w') as ef:
            ef.write(events)

        for fn in empty_files:
            open(p(fn), 'w').close()
    else:
        print(config)


def command_events(args):
    def setup(parser):
        pass

    parser, options, args = cl_parse('events', args, setup)
    if len(args) != 1:
        help_and_die(parser, 'missing arguments')

    config_path = args[0]
    config = grond.read_config(config_path)

    print('Available Events:')
    for event_name in grond.get_event_names(config):
        print('* %s' % event_name)


def command_check(args):
    def setup(parser):
        parser.add_option(
            '--target-ids', dest='target_string_ids', metavar='TARGET_IDS',
            help='process only selected targets. TARGET_IDS is a '
                 'comma-separated list of target IDs. Target IDs have the '
                 'form SUPERGROUP.GROUP.NETWORK.STATION.LOCATION.CHANNEL.')

        parser.add_option(
            '--plot', dest='show_plot', action='store_true',
            help='plot sample synthetics and data.')

        parser.add_option(
            '--waveforms', dest='show_waveforms', action='store_true',
            help='show raw, restituted, projected, and processed waveforms')

        parser.add_option(
            '--nrandom', dest='n_random_synthetics', metavar='N', type='int',
            default=10,
            help='set number of random synthetics to forward model (default: '
                 '10). If set to zero, create synthetics for the reference '
                 'solution.')

    parser, options, args = cl_parse('check', args, setup)
    if len(args) < 2:
        help_and_die(parser, 'missing arguments')

    config_path = args[0]
    event_names = args[1:]

    config = grond.read_config(config_path)

    target_string_ids = None
    if options.target_string_ids:
        target_string_ids = options.target_string_ids.split(',')

    grond.check(
        config,
        event_names=event_names,
        target_string_ids=target_string_ids,
        show_plot=options.show_plot,
        show_waveforms=options.show_waveforms,
        n_random_synthetics=options.n_random_synthetics)


def command_go(args):
    def setup(parser):
        parser.add_option(
            '--force', dest='force', action='store_true',
            help='overwrite existing run directory')
        parser.add_option(
            '--status', dest='status', default='state',
            help='status output selection (choices: state, matrix)')
        parser.add_option(
            '--parallel', dest='nparallel', type='int', default=1,
            help='set number of events to process in parallel')

    parser, options, args = cl_parse('go', args, setup)

    if len(args) == 1:
        config_path = args[0]
        config = grond.read_config(config_path)

        for event_name in grond.get_event_names(config):
            print(event_name)

        help_and_die(parser, 'missing arguments')

    if len(args) < 2:
        help_and_die(parser, 'missing arguments')

    config_path = args[0]
    event_names = args[1:]

    config = grond.read_config(config_path)
    if options.status == 'quiet':
        status = ()
    else:
        status = tuple(options.status.split(','))

    grond.go(
        config,
        event_names=event_names,
        force=options.force,
        status=status,
        nparallel=options.nparallel)


def command_forward(args):
    def setup(parser):
        pass

    parser, options, args = cl_parse('forward', args, setup)
    if len(args) < 1:
        help_and_die(parser, 'missing arguments')

    event_names = args[1:]

    if not event_names:
        help_and_die(parser, 'no event names given')

    run_path = args[0]
    grond.forward(
        run_path,
        event_names=event_names)


def command_harvest(args):
    def setup(parser):
        parser.add_option(
            '--force', dest='force', action='store_true',
            help='overwrite existing harvest directory')
        parser.add_option(
            '--neach', dest='neach', type=int, default=10,
            help='take NEACH best samples from each chain (default: 10)')
        parser.add_option(
            '--weed', dest='weed', type=int, default=0,
            help='weed out bootstrap samples with bad global performance. '
                 '0: no weeding (default), '
                 '1: only bootstrap chains where all NEACH best samples '
                 'global misfit is less than the global average misfit of all '
                 'NEACH best in all chains plus one standard deviation are '
                 'included in the harvest ensemble, '
                 '2: same as 1 but additionally individual samples are '
                 'removed if their global misfit is greater than the global '
                 'average misfit of all NEACH best in all chains, '
                 '3: harvesting is done on the global chain only, bootstrap '
                 'chains are excluded')

    parser, options, args = cl_parse('harvest', args, setup)
    if len(args) != 1:
        help_and_die(parser, 'no rundir')

    run_path, = args
    grond.harvest(
        run_path,
        force=options.force,
        nbest=options.neach,
        weed=options.weed)


def command_plot(args):

    def setup(parser):
        parser.add_option(
            '--save', dest='save', action='store_true', default=False,
            help='save figures to files')

        parser.add_option(
            '--format', '--formats', dest='formats', default='pdf',
            help='comma-separated list of ouptut formats (default: pdf)')

        parser.add_option(
            '--dpi', '--dpi', dest='dpi', type=float, default=120.,
            help='DPI setting for raster formats (default=120)')

    # plotnames_avail = plot.available_plotnames()
    # made explicit to avoid import of pyplot before backend can be chosen
    plotnames_avail = [
        'bootstrap',
        'sequence',
        'contributions',
        'jointpar',
        'histogram',
        'hudson',
        'fits',
        'fits_statics',
        'fits_ensemble',
        'solution',
        'location']

    details = '''Available <plotnames> are: %s, or "all". Multiple plots are
selected by specifying a comma-separated list.''' % (
        ', '.join('"%s"' % x for x in plotnames_avail))

    parser, options, args = cl_parse('plot', args, setup, details)

    if len(args) != 2:
        help_and_die(parser, 'two arguments required')

    if args[0] == 'all':
        plotnames = plotnames_avail
    else:
        plotnames = args[0].split(',')

    formats = options.formats.split(',')
    dirname = args[1]

    if options.save:
        import matplotlib
        matplotlib.use('Agg')

    from grond import plot

    try:
        plot.plot_result(
            dirname, plotnames,
            save=options.save, formats=formats, dpi=options.dpi)

    except grond.GrondError as e:
        die(str(e))


def command_movie(args):

    def setup(parser):
        pass

    parser, options, args = cl_parse('movie', args, setup)

    if len(args) != 4:
        help_and_die(parser, 'four arguments required')

    run_path, xpar_name, ypar_name, movie_filename_template = args

    from grond import plot

    movie_filename = movie_filename_template % {
        'xpar': xpar_name,
        'ypar': ypar_name}

    try:
        plot.make_movie(run_path, xpar_name, ypar_name, movie_filename)

    except grond.GrondError as e:
        die(str(e))


def command_export(args):

    def setup(parser):
        parser.add_option(
            '--type', dest='type', metavar='TYPE',
            choices=('event', 'event-yaml', 'source', 'vector'),
            help='select type of objects to be exported. Choices: '
                 '"event" (default), "event-yaml", "source", "vector".')

        parser.add_option(
            '--parameters', dest='parameters', metavar='PLIST',
            help='select parameters to be exported. PLIST is a '
                 'comma-separated list where each entry has the form '
                 '"<parameter>[.<measure>]". Available measures: "best", '
                 '"mean", "std", "minimum", "percentile16", "median", '
                 '"percentile84", "maximum".')

        parser.add_option(
            '--output', dest='filename', metavar='FILE',
            help='write output to FILE')

    parser, options, args = cl_parse('export', args, setup)
    if len(args) < 2:
        help_and_die(parser, 'arguments required')

    what = args[0]

    dirnames = args[1:]

    what_choices = ('best', 'mean', 'ensemble', 'stats')

    if what not in what_choices:
        help_and_die(
            parser,
            'invalid choice: %s (choose from %s)' % (
                repr(what), ', '.join(repr(x) for x in what_choices)))

    if options.parameters:
        pnames = options.parameters.split(',')
    else:
        pnames = None

    try:
        grond.export(
            what,
            dirnames,
            filename=options.filename,
            type=options.type,
            pnames=pnames)

    except grond.GrondError as e:
        die(str(e))


def command_report(args):

    import grond.report

    def setup(parser):
        pass

    parser, options, args = cl_parse('report', args, setup)
    if len(args) < 1:
        help_and_die(parser, 'arguments required')

    if all(op.isdir(rundir) for rundir in args):
        rundirs = args
        try:
            for rundir in rundirs:
                grond.report.report(rundir)

        except grond.GrondError as e:
            die(str(e))

    else:
        if len(args) < 2:
            help_and_die(parser, 'arguments required')

        config_path = args[0]
        event_names = args[1:]

        try:
            config = grond.read_config(config_path)
            for event_name in event_names:
                grond.report.report(config_and_event_name=(config, event_name))

        except grond.GrondError as e:
            die(str(e))


def command_report_index(args):

    import grond.report

    def setup(parser):
        pass

    parser, options, args = cl_parse('report-index', args, setup)
    if len(args) < 1:
        help_and_die(parser, 'arguments required')

    report_base_path = args[0]

    try:
        grond.report.report_index(report_base_path)

    except grond.GrondError as e:
        die(str(e))


def command_qc_polarization(args):

    def setup(parser):
        parser.add_option(
            '--time-factor-pre', dest='time_factor_pre', type=float,
            metavar='NUMBER',
            default=0.5,
            help='set duration to extract before synthetic P phase arrival, '
                 'relative to 1/fmin. fmin is taken from the selected target '
                 'group in the config file (default=%default)')
        parser.add_option(
            '--time-factor-post', dest='time_factor_post', type=float,
            metavar='NUMBER',
            default=0.5,
            help='set duration to extract after synthetic P phase arrival, '
                 'relative to 1/fmin. fmin is taken from the selected target '
                 'group in the config file (default=%default)')
        parser.add_option(
            '--distance-min', dest='distance_min', type=float,
            metavar='NUMBER',
            help='minimum event-station distance [m]')
        parser.add_option(
            '--distance-max', dest='distance_max', type=float,
            metavar='NUMBER',
            help='maximum event-station distance [m]')
        parser.add_option(
            '--depth-min', dest='depth_min', type=float,
            metavar='NUMBER',
            help='minimum station depth [m]')
        parser.add_option(
            '--depth-max', dest='depth_max', type=float,
            metavar='NUMBER',
            help='maximum station depth [m]')
        parser.add_option(
            '--picks', dest='picks_filename',
            metavar='FILENAME',
            help='add file with P picks in Snuffler marker format')
        parser.add_option(
            '--save', dest='output_filename',
            metavar='FILENAME.FORMAT',
            help='save output to file FILENAME.FORMAT')
        parser.add_option(
            '--dpi', dest='output_dpi', type=float, default=120.,
            metavar='NUMBER',
            help='DPI setting for raster formats (default=120)')

    parser, options, args = cl_parse('qc-polarization', args, setup)
    if len(args) != 3:
        help_and_die(parser, 'missing arguments')

    if options.output_filename:
        import matplotlib
        matplotlib.use('Agg')

    import grond.qc

    config_path, event_name, target_config_name = args

    config = grond.read_config(config_path)

    ds = config.get_dataset(event_name)

    engine = config.engine_config.get_engine()

    nsl_to_time = None
    if options.picks_filename:
        markers = marker.load_markers(options.picks_filename)
        marker.associate_phases_to_events(markers)

        nsl_to_time = {}
        for m in markers:
            if isinstance(m, marker.PhaseMarker):
                ev = m.get_event()
                if ev is not None and ev.name == event_name:
                    nsl_to_time[m.one_nslc()[:3]] = m.tmin

        if not nsl_to_time:
            help_and_die(
                parser,
                'no markers associated with event "%s" found in file "%s"' % (
                    event_name, options.picks_filename))

    target_config_names_avail = []
    for target_config in config.target_configs:
        name = '%s.%s' % (target_config.super_group, target_config.group or '')
        if name == target_config_name:
            imc = target_config.inner_misfit_config
            fmin = imc.fmin
            fmax = imc.fmax
            ffactor = imc.ffactor

            store = engine.get_store(target_config.store_id)
            timing = '{cake:P|cake:p|cake:P\\|cake:p\\}'

            grond.qc.polarization(
                ds, store, timing, fmin=fmin, fmax=fmax, ffactor=ffactor,
                time_factor_pre=options.time_factor_pre,
                time_factor_post=options.time_factor_post,
                distance_min=options.distance_min,
                distance_max=options.distance_max,
                depth_min=options.depth_min,
                depth_max=options.depth_max,
                nsl_to_time=nsl_to_time,
                output_filename=options.output_filename,
                output_dpi=options.output_dpi)

            return

        target_config_names_avail.append(name)

        die('no target_config named "%s" found. Available: %s' % (
            target_config_name, ', '.join(target_config_names_avail)))


if __name__ == '__main__':
    main()
