from __future__ import print_function
import logging
import os.path as op
import shutil
import os
import tarfile
import threading
import signal
import time

from http.server import HTTPServer, SimpleHTTPRequestHandler

from pyrocko import guts, util
from pyrocko.model import Event
from pyrocko.guts import Object, String, Unicode, Bool

from grond.meta import HasPaths, Path, expand_template, GrondError

from grond import core, environment
from grond.problems import ProblemInfoNotAvailable, ProblemDataNotAvailable
from grond.version import __version__
from grond import info
from grond.plot import PlotConfigCollection, get_all_plot_classes
from grond.run_info import RunInfo

guts_prefix = 'grond'
logger = logging.getLogger('grond.report')


class ReportIndexEntry(Object):
    path = String.T()
    problem_name = String.T()
    event_reference = Event.T(optional=True)
    event_best = Event.T(optional=True)
    grond_version = String.T(optional=True)
    run_info = RunInfo.T(optional=True)


class ReportConfig(HasPaths):
    report_base_path = Path.T(default='report')
    entries_sub_path = String.T(
        default='${event_name}/${problem_name}')
    title = Unicode.T(
        default=u'Grond Report',
        help='Title shown on report overview page.')
    description = Unicode.T(
        default=u'This interactive document aggregates earthquake source '
                u'inversion results from optimisations performed with Grond.',
        help='Description shown on report overview page.')
    plot_config_collection = PlotConfigCollection.T(
        help='Configurations for plots to be included in the report.')
    make_archive = Bool.T(
        default=True,
        help='Set to `false` to prevent creation of compressed archive.')


class ReportInfo(Object):
    title = Unicode.T(optional=True)
    description = Unicode.T(optional=True)
    version_info = info.VersionInfo.T()
    have_archive = Bool.T(optional=True)


def read_config(path):
    get_all_plot_classes()  # make sure all plot modules are imported
    try:
        config = guts.load(filename=path)
    except OSError:
        raise GrondError(
            'Cannot read Grond report configuration file: %s' % path)

    if not isinstance(config, ReportConfig):
        raise GrondError(
            'Invalid Grond report configuration in file "%s".' % path)

    config.set_basepath(op.dirname(path) or '.')
    return config


def write_config(config, path):
    try:
        basepath = config.get_basepath()
        dirname = op.dirname(path) or '.'
        config.change_basepath(dirname)
        guts.dump(
            config,
            filename=path,
            header='Grond report configuration file, version %s' % __version__)

        config.change_basepath(basepath)

    except OSError:
        raise GrondError(
            'Cannot write Grond report configuration file: %s' % path)


def iter_report_entry_dirs(report_base_path):
    for path, dirnames, filenames in os.walk(report_base_path):
        for dirname in dirnames:
            dirpath = op.join(path, dirname)
            stats_path = op.join(dirpath, 'problem.yaml')
            if op.exists(stats_path):
                yield dirpath


def copytree(src, dst):
    names = os.listdir(src)
    if not op.exists(dst):
        os.makedirs(dst)

    for name in names:
        srcname = op.join(src, name)
        dstname = op.join(dst, name)
        if op.isdir(srcname):
            copytree(srcname, dstname)
        else:
            shutil.copy(srcname, dstname)


def report(env, report_config=None, update_without_plotting=False,
           make_index=True, make_archive=True, nthreads=0):

    if report_config is None:
        report_config = ReportConfig()
        report_config.set_basepath('.')

    event_name = env.get_current_event_name()
    problem = env.get_problem()
    logger.info('Creating report entry for run "%s"...' % problem.name)

    optimiser = env.get_optimiser()
    optimiser.set_nthreads(nthreads)

    fp = report_config.expand_path
    entry_path = expand_template(
        op.join(
            fp(report_config.report_base_path),
            report_config.entries_sub_path),
        dict(
            event_name=event_name,
            problem_name=problem.name))

    if op.exists(entry_path) and not update_without_plotting:
        shutil.rmtree(entry_path)

    try:
        problem.dump_problem_info(entry_path)

        guts.dump(env.get_config(),
                  filename=op.join(entry_path, 'config.yaml'),
                  header=True)

        util.ensuredir(entry_path)
        plots_dir_out = op.join(entry_path, 'plots')
        util.ensuredir(plots_dir_out)

        event = env.get_dataset().get_event()
        guts.dump(event, filename=op.join(entry_path, 'event.reference.yaml'))

        try:
            rundir_path = env.get_rundir_path()

            core.export(
                'stats', [rundir_path],
                filename=op.join(entry_path, 'stats.yaml'))

            core.export(
                'best', [rundir_path],
                filename=op.join(entry_path, 'event.solution.best.yaml'),
                type='event-yaml')

            core.export(
                'mean', [rundir_path],
                filename=op.join(entry_path, 'event.solution.mean.yaml'),
                type='event-yaml')

            core.export(
                'ensemble', [rundir_path],
                filename=op.join(entry_path, 'event.solution.ensemble.yaml'),
                type='event-yaml')

        except (environment.NoRundirAvailable, ProblemInfoNotAvailable,
                ProblemDataNotAvailable):

            pass

        if not update_without_plotting:
            from grond import plot
            pcc = report_config.plot_config_collection.get_weeded(env)
            plot.make_plots(
                env,
                plots_path=op.join(entry_path, 'plots'),
                plot_config_collection=pcc)

        try:
            run_info = env.get_run_info()
        except environment.NoRundirAvailable:
            run_info = None

        rie = ReportIndexEntry(
            path='.',
            problem_name=problem.name,
            grond_version=problem.grond_version,
            run_info=run_info)

        fn = op.join(entry_path, 'event.solution.best.yaml')
        if op.exists(fn):
            rie.event_best = guts.load(filename=fn)

        fn = op.join(entry_path, 'event.reference.yaml')
        if op.exists(fn):
            rie.event_reference = guts.load(filename=fn)

        fn = op.join(entry_path, 'index.yaml')
        guts.dump(rie, filename=fn)

    except Exception as e:
        logger.warn(
            'Failed to create report entry, removing incomplete subdirectory: '
            '%s' % entry_path)
        raise e

        if op.exists(entry_path):
            shutil.rmtree(entry_path)

    logger.info('Done creating report entry for run "%s".' % problem.name)

    if make_index:
        report_index(report_config)

    if make_archive:
        report_archive(report_config)


def report_index(report_config=None):
    if report_config is None:
        report_config = ReportConfig()

    report_base_path = report_config.report_base_path
    entries = []
    for entry_path in iter_report_entry_dirs(report_base_path):

        fn = op.join(entry_path, 'index.yaml')
        if not os.path.exists(fn):
            logger.warn('Skipping indexing of incomplete report entry: %s'
                        % entry_path)

            continue

        logger.info('Indexing %s...' % entry_path)

        rie = guts.load(filename=fn)
        report_relpath = op.relpath(entry_path, report_base_path)
        rie.path = report_relpath
        entries.append(rie)

    guts.dump_all(
        entries,
        filename=op.join(report_base_path, 'report_list.yaml'))

    guts.dump(
        ReportInfo(
            title=report_config.title,
            description=report_config.description,
            version_info=info.version_info(),
            have_archive=report_config.make_archive),
        filename=op.join(report_base_path, 'info.yaml'))

    app_dir = op.join(op.split(__file__)[0], 'app')
    copytree(app_dir, report_base_path)
    logger.info('Created report in %s/index.html' % report_base_path)


def report_archive(report_config):
    if report_config is None:
        report_config = ReportConfig()

    if not report_config.make_archive:
        return

    report_base_path = report_config.report_base_path

    logger.info('Generating report\'s archive...')
    with tarfile.open(op.join(report_base_path, 'grond-report.tar.gz'),
                      mode='w:gz') as tar:
        tar.add(report_base_path, arcname='grond-report')


def serve_ip(host):
    if host == 'localhost':
        ip = '127.0.0.1'
    elif host == 'default':
        import socket
        ip = [
            (s.connect(('4.4.4.4', 80)), s.getsockname()[0], s.close())
            for s in [socket.socket(socket.AF_INET, socket.SOCK_DGRAM)]][0][1]
    elif host == '*':
        ip = ''
    else:
        ip = host

    return ip


class ReportHandler(SimpleHTTPRequestHandler):

    def _log_error(self, fmt, *args):
        logger.error(fmt % args)

    def _log_message(self, fmt, *args):
        logger.debug(fmt % args)

    def end_headers(self):
        self.send_header('Cache-Control', 'no-cache')
        SimpleHTTPRequestHandler.end_headers(self)


ReportHandler.extensions_map.update({
    '.yaml': 'application/x-yaml',
    '.yml': 'application/x-yaml'})


g_terminate = False


def serve_report(
        addr=('127.0.0.1', 8383),
        report_config=None,
        fixed_port=False,
        open=False):

    if report_config is None:
        report_config = ReportConfig()

    path = report_config.expand_path(report_config.report_base_path)
    os.chdir(path)

    host, port = addr
    if fixed_port:
        ports = [port]
    else:
        ports = range(port, port+20)

    httpd = None
    for port in ports:
        try:
            httpd = HTTPServer((host, port), ReportHandler)
            break
        except OSError as e:
            logger.warn(str(e))

    if httpd:
        logger.info(
            'Starting report web service at http://%s:%d' % (host, port))

        thread = threading.Thread(None, httpd.serve_forever)
        thread.start()

        if open:
            import webbrowser
            if open:
                webbrowser.open('http://%s:%d' % (host, port))

        def handler(signum, frame):
            global g_terminate
            g_terminate = True

        signal.signal(signal.SIGINT, handler)
        signal.signal(signal.SIGTERM, handler)

        while not g_terminate:
            time.sleep(0.1)

        signal.signal(signal.SIGINT, signal.SIG_DFL)
        signal.signal(signal.SIGTERM, signal.SIG_DFL)

        logger.info('Stopping report web service...')

        httpd.shutdown()
        thread.join()

        logger.info('... done')

    else:
        logger.error('Failed to start web service.')


__all__ = '''
    report
    report_index
    report_archive
    ReportConfig
    ReportIndexEntry
    ReportInfo
    serve_ip
    serve_report
    read_config
    write_config
'''.split()
