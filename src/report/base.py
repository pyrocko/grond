import logging
import os.path as op
import shutil
import os

from http.server import HTTPServer, SimpleHTTPRequestHandler

from pyrocko import guts, util
from pyrocko.model import Event
from pyrocko.guts import Object, String

from grond.meta import Path, expand_template

from grond import core, environment

guts_prefix = 'grond'
logger = logging.getLogger('grond.report')


class ReportIndexEntry(Object):
    path = String.T()
    problem_name = String.T()
    event_reference = Event.T(optional=True)
    event_best = Event.T(optional=True)


class ReportConfig(Object):
    reports_base_path = Path.T(default='reports')
    report_sub_path = String.T(
        default='${event_name}/${problem_name}')


def report(env, report_config=None, update_without_plotting=False):
    if report_config is None:
        report_config = ReportConfig()

    event_name = env.get_current_event_name()
    problem = env.get_problem()
    logger.info('Creating report for event %s...' % event_name)

    report_path = expand_template(
        op.join(
            report_config.reports_base_path,
            report_config.report_sub_path),
        dict(
            event_name=event_name,
            problem_name=problem.name))

    if op.exists(report_path) and not update_without_plotting:
        shutil.rmtree(report_path)

    problem.dump_problem_info(report_path)

    util.ensuredir(report_path)
    plots_dir_out = op.join(report_path, 'plots')
    util.ensuredir(plots_dir_out)

    event = env.get_dataset().get_event()
    guts.dump(event, filename=op.join(report_path, 'event.reference.yaml'))

    try:
        rundir_path = env.get_rundir_path()

        core.export(
            'stats', [rundir_path],
            filename=op.join(report_path, 'stats.yaml'))

        core.export(
            'best', [rundir_path],
            filename=op.join(report_path, 'event.solution.best.yaml'),
            type='event-yaml')

        core.export(
            'mean', [rundir_path],
            filename=op.join(report_path, 'event.solution.mean.yaml'),
            type='event-yaml')

        core.export(
            'ensemble', [rundir_path],
            filename=op.join(report_path, 'event.solution.ensemble.yaml'),
            type='event-yaml')

    except environment.NoRundirAvailable:
        pass

    if not update_without_plotting:
        from grond import plot
        plot.make_plots(env, plots_path=op.join(report_path, 'plots'))

    rie = ReportIndexEntry(path='.', problem_name=problem.name)

    fn = op.join(report_path, 'event.solution.best.yaml')
    if op.exists(fn):
        rie.event_best = guts.load(filename=fn)

    fn = op.join(report_path, 'event.reference.yaml')
    if op.exists(fn):
        rie.event_reference = guts.load(filename=fn)

    fn = op.join(report_path, 'index.yaml')
    guts.dump(rie, filename=fn)

    report_index(report_config)


def report_index(report_config=None):
    if report_config is None:
        report_config = ReportConfig()

    reports_base_path = report_config.reports_base_path
    reports = []
    for report_path in iter_report_dirs(reports_base_path):
        logger.info('Indexing %s...' % report_path)

        fn = op.join(report_path, 'index.yaml')
        rie = guts.load(filename=fn)
        report_relpath = op.relpath(report_path, reports_base_path)
        rie.path = report_relpath
        reports.append(rie)

    guts.dump_all(
        reports,
        filename=op.join(reports_base_path, 'report_list.yaml'))

    app_dir = op.join(op.split(__file__)[0], 'app')
    copytree(app_dir, reports_base_path)
    logger.info('Created report in %s/index.html' % reports_base_path)


def iter_report_dirs(reports_base_path):
    for path, dirnames, filenames in os.walk(reports_base_path):
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


class ReportHandler(SimpleHTTPRequestHandler):

    def log_error(self, fmt, *args):
        logger.error(fmt % args)

    def log_message(self, fmt, *args):
        logger.debug(fmt % args)


def serve_report(host=('127.0.0.1', 8383), report_config=None):
    if report_config is None:
        report_config = ReportConfig()

    logger.info('Starting report webserver at http://%s:%d...' % host)

    os.chdir(report_config.reports_base_path)
    httpd = HTTPServer(host, ReportHandler)
    httpd.serve_forever()


__all__ = '''
    report
    report_index
    ReportConfig
    ReportIndexEntry
    serve_report
'''.split()
