import logging
import os.path as op
import shutil
import os

from pyrocko import guts, util
from pyrocko.guts import Object, String

from grond.meta import Path, expand_template

from grond import core

guts_prefix = 'grond'
logger = logging.getLogger('grond.report')


class ReportEntry(Object):
    path = String.T()


class ReportConfig(Object):
    reportdir_template = Path.T(
        default='reports/${event_name}/${problem_name}')


def report(env, report_config=None):
    if report_config is None:
        report_config = ReportConfig()

    event_name = env.get_event_name()
    problem = env.get_problem()
    logger.info('Creating report for event %s...' % event_name)

    reportdir = expand_template(
        report_config.reportdir_template,
        dict(
            event_name=event_name,
            problem_name=problem.name))

    if op.exists(reportdir):
        shutil.rmtree(reportdir)

    problem.dump_problem_info(reportdir)

    util.ensuredir(reportdir)
    plots_dir_out = op.join(reportdir, 'plots')
    util.ensuredir(plots_dir_out)

    rundir_path = env.get_rundir_path()

    core.export(
        'stats', [rundir_path], filename=op.join(reportdir, 'stats.yaml'))

    from grond import plot
    plot.make_plots(env, plots_path=op.join(reportdir, 'plots'))

    report_index('reports/')


def report_index(report_base_path):
    reports = []
    for report_path in iter_report_dirs(report_base_path):
        logger.info('Indexing %s...' % report_path)
        report_relpath = op.relpath(report_path, report_base_path)
        reports.append(ReportEntry(
            path=report_relpath))

    guts.dump_all(
        reports,
        filename=op.join(report_base_path, 'report_list.yaml'))

    app_dir = op.join(op.split(__file__)[0], 'app')
    copytree(app_dir, report_base_path)


def iter_report_dirs(report_base_path):
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


__all__ = '''
    report
    report_index
    ReportConfig
    ReportEntry
'''.split()
