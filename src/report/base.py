import os.path as op
import shutil
import os

from pyrocko import guts, util
from pyrocko.guts import Object, List, String

from grond.meta import Path, expand_template

from grond import core
from grond.problems.base import load_problem_info

guts_prefix = 'grond'


class ReportConfig(Object):
    reportdir_template = Path.T(
        default='reports/${event_name}/${problem_name}')


class ReportPlot(Object):
    name = String.T()
    file_names = List.T(String.T())


class ReportPlots(Object):
    plots = List.T(ReportPlot.T())


def report(
        rundir=None, config_and_event_name=None, report_config=None,
        skip_plotting=False):

    if report_config is None:
        report_config = ReportConfig()

    from grond.config import read_config
    from grond.core import check
    from grond.plot import plot_result, available_plotnames

    if config_and_event_name is None:
        config = read_config(op.join(rundir, 'config.yaml'))
        problem = load_problem_info(rundir)
        event_name = problem.base_source.name

    else:
        config, event_name = config_and_event_name
        ds = config.get_dataset(event_name)
        event = ds.get_event()
        problem = config.get_problem(event)

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

    fns = {}
    if 'target_check' in report_config.plot_names:
        fns.update(check(
            config,
            event_names=[event_name],
            save_plot_path=op.join(plots_dir_out)))

    if rundir:
        core.export(
            'stats', [rundir], filename=op.join(reportdir, 'stats.yaml'))

        avail = available_plotnames()
        fns.update(plot_result(
            rundir,
            [name for name in report_config.plot_names if name in avail],
            save_path=plots_dir_out, formats=('png',)))

    rps = []
    for plot_name in sorted(
            fns.keys(),
            key=lambda x: report_config.plot_names.index(x)):

        plot_paths = fns.get(plot_name, [])
        plot_basenames = []
        for plot_path in plot_paths:
            plot_basename = op.basename(plot_path)
            plot_basenames.append(plot_basename)

        rp = ReportPlot(
            name=plot_name,
            file_names=plot_basenames)

        rps.append(rp)

    guts.dump_all(rps, filename=op.join(reportdir, 'plots.yaml'))


def iter_report_dirs(report_base_path):
    for path, dirnames, filenames in os.walk(report_base_path):
        for dirname in dirnames:
            dirpath = op.join(path, dirname)
            stats_path = op.join(dirpath, 'problem.yaml')
            if op.exists(stats_path):
                yield dirpath


class ReportEntry(Object):
    path = String.T()


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


def report_index(report_base_path):
    reports = []
    for report_path in iter_report_dirs(report_base_path):
        report_relpath = op.relpath(report_path, report_base_path)
        reports.append(ReportEntry(
            path=report_relpath))

    guts.dump_all(
        reports,
        filename=op.join(report_base_path, 'report_list.yaml'))

    app_dir = op.join(op.split(__file__)[0], 'app')
    copytree(app_dir, report_base_path)


__all__ = '''
    report
    report_index
    ReportConfig
    ReportPlot
    ReportEntry
'''.split()
