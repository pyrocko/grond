import os.path as op
import glob
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
    plot_names = List.T(
        String.T(),
        default=['solution', 'fits', 'jointpar', 'hudson', 'sequence',
                 'fits_ensemble', 'contributions', 'bootstrap'])


class ReportPlot(Object):
    name = String.T()
    file_names = List.T(String.T())


class ReportPlots(Object):
    plots = List.T(ReportPlot.T())


def report(rundir, report_config=None):

    if report_config is None:
        report_config = ReportConfig()

    config = guts.load(
        filename=op.join(rundir, 'config.yaml'))

    config.set_basepath(rundir)

    problem = load_problem_info(rundir)

    reportdir = expand_template(
        report_config.reportdir_template,
        dict(
            event_name=problem.base_source.name,
            problem_name=problem.name))

    util.ensuredir(reportdir)
    core.export('stats', [rundir], filename=op.join(reportdir, 'stats.yaml'))

    plots_dir_out = op.join(reportdir, 'plots')
    util.ensuredir(plots_dir_out)
    rps = []
    for plot_name in report_config.plot_names:
        plot_path_pat = op.join(rundir, 'plots', plot_name + '-*.png')
        plot_paths = sorted(glob.glob(plot_path_pat))
        plot_basenames = []
        for plot_path in plot_paths:
            plot_basename = op.basename(plot_path)
            plot_path_out = op.join(plots_dir_out, plot_basename)
            shutil.copy(plot_path, plot_path_out)
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
            stats_path = op.join(dirpath, 'stats.yaml')
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
