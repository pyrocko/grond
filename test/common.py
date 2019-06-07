from __future__ import absolute_import

import logging
import os
import shutil
import sys
from io import StringIO

from pyrocko import util
from grond.apps.grond import main
from grond import Environment
from grond.meta import expand_template

url = 'http://data.pyrocko.org/testing/grond/'

logger = logging.getLogger('grond.test')
op = os.path


def test_data_path(fn):
    return op.abspath(os.path.join(os.path.split(__file__)[0], 'data', fn))


def get_test_data(path):
    fn = test_data_path(path)
    if not op.exists(fn):
        util.ensuredirs(fn)
        if path.endswith('/'):
            util.download_dir(op.join(url, path), fn)
        else:
            util.download_file(op.join(url, path), fn)

    return fn


def put_test_data(path, dest=None):
    if dest is None:
        dest = path
    fn = get_test_data(path)
    util.ensuredirs(dest)
    shutil.copytree(fn, dest)


def link_test_data(path, dest=None):
    if dest is None:
        dest = path
    fn = get_test_data(path)
    util.ensuredirs(dest)
    os.symlink(fn, dest.rstrip(os.sep))


class GrondExit(Exception):
    def __init__(self, res):
        Exception.__init__(self, str(res))
        self.result = res


class Capture(object):
    def __init__(self, tee=False):
        self.file = StringIO()
        self.tee = tee

    def __enter__(self):
        self.orig_stdout = sys.stdout
        self.orig_exit = sys.exit
        sys.stdout = self

        def my_exit(res):
            raise GrondExit(res)

        sys.exit = my_exit

    def __exit__(self, *args):
        sys.stdout = self.orig_stdout
        sys.exit = self.orig_exit

    def write(self, data):
        self.file.write(data)
        if self.tee:
            self.orig_stdout.write(data)

    def writelines(self, lines):
        for l in lines:
            self.write(l)

    def flush(self):
        self.file.flush()

    def isatty(self):
        return False

    def getvalue(self):
        return self.file.getvalue()


def grond(*args, **kwargs):
    # tee = True
    tee = kwargs.get('tee', False)
    logger.info('Calling: grond %s' % ' '.join(args))
    cap = Capture(tee=tee)
    with cap:
        main(['grond'] + list(args))

    return cap.getvalue()


def assert_grond_usage(*args):
    res = None
    try:
        grond(*args)
    except GrondExit as e:
        res = e.result

    assert res.startswith('Usage')


def get_playground_dir():
    playground_dir = 'test_playground'
    util.ensuredir(playground_dir)
    return playground_dir


def get_rundir_paths(config_path, event_names):
    env = Environment([config_path] + event_names)
    conf = env.get_config()

    rundir_paths = []
    for event_name in event_names:
        env.set_current_event_name(event_name)
        problem_name = env.get_problem().name
        rundir_paths.append(expand_template(
            conf.rundir_template,
            dict(problem_name=problem_name)))

    return rundir_paths


class chdir(object):

    def __init__(self, path):
        self._path = path

    def __enter__(self):
        self._oldwd = os.getcwd()
        os.chdir(self._path)

    def __exit__(self, *args):
        os.chdir(self._oldwd)


def run_in_project(
        main,
        project_dir_source,
        project_dir,
        event_name,
        config_path):

    playground_dir = get_playground_dir()
    with chdir(playground_dir):

        if os.path.exists(project_dir):
            shutil.rmtree(project_dir)

        put_test_data(project_dir_source + '/', project_dir)

        with chdir(project_dir):

            link_test_data(
                'events/%s/' % event_name, 'data/events/%s/' % event_name)

            env = Environment([config_path, event_name])
            conf = env.get_config()

            store_ids = conf.get_elements('target_groups[:].store_id')
            for store_id in store_ids:
                store_path = 'gf_stores/%s/' % store_id
                if not os.path.exists(store_path):
                    link_test_data(store_path)

            problem_name = env.get_problem().name
            rundir_path = expand_template(
                conf.rundir_template,
                dict(problem_name=problem_name))

            return main(env, rundir_path)
