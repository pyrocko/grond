import re
import os.path as op
from pyrocko import guts, gf
from pyrocko.guts import Bool, List, String

from .meta import Path, HasPaths, GrondError
from .dataset import DatasetConfig
from .analysers.base import AnalyserConfig
from .analysers.target_balancing import TargetBalancingAnalyserConfig
from .problems.base import ProblemConfig
from .optimisers.base import OptimiserConfig
from .targets.base import TargetGroup
from .version import __version__

guts_prefix = 'grond'


def color_diff(diff):
    green = '\x1b[32m'
    red = '\x1b[31m'
    blue = '\x1b[34m'
    dim = '\x1b[2m'
    reset = '\x1b[0m'

    for line in diff:
        if line.startswith('+'):
            yield green + line + reset
        elif line.startswith('-'):
            yield red + line + reset
        elif line.startswith('^'):
            yield blue + line + reset
        elif line.startswith('@'):
            yield dim + line + reset
        else:
            yield line


class EngineConfig(HasPaths):
    gf_stores_from_pyrocko_config = Bool.T(
        default=True,
        help='Load the GF stores from ~/.pyrocko/config')
    gf_store_superdirs = List.T(
        Path.T(),
        help='List of path hosting collection of Green\'s function stores.')
    gf_store_dirs = List.T(
        Path.T(),
        help='List of Green\'s function stores')

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
    rundir_template = Path.T(
        help='Rundir for the optimisation, supports templating'
             ' (eg. ${event_name})')
    dataset_config = DatasetConfig.T(
        help='Dataset configuration object')
    target_groups = List.T(
        TargetGroup.T(),
        help='List of ``TargetGroup``s')
    problem_config = ProblemConfig.T(
        help='Problem config')
    analyser_configs = List.T(
        AnalyserConfig.T(),
        default=[TargetBalancingAnalyserConfig.D()],
        help='List of problem analysers')
    optimiser_config = OptimiserConfig.T(
        help='The optimisers configuration')
    engine_config = EngineConfig.T(
        default=EngineConfig.D(),
        help=':class:`pyrocko.gf.LocalEngine` configuration')
    event_names = List.T(
        String.T(),
        help='Restrict application to given event names. If empty, all events '
             'found through the dataset configuration are considered.')
    event_names_exclude = List.T(
        String.T(),
        help='Event names to be excluded')

    def __init__(self, *args, **kwargs):
        HasPaths.__init__(self, *args, **kwargs)

    def get_event_names(self):
        if self.event_names:
            names = self.event_names
        else:
            names = self.dataset_config.get_event_names()

        return [name for name in names if name not in self.event_names_exclude]

    @property
    def nevents(self):
        return len(self.dataset_config.get_events())

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
        problem = self.problem_config.get_problem(
            event, self.target_groups, targets)
        self.setup_modelling_environment(problem)
        return problem

    def get_elements(self, ypath):
        return list(guts.iter_elements(self, ypath))

    def set_elements(self, ypath, value):
        guts.set_elements(self, ypath, value, regularize=True)

    def clone(self):
        return guts.clone(self)


def read_config(path):
    try:
        config = guts.load(filename=path)
    except OSError:
        raise GrondError(
            'Cannot read Grond configuration file: %s' % path)

    if not isinstance(config, Config):
        raise GrondError('Invalid Grond configuration in file "%s".' % path)

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
            header='Grond configuration file, version %s' % __version__)

        config.change_basepath(basepath)

    except OSError:
        raise GrondError(
            'Cannot write Grond configuration file: %s' % path)


def diff_configs(path1, path2):
    import sys
    import difflib
    from pyrocko import guts_agnostic as aguts

    t1 = aguts.load(filename=path1)
    t2 = aguts.load(filename=path2)

    s1 = aguts.dump(t1)
    s2 = aguts.dump(t2)

    result = list(difflib.unified_diff(
        s1.splitlines(1), s2.splitlines(1),
        'left', 'right'))

    if sys.stdout.isatty():
        sys.stdout.writelines(color_diff(result))
    else:
        sys.stdout.writelines(result)


class YPathError(GrondError):
    pass


def parse_yname(yname):
    ident = r'[a-zA-Z][a-zA-Z0-9_]*'
    rint = r'-?[0-9]+'
    m = re.match(
        r'^(%s)(\[((%s)?(:)(%s)?|(%s))\])?$'
        % (ident, rint, rint, rint), yname)

    if not m:
        raise YPathError('Syntax error in component: "%s"' % yname)

    d = dict(
        name=m.group(1))

    if m.group(2):
        if m.group(5):
            istart = iend = None
            if m.group(4):
                istart = int(m.group(4))
            if m.group(6):
                iend = int(m.group(6))

            d['slice'] = (istart, iend)
        else:
            d['index'] = int(m.group(7))

    return d


def _decend(obj, ynames):
    if ynames:
        for sobj in iter_get_obj(obj, ynames):
            yield sobj
    else:
        yield obj


def iter_get_obj(obj, ynames):
    yname = ynames.pop(0)
    d = parse_yname(yname)
    if d['name'] not in obj.T.propnames:
        raise AttributeError(d['name'])

    obj = getattr(obj, d['name'])

    if 'index' in d:
        sobj = obj[d['index']]
        for ssobj in _decend(sobj, ynames):
            yield ssobj

    elif 'slice' in d:
        for i in range(*slice(*d['slice']).indices(len(obj))):
            sobj = obj[i]
            for ssobj in _decend(sobj, ynames):
                yield ssobj
    else:
        for sobj in _decend(obj, ynames):
            yield sobj


__all__ = '''
    EngineConfig
    Config
    read_config
    write_config
    diff_configs
'''.split()
