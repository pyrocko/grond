import sys
import copy
import difflib
from pyrocko import guts_agnostic as aguts
from .config import color_diff


def rename_attribute(old, new):
    def func(path, obj):
        if old in obj:
            obj.rename_attribute(old, new)

    return func


def rename_class(new):
    def func(path, obj):
        obj._tagname = new

    return func


def drop_attribute(old):
    def func(path, obj):
        if old in obj:
            obj.drop_attribute(old)

    return func


def set_attribute(k, v):
    def func(path, obj):
        obj[k] = v

    return func


def upgrade_config_file(fn, diff=True):
    rules = [
        ('grond.HighScoreOptimizerConfig',
            rename_class('grond.HighScoreOptimiserConfig')),
        ('grond.Config',
            rename_attribute('optimizer_config', 'optimiser_config')),
        ('grond.CMTProblemConfig',
            drop_attribute('apply_balancing_weights')),
        ('grond.DoubleDCProblemConfig',
            drop_attribute('apply_balancing_weights')),
        ('grond.RectangularProblemConfig',
            drop_attribute('apply_balancing_weights')),
        ('grond.Config',
            drop_attribute('analyser_config')),
        ('grond.Config',
            set_attribute(
                'analyser_configs',
                aguts.load(string='''
                    - !grond.TargetBalancingAnalyserConfig
                      niterations: 1000
                    '''))),

    ]

    def apply_rules(path, obj):
        for tagname, func in rules:
            if obj._tagname == tagname:
                func(path, obj)

    t1 = aguts.load(filename=fn)
    t2 = copy.deepcopy(t1)

    aguts.apply_tree(t2, apply_rules)

    s1 = aguts.dump(t1)
    s2 = aguts.dump(t2)

    if diff:
        result = list(difflib.unified_diff(
            s1.splitlines(1), s2.splitlines(1),
            'normalized old', 'normalized new'))

        if sys.stdout.isatty():
            sys.stdout.writelines(color_diff(result))
        else:
            sys.stdout.writelines(result)
    else:
        print(aguts.dump(t2, header=True))


if __name__ == '__main__':
    fn = sys.argv[1]
    upgrade_config_file(fn)
