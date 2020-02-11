from __future__ import print_function

import sys
import copy
import difflib
from pyrocko import guts_agnostic as aguts
from .config import color_diff
import grond


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


def to_aguts(obj):
    return aguts.load(string=obj.dump())


def upgrade_solver(path, obj):

    stages = [
        grond.UniformSamplerPhase(
            niterations=obj['niter_uniform'])]

    if obj.get('niter_transition', 0) != 0:
        stages.append(
            grond.DirectedSamplerPhase(
                niterations=obj['niter_transition'],
                scatter_scale_begin=obj.get('scatter_scale_transition', 2.0),
                scatter_scale_end=obj.get('scatter_scale', 1.0),
                sampler_distribution=obj.get('sampler_distribution'),
                standard_deviation_estimator=obj.get(
                    'standard_deviation_estimator',
                    'median_density_single_chain')))

    if obj.get('niter_explorative', 10000) != 0:
        stages.append(
            grond.DirectedSamplerPhase(
                niterations=obj.get('niter_explorative', 10000),
                scatter_scale_begin=obj.get('scatter_scale', 1.0),
                scatter_scale_end=obj.get('scatter_scale', 1.0),
                sampler_distribution=obj.get(
                    'sampler_distribution', 'multivariate_normal'),
                standard_deviation_estimator=obj.get(
                    'standard_deviation_estimator',
                    'median_density_single_chain'),
                starting_point='excentricity_compensated' if obj.get(
                    'compensate_excentricity', True) else 'random'))

    if obj.get('niter_non_explorative', 0) != 0:
        stages.append(
            grond.DirectedSamplerPhase(
                niterations=obj['niter_non_explorative'],
                scatter_scale_begin=obj.get('scatter_scale', 1.0),
                scatter_scale_end=obj.get('scatter_scale', 1.0),
                sampler_distribution=obj.get(
                    'sampler_distribution',
                    'multivariate_normal'),
                standard_deviation_estimator=obj.get(
                    'standard_deviation_estimator',
                    'median_density_single_chain'),
                starting_point='mean'))

    optimiser_config = grond.HighScoreOptimiserConfig(
        sampler_phases=stages)

    obj.replace(to_aguts(optimiser_config))


def upgrade_group(path, obj):
    if 'super_group' in obj or 'group' in obj:
        obj['path'] = obj['super_group'] + '.' + obj['group']
        obj['normalisation_family'] = obj['super_group']
        obj.drop_attribute('super_group')
        obj.drop_attribute('group')


def upgrade_nbootstrap(path, obj):
    if 'nbootstrap' in obj['problem_config']:
        obj['optimiser_config']['nbootstrap'] = obj['problem_config']['nbootstrap']
        obj['problem_config'].drop_attribute('nbootstrap')


def upgrade_config_file(fn, diff=True):
    rules = [
        ('grond.HighScoreOptimizerConfig',
            rename_class('grond.HighScoreOptimiserConfig')),
        ('grond.Config',
            rename_attribute('target_configs', 'target_groups')),
        ('grond.TargetConfig',
            rename_attribute('inner_misfit_config', 'misfit_config')),
        ('grond.TargetConfig', upgrade_group),
        ('grond.TargetConfig',
            rename_class('grond.WaveformTargetGroup')),
        ('grond.InnerMisfitConfig',
            rename_class('grond.WaveformMisfitConfig')),
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
        ('grond.Config',
            rename_attribute('solver_config', 'optimiser_config')),
        ('grond.SolverConfig', upgrade_solver),
        ('grond.Config', upgrade_nbootstrap),
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
