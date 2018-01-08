import logging

from pyrocko import gf

from ..base import MisfitTarget, MisfitConfig, MisfitResult, TargetGroup

guts_prefix = 'grond'
logger = logging.getLogger('grond.targets.gnss.target')


class GNSSTargetGroup(TargetGroup):
    pass


class GNSSMisfitResult(MisfitResult):
    pass


class GNSSMisfitConfig(MisfitConfig):
    pass


class GNSSMisfitTarget(gf.StaticTarget, MisfitTarget):
    pass
