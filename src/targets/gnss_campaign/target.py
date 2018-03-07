import logging

from pyrocko import gf
from pyrocko.guts import Object

from ..base import MisfitTarget, MisfitResult, TargetGroup

guts_prefix = 'grond'
logger = logging.getLogger('grond.targets.gnss.target')


class GNSSTargetGroup(TargetGroup):
    pass


class GNSSMisfitResult(MisfitResult):
    pass


class GNSSMisfitConfig(Object):
    pass


class GNSSMisfitTarget(gf.StaticTarget, MisfitTarget):
    pass
