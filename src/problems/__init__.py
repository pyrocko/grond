'''
Definition of the objective function and source model parameter space.
'''
import warnings

from .base import *  # noqa
from .cmt.problem import *  # noqa
from .rectangular.problem import *  # noqa
from .double_dc.problem import *  # noqa
from .volume_point.problem import *  # noqa

try:
    from .vlvd.problem import *  # noqa
except ImportError:
    warnings.warn(
        'could not import pyrocko.gf.VLVDSource. Update pyrocko to enable'
        ' inversion of VLVDProblem.', ImportWarning)
