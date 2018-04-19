# add all modules which register plots here
import grond.optimizers.plot
import grond.optimizers.highscore.plot
import grond.targets.waveform.plot
import grond.targets.gnss_campaign.plot
import grond.targets.satellite.plot
import grond.problems.plot  # noqa

from grond.plot.config import *  # noqa
from grond.plot.main import *  # noqa
from grond.plot.plotter import *  # noqa
from grond.plot.registry import *  # noqa
