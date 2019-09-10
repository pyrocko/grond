import matplotlib
matplotlib.use('Agg')
from pyrocko import util
util.force_dummy_progressbar = True
util.setup_logging('grondtest', 'info')

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
