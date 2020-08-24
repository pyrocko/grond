import matplotlib
matplotlib.use('Agg')
from pyrocko import util  # noqa
util.force_dummy_progressbar = True
util.setup_logging('grondtest', 'info')

import warnings  # noqa
warnings.simplefilter(action='ignore', category=FutureWarning)
