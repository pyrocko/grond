import logging
import tempfile
import os
import shutil

from grond import dataset
from grond import problem  # noqa
from grond import targets  # noqa
from grond import optimizers  # noqa

from pyrocko import util

url = 'http://data.pyrocko.org/testing/grond/'

logger = logging.getLogger('grond.test')
op = os.path

test_data = op.join(op.abspath(__file__), 'data')


def get_test_data(self, path):
    fn = op.join(test_data, path)

    if not os.exists(test_data):
        os.mkdir(test_data)

    if not os.exists(fn):
        if path.endswith('/'):
            util.download_dir(op.join(url, path), fn)
        else:
            util.download_file(op.join(url, path), fn)

    return fn


def get_dateset_config():

    def create_test_dataset(self, *args, **kwargs):
        dataset.DatasetConfig.__init__(self, *args, **kwargs)

        tmpdir = tempfile.mkdtemp(prefix='grond')
        dirs = ['gf_store/', 'waveforms/', 'insar/', 'gps/', 'meta/']

        for d in dirs:
            os.mkdir(op.join(tmpdir, d))
            shutil.copy(os.get_test_data(d), tmpdir)

        self.path_prefix = tmpdir
        self.waveform_paths = ['waveforms/']
        self.stations_stationxml_paths = ['meta/station.xml']
        self.responses_stationxml_paths = ['meta/station.xml']
        self.picks_path = ['meta/picks.txt']

        self.set_basepath(tmpdir)
        logger.info('Initialising new dataset')

    def delete_test_dataset(self):
        os.rmdir(self.tmpdir)

    DatasetConfig = dataset.DatasetConfig
    DatasetConfig.__init__ = create_test_dataset
    return DatasetConfig


def get_grond_target_config():
    pass


def get_grond_configs():
    pass
