import multiprocessing as _mp
import signal as _signal

from .server import Baraddur, BaraddurConfig  # noqa


class BaraddurProcess(_mp.Process):
    def __init__(self, *args, **kwargs):
        self.server = Baraddur(*args, **kwargs)
        self.shutdown_signal = _mp.Queue(1)
        _mp.Process.__init__(self)

    def run(self):
        _signal.signal(_signal.SIGINT, self.server.stop)
        self.server.start(signal=self.shutdown_signal)

    def stop(self):
        self.shutdown_signal.put(True)


if __name__ == '__main__':
    p = BaraddurProcess(project_dir='/home/marius/Development/testing/grond')
    p.start()
