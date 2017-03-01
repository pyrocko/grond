import multiprocessing as _mp
import signal as _signal

from .server import Baraddur  # noqa


class BaraddurProcess(_mp.Process):
    def __init__(self, *args, **kwargs):
        self.server = Baraddur(*args, **kwargs)
        self.shutdown_signal = _mp.Queue(1)
        _mp.Process.__init__(self)
        _signal.signal(_signal.SIGINT, self.server.stop)

    def run(self):
        self.server.start(signal=self.shutdown_signal)

    def stop(self):
        self.shutdown_signal.put(True)


if __name__ == '__main__':
    p = BaraddurProcess('/home/marius/Development/testing/grond/rundir')
    p.start()
