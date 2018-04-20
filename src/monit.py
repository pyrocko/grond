import time
import logging
import threading

from .problems.base import ModelHistory, load_optimizer_info

logger = logging.getLogger('grond.monit')


class GrondMonit(object):

    def __init__(self, rundir='.', interval=1.):
        logger.info('Waiting for %s to come alive...' % rundir)
        self.interval = interval
        self.rundir = rundir

        self.history = ModelHistory.follow(rundir)
        self.optimizer = load_optimizer_info(rundir)
        self.problem = self.history.problem

        self.starttime = time.time()

        self.iter_per_second = 0.
        self.optimizer_status = None
        self.last_update = None

        self._prev_iiter = 0

    @property
    def runtime(self):
        return timedelta(seconds=round(time.time() - self.starttime))

    @property
    def runtime_remaining(self):
        if self.iter_per_second == 0.:
            return timedelta()
        return timedelta(seconds=round((self.niter - self.iiter)
                         / self.iter_per_second))

    def update(self):
        self.history.update()

        self.optimizer_status = self.optimizer.get_status(self.history)

        if self._prev_iiter:
            self.iter_per_second = \
                (self.optimizer_status['iiter'] - self._prev_iiter) / \
                (time.time() - self.last_update)
            self._prev_iiter = self.optimizer_status['iiter']

        self.optimizer_status['iter_per_second'] = self.iter_per_second
        self.last_update = time.time()

    def print_status(self):
        print(self.optimizer_status)

    def run(self):
        while True:
            time.sleep(self.interval)
            self.update()
            self.print_status()

    @classmethod
    def start(cls, rundir='.'):
        monitor = cls(rundir)
        monit_thread = threading.Thread(target=monitor.run)
        print('starting thread!')
        monit_thread.start()
        return monit_thread
