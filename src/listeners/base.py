from threading import Thread


class Listener(Thread):

    def __init__(self):
        Thread.__init__(self)

    def progress_start(self, name, niter):
        raise NotImplementedError()

    def progress_finish(self, name):
        raise NotImplementedError()

    def progress_update(self, name, iiter):
        raise NotImplementedError()

    def state(self, state):
        raise NotImplementedError()
