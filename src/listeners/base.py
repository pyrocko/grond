class Listener(object):

    def progress_start(self, name, niter):
        raise NotImplementedError()

    def progress_finish(self, name):
        raise NotImplementedError()

    def progress_update(self, name, iiter):
        raise NotImplementedError()

    def state(self, state):
        raise NotImplementedError()
