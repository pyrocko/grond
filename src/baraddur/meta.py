import numpy as num


class ColorCycler(list):

    def __init__(self, *args, **kwargs):
        list.__init__(self, *args, **kwargs)
        self.index = -1

    def __next__(self):
        self.index += 1
        if self.index >= len(self):
            self.index = 0
        return self[self.index]


def makeColorGradient(misfits, fr=1., fg=.5, fb=1.,
                      pr=0, pg=2.5, pb=4):
    misfits /= misfits.max()
    r = num.sin(fr * misfits + pr) * 127 + 128
    g = num.sin(fg * misfits + pg) * 127 + 128
    b = num.sin(fb * misfits + pb) * 127 + 128
    return ['#%02x%02x%02x' % (r[i], g[i], b[i]) for i in range(misfits.size)]
