import numpy as num
from grond import meta


def test_nanmedian():

    a = num.array([[1, 2, 3, 4, 5],
                   [1, None, 2, None, 3],
                   [None, 1, 2, 3, None]], dtype=num.float)

    res = meta.nanmedian(a, axis=0)
    for ix in range(a.shape[1]):
        b = a[:, ix]
        assert meta.nanmedian(b) == res[ix]

    res = meta.nanmedian(a, axis=1)
    for iy in range(a.shape[0]):
        b = a[iy, :]
        assert meta.nanmedian(b) == res[iy]
