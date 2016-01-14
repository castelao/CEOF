# -*- coding: utf-8 -*-

import numpy as np


from ceof.ceof_gap import eof

def test_eof_shapes():
    x = np.random.random((3, 3))
    p, v, e = eof(x)
    assert p.shape == (3, 3)
    assert v.shape == (3,)
    assert e.shape == (3, 3)

    # Let's suppose 3 sensors with a time series of 10 measurements
    x = np.random.random((10, 3))
    p, v, e = eof(x)
    # It should return 3 modes
    assert p.shape == (10, 3)
    assert v.shape == (3,)
    assert e.shape == (3, 3)


def eof_reconstruct(x):
    """ Decompose and reconstruct the original data
    """
    p, v, e = eof(x)
    y = np.zeros_like(x)
    for n in v.size():
        y += (p[:, n][:,None] * e[:, n])

    return y


def test_eof_recons():
    x = np.random.random((4,3))
    assert np.all_close(x, eof_reconstruct(x))

    # FIXME: Create another test for this signal
    # y = np.sin(2 * pi * t)[:, None] + (0.1 * np.sin(2 * pi * t/2))[:, None] * x
