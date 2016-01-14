# -*- coding: utf-8 -*-

import numpy as np


from ceof.ceof_gap import eof

def test_eof_shapes():
    x = np.random.random((3, 3))
    p, v, e = eof(x)
    assert p.shape == (3, 3)
    assert v.shape == (3,)
    assert e.shape == (3, 3)

    x = np.random.random((10, 3))
    p, v, e = eof(x)
    assert p.shape == (10, 3)
    assert v.shape == (3,)
    assert e.shape == (3, 3)


def eof_reconstruct(x):
    """ Decompose and reconstruct the original data
    """
    p, v, e = eof(x)
    y = np.zeros_like(x)
    for n in range(len(v)):
        y += (p[:, n][:,None] * e[:, n])

    return y


def test_eof_recons():
    x = np.random.random((4,3))
    assert np.allclose(x, eof_reconstruct(x))

    # FIXME: Create another test for this signal
    # y = np.sin(2 * pi * t)[:, None] + (0.1 * np.sin(2 * pi * t/2))[:, None] * x
