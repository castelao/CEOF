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
