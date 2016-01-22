# -*- coding: utf-8 -*-

import numpy as np


from ceof.ceof_gap import eof


def check_eof_shape(x):
    """ Check if the outputs of eof have the correct shape

        This first dimension of the input will be the PC dimension
    """
    p, v, e = eof_from_svd(x)
    n, m = x.shape
    nmodes = min(n, m)
    assert p.shape == (n, nmodes)
    assert v.shape == (nmodes,)
    assert e.shape == (m, nmodes)


def test_eof_shapes():
    check_eof_shape(np.random.random((3, 3)))

    check_eof_shape(np.random.random((10, 3)))

    check_eof_shape(np.random.random((3, 4)))


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


def test_sorted_lambda():
    """ Check if eigenvalues are sorted
    """
    for i in range(10):
        X = np.random.random((10,4))
        pcs, exp_var, eofs = eof_with_gaps(X)
        assert (np.diff(exp_var) <= 0).all()
