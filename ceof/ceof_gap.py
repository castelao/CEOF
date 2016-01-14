import numpy as np


def cov2D_gap(x):
    """ Covariance matrix allowing gaps

    """
    x = np.asanyarray(x)

    assert x.ndim == 2

    N, J = x.shape

    y = ma.masked_all((J, J), x.dtype)

    for i in range(J):
        for j in range(J):
            # FIXME: Double check if I should use mean or N*mean
            y[i, j] = (x[:,i] * x[:,j]).mean()

    return y
