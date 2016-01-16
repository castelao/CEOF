import numpy as np
from numpy import ma


def eof_from_svd(data):
    """
    expects: data[time,positions]

    """
    U, s, V = np.linalg.svd(data, full_matrices=False)  # SVD analisys
    S = np.diag(s)
    exp_var = s**2 / np.sum(s**2)  # explained variance of each mode for y1
    PC = np.dot(U, S)

    return PC, exp_var, V.T


def eof_from_eig(x):
    """ Didatic eof decomposition using eig()
    """
    R = np.dot(x.T, x)
    L, C = np.linalg.eig(R)
    eofs = C.T

    exp_var = L/L.sum()

    nmodes = L.size
    pcs = np.empty((x.shape[0], nmodes))
    for n in range(nmodes):
        pcs[:, n] = np.dot(x, C[:, n])

    return pcs, exp_var, eofs


def eof2D_gap(x):
    """ EOF decomposition of 2D array that allows some gaps.

        It's important to keep in mind that it assumes that the gaps does not
          compromise the estimate of the covariance between the data series
          of two positions.
    """
    R = cov2D_gap(x)
    L, C = np.linalg.eig(R)
    eofs = C.T

    exp_var = L/L.sum()

    nmodes = L.size
    pcs = np.empty((x.shape[0], nmodes))
    for n in range(nmodes):
        pcs[:, n] = np.dot(x, C[:, n])

    return pcs, exp_var, eofs


def cov2D_gap(x):
    """ Covariance matrix allowing gaps


        ATENTION: Might be a good idea to set a constraint on max gaps
          allowed, like at least 98% of the data?
    """
    x = np.asanyarray(x)

    assert x.ndim == 2

    N, J = x.shape

    y = ma.masked_all((J, J), x.dtype)

    for i in range(J):
        for j in range(J):
            # FIXME: Double check if I should use mean or N*mean
            #   i.e. does it take the sum or the mean? I think it is the sum.
            y[i, j] = N*(x[:, i] * x[:, j]).mean()

    return y


def eof_gap(x):
    """ EOF decomposition of nD array that allows some gaps.

        I'm not happy with this. Improve it at some point
    """
    x = np.asanyarray(x)

    if x.ndim == 2:
        return eof2D_gap(x)

    # First reshape to 2D
    S = x.shape
    x = x.reshape(S[0], x.size/S[0])
    pcs, exp_var, eofs = eof2D_gap(x)

    # Remap into the original form. So EOF would have the same dimensions
    #  of x, but without the first one that gone into pcs.
    x = x.reshape(S)
    nmodes = eofs.shape[-1]
    eofs = eofs.reshape(list(S[1:]) + [nmodes])

    return pcs, exp_var, eofs
