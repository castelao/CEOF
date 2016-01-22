# -*- coding: utf-8 -*-

"""

    Core functions for EOF decomposition
"""

import numpy as np
from numpy import ma


def eof_from_svd(data):
    """
    Input: 2D array

    Output: PCs, explainned variance, EOFs

    The output PC will be composed over the first dimension of input data,
      while the EOF will be composed from the second dimension. For example,
      for a series of sensors, a typical analysis is to have the time along
      the first dimension (rows) and each sensor on the second dimension
      (columns).
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

    # I learnned from Eric Firing that L is not necessarily sorted
    isort = np.argsort(L)[::-1]
    L = L[isort]
    C = C[:, isort]

    eofs = C.T

    exp_var = L/L.sum()

    nmodes = L.size
    pcs = np.empty((x.shape[0], nmodes))
    for n in range(nmodes):
        pcs[:, n] = np.dot(x, C[:, n])

    return pcs, exp_var, eofs


def eof_with_gaps(x):
    """ EOF decomposition of 2D array that allows some gaps.

        It's important to keep in mind that it assumes that the gaps does not
          compromise the estimate of the covariance between the data series
          of two positions.
    """
    R = cov_with_gaps(x)
    L, C = np.linalg.eig(R)

    # I learnned from Eric Firing that L is not necessarily sorted
    isort = np.argsort(L)[::-1]
    L = L[isort]
    C = C[:, isort]

    eofs = C.T

    exp_var = L/L.sum()

    nmodes = L.size
    pcs = np.empty((x.shape[0], nmodes))
    for n in range(nmodes):
        pcs[:, n] = np.dot(x, C[:, n])

    return pcs, exp_var, eofs


def cov_with_gaps(x):
    """ Covariance matrix allowing gaps


        ATENTION: Might be a good idea to set a constraint on max gaps
          allowed, like at least 98% of the data?
    """
    assert x.ndim == 2

    N, J = x.shape

    y = np.empty((J, J), x.dtype)
    # While I don't apply a criteria if the cov can be estimated, it will
    #   run for all points, hence a regular np.array can do it.
    # y = ma.masked_all((J, J), x.dtype)

    for i in range(J):
        for j in range(J):
            # FIXME: Double check if I should use mean or N*mean
            #   i.e. does it take the sum or the mean? I think it is the sum.
            y[i, j] = N*(x[:, i] * x[:, j]).mean()

    return y


def eof_decompose(x):
    """

    """
    x = np.asanyarray(x)


    # Probably improve this decision test
    if (type(x) is ma.MaskedArray) and (x.mask.any()):
        return eof_with_gaps(x)
    else:
        return eof_from_svd(x)


def eof_reconstruct(pcs, eofs, nmodes=None):
    """ Reconstruct the input data from pcs + eofs
    """
    # At this point works only with 2D arrays.
    assert pcs.ndim == 2
    assert eofs.ndim == 2
    # It assumes that the modes are in last dimension.
    assert pcs.shape(1) == eofs.shape(1)

    
    assert (nmodes is None) or (nmodes <= pcs.shape(1))

    if nmodes is None:
        data = np.dot(pcs, eofs.T)
    else:
        data = np.dot(pcs[:,:nmodes], eofs[:,:nmodes].T)

    return data
