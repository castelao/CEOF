import numpy as np
from scipy import signal

def eof(data):
    """
    expects: data[time,positions]

    """
    data = signal.detrend(data, axis=0)
    U, s, V = np.linalg.svd(data, full_matrices=False) # SVD analisys
    S  = np.diag(s)
    exp_var = s**2/np.sum(s**2) # explained variance of each mode for y1 
    PC = np.dot(U,S) 
    return PC, exp_var, V


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
            #   i.e. does it take the sum or the mean? I think it is the sum.
            y[i, j] = N*(x[:,i] * x[:,j]).mean()

    return y
