import numpy as np
from scipy import signal

def eof(data):
    """
    expects: data[time,positions]

    """
    data = signal.detrend(data, axis=0)
    U, s, V = np.linalg.svd(data, full_matrices=True) # SVD analisys
    S = np.zeros((data).shape) # We need to complete with zeros
    N = S.shape[0]
    S[:N, :N] = np.diag(s)
    exp_var = s**2/np.sum(s**2) # explained variance of each mode for y1 
    PC = np.dot(S,V) 
    return PC, exp_var, U[:,0]

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
