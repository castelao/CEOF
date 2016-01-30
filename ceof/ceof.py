#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Class to deal with Complex EOF
"""

from UserDict import UserDict

import numpy
import numpy as np
from numpy import ma
import scipy.fftpack
import scipy.signal
from scipy import linalg

try:
    from pyclimate.svdeofs import svdeofs, getvariancefraction
except:
    print("pyclimate is not available!")

from utils import scaleEOF

"""

    Possible cases:

       - EOF
         - 2D
           - array
           - Masked Array
         - ND (EOF 2D)
           - array
           - Masked Array
       - CEOF
         - 2D (EOF 2D)
           - array
           - Masked Array
         - ND (EOF ND)
           - array
           - Masked Array

"""

"""

    Core functions for EOF decomposition
"""

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


def eof_with_gaps(x):
    """ EOF decomposition of 2D array that allows some gaps.

        It's important to keep in mind that it assumes that the gaps does not
          compromise the estimate of the covariance between the data series
          of two positions.
    """
    x = np.asanyarray(x)
    assert x.ndim == 2, "Input must be 2D"

    #Cov = np.dot(x.T, x)
    Cov = cov_with_gaps(x)
    L, C = np.linalg.eig(Cov)

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
    assert x.ndim == 2, "Input must be 2D"

    #assert x.any(axis=0).all(), \
    #        "Must have at least one valid data along axis=0"

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

    # If more than 2 dimensions, reshape and than re-reshape
    if x.ndim > 2:
        S = x.shape
        pcs, exp_var, eofs = eof_decompose(
                x.reshape((S[0], x.size/S[0]), order="C"))

        nmodes = eofs.shape[-1]
        eofs = eofs.reshape(list(S[1:]) + [nmodes], order="C")

        return pcs, exp_var, eofs
    # Below here it's expected only 2D input x.

    mask = ~np.isfinite(np.array(x)) | ma.getmaskarray(x)

    # If all data is valid, use linalg.svd is faster
    #if (np.isfinite(x).all()) & (not ma.getmaskarray(x).any()):
    if ~mask.any():
        return eof_from_svd(x)

    # If is there any column wihtout any valid data, remove those and
    #   try again.
    if mask.all(axis=0).any():
        ind = ~mask.all(axis=0)
        pcs, exp_var, tmp = eof_decompose(x[:,ind])
        eofs = ma.masked_all((ind.size, pcs.shape[1]), dtype=tmp.dtype)
        eofs[ind] = tmp
        return pcs, exp_var, eofs

    # At this point all columns have some valid data, but not all
    return eof_with_gaps(ma.fix_invalid(x))


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


def ceof_scalar2D(x):
    """ Estimate the complex EOF on a 2D array.

        Time should be the first dimension, so that the PC (eigenvalues) will
          be in respect to the first dimension.

    """
    assert type(data) is np.ndarray, \
        "ceof_scalar2D requires an ndarray but got: %s" % type(data)
    assert np.isfinite(data).all(), \
        "ceof_scalar2D requires a full valid values array"

    # ---- Creating the complex field using Hilbert transform
    input_H = numpy.empty(data.shape, dtype=data.dtype)
    for i in range(data.shape[1]):
        input_H[:,i] = scipy.fftpack.hilbert(data[:,i])

    U = data + 1j*input_H

    pcs, lambdas, eofs = svdeofs(U)

    return pcs, lambdas, eofs


def CEOF_2D(data, cfg=None):
    """ Complex EOF of a scalar 2D array
    """
    assert data.ndim == 2, "CEOF_2D requires a 2D ndarray input"

    if cfg is None:
        cfg = {'cumvar':1, 'normalize':'pc_median'}
    assert type(cfg) is dict, "cfg must be a dictionary"

    # ---- Creating the complex field using Hilbert transform
    #if self['input'].mask.any():
    #    print "There are masked values in U at CEOF_2D()"

    pcs, lambdas, eofs = ceof_scalar2D(data)

    expvar = getvariancefraction(lambdas)

    # Define how many modes will be returned by the explainned variance.
    # cumvar = 1 means 100%, i.e. all modes
    if cfg['cumvar'] == 1:
        nmodes = len(lambdas)
    else:
        # This doesn't work. Need to improve it.
        nmodes = (np.arange(len(expvar))[np.cumsum(expvar)>=cfg['cumvar']])[0]

    if 'maxnmodes' in cfg:
        nmodes = min(nmodes, cfg['maxnmodes'])

    print "Considering the first %s of %s modes." % (nmodes,len(lambdas))

    # ---- Normalize -----------------------------------------------------
    if 'normalize' in  cfg:
        pcs, eofs = scaleEOF(pcs, eofs, scaletype=cfg['normalize'])

    output = {}
    output['eofs'] = eofs[:,:nmodes]
    output['pcs'] = pcs[:,:nmodes]
    output['lambdas'] = lambdas[:nmodes]
    output['variancefraction'] = expvar[:nmodes]

    return output


#def CEOF_ND(data, cfg=None):

class CEOF(UserDict):
    """
    """
    def __init__(self, input, metadata={}, logger=None, **keywords):
        """
            Time should be the first dimension, i.e. axis=0
        """
        self.input = input.copy()
        self.data = input.copy()
        self.metadata = metadata

        self.go()

        # Save ceof
        if outputfilename is not None:
            save_ceof(self.data, outputfilename, self.nmodes)

        return


    def select_data(self, var, polygon_coordinates):
        """
        """
        #var = 'ssh'
        T,I,J = self.data[var].shape
        tmp=numpy.ones((J,K))==1
        for i in range(I):
            for j in range(J):
                   tmp[i,j] = ((self.data[var].mask)[:,i,j]).all()==False


        from shapely.geometry import Polygon
        from shapely.geometry import Point
        polygon = Polygon(polygon_coordinates)

        ind = ind&tmp
        return ind


    def go(self):
        var = self.metadata['ceof']['var']

        if ('Lat' not in self.keys()) or ('Lon' not in self.keys()):
            self.data['Lon'], self.data['Lat'] = numpy.meshgrid(self.data['lon'],self.data['lat'])

        # ---- Normalize -----------------------------------------------------
        #self.data['ssh']=self.data['ssh']-self.data['ssh'].mean()
        # --------------------------------------------------------------------
        ind = ma.getmaskarray(self.data[var]).any(axis=0) == False

        I, J, K = self.data[var].shape
    
        if 'ceof_coord' in self.metadata:
            coord = self.metadata['ceof_coord']
            assert type(coord) == list
            from shapely.geometry import Point, Polygon
            polygon = Polygon(coord)
            for j in range(J):
                for k in range(K):
                        ind[j,k] = ind[j,k] & polygon.intersects(
                                Point(self.data['Lon'][j,k],
                                    self.data['Lat'][j,k]))
  
        #ind = ind&tmp
    
        N = int((numpy.ones(ind.shape)[ind]).sum())
        grid_index = ma.masked_all((N,2), dtype='i')
        n = -1
        for j in range(J):
            for k in range(K):
                if ind[j, k]:
                    n += 1
                    grid_index[n, :] = [j, k]
    
        self.grid_index = grid_index
        data2D = numpy.zeros((I,N), dtype=self.data[var].dtype)
        for n, ind in enumerate(self.grid_index):
            data2D[:, n] = self.data[var][:,ind[0], ind[1]]

        print "Running CEOF_2D()"
        output = CEOF_2D(data2D, cfg=self.metadata['ceof'])

        nmodes = len(output['lambdas'])

        for k in [k for k in output.keys() if k is not 'ceof']:
            self.data[k] = output[k]

        self.data['eofs'] = ma.masked_all((J, K, nmodes), dtype='complex128')
        for n,ind in enumerate(self.grid_index):
            self.data['eofs'][ind[0],ind[1],:] = output['eofs'][n,:]

	# ----
        self.data['L_x'] = wavelenght_from_ceof(self.data['eofs'],
                self.data['lat'], self.data['lon'])
