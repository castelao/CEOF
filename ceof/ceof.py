#!/usr/bin/env python
# -*- coding: Latin-1 -*-

""" Class to deal with Complex EOF
"""

from UserDict import UserDict

import numpy
import numpy as np
from numpy import ma
import scipy.fftpack

try:
    from pyclimate.svdeofs import svdeofs, getvariancefraction
except:
    print("pyclimate is not available!")

from utils import scaleEOF


def ceof_scalar2D(data):
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


class CEOF_Filter():
    """ Unfinished
    
        This will be a class to filter using EOF. One define the number of modes
	  or the desired variability to be explainned. The field is decomposed
	  by CEOF, than the field is reconstructed considering only the n first
	  modes.



    """
    def __init__(self,input,metadata={'variancefraction_explainned':0.95}):
        """
	"""
	N = pcs.shape[1]    # Now it is all the modes.
        T = pcs.shape[0]
	I = eofs.shape[0]
	J = eofs.shape[1]
        data_filtered = numpy.zeros((T, I, J))

        eof_amp=(eof.real**2+eof.imag**2)**0.5
        eof_phase=numpy.arctan2(eof.imag,eof.real)
        pc_amp = (numpy.real(pc)**2+numpy.imag(pc)**2)**0.5
        pc_phase = numpy.arctan2(numpy.imag(pc),numpy.real(pc))	

	for t in range(T):
            for n in range(N):
	        data_filtered[t] = data_filtered[t] + eof_amp[:,:,n]*pc_amp[t,n]*numpy.cos(eof_phase[:,:,n]+pc_phase[t,n])
        return



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
