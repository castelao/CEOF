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

def CEOF_2D_limited(input,metadata={'variancefraction_explainned':0.95}):
        #self.data['input']=input.copy
        ceof = CEOF_2D(input)
        if 'nmodes' not in metadata:
            nmodes=len(ceof['variancefraction'])
        if 'variancefraction_explainned' in metadata:
            nmodes=(numpy.ones(ceof['variancefraction'].shape)[numpy.cumsum(ceof['variancefraction'])<=metadata['variancefraction_explainned']]).sum().astype('i')
        if 'nmodes_max' in metadata:
            if metadata['nmodes_max']<nmodes:
                nmodes=metadata['nmodes_max']
        nt,ni,nj=input.shape
        filtred=numpy.zeros((nt,ni,nj))
        for t in range(nt):
            #for n in range(20):
            for n in range(ceof['pcs'].shape[1]):
                filtred[t,:,:]+=ceof['eofs'][:,:,n].real*x['pcs'][t,n].real
        #self.data['filtred']=filtred
        return filtred


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
        return

    def filter(self,var,l,type,l2=None):
        #from maud import window_mean
        from maud import window_1Dmean_grid, get_halfpower_period
        from datetime import timedelta
        
        if len((set(numpy.diff(self.data['datetime'])))) !=1:
            print "Class incomplete. Can't deal with a non regular time series"
            return

        dt=self.data['datetime'][1]-self.data['datetime'][0]
        if type == 'bandpass':
            tscale = dt.days+dt.seconds/86400.
            ll = (l.days+l.seconds/86400.)/tscale
            ll2 = (l2.days+l2.seconds/86400.)/tscale
            #lowpass = window_mean.window_1Dmean_grid(self.data[var], ll/2., method='hann', axis=0)
            lowpass = window_1Dmean_grid(self.data[var], ll/2., method='hann', axis=0)
            output = window_1Dmean_grid(lowpass, ll2/2., method='hann', axis=0)
            output = lowpass - output

            print "ATENTION!!!! Improve this here!!!"
            self.halfpower_period = "20-120"
        else:
            ll=(l.days+l.seconds/86400.)/(dt.days+dt.seconds/86400.)

            if ll<1:
                print "This filter will have no effect. Data series have not enough resolution."
                return

            lowpass = window_1Dmean_grid(self.data[var],ll/2.,method='hann',axis=0)
            #lowpass=window_mean.window_1Dmean(self.data[var],ll,method='hanning',axis=0)
            if type=='lowpass':
                output=lowpass
            elif type=='highpass':
                #self.data[var]=x_highpass=self.data[var]-(lowpass-lowpass.mean())
                output=self.data[var]-(lowpass)
                #output=self.data[var]-lowpass
            else:
                print "On function filter, type must be lowpass or highpass"

            halfpower_period = get_halfpower_period(self.data[var],
                    output, dt=dt)
            print "Filter half window size: %s" % l
            print "Half Power Period: %s" % halfpower_period
            self.halfpower_period = halfpower_period
        
        # ----
	## I should move this to inside the window_mean_1D_grid
	#nt,ni,nj = self.data[var].shape
	#gain = ma.masked_all((nt,ni,nj))
	#for i in range(ni):
	#    for j in range(nj):
	#        if output[:,i,j].mask.all()==False:
        #            gain[:,i,j] = numpy.absolute(numpy.fft.fft(output[:,i,j]-output[:,i,j].mean())) / numpy.absolute(numpy.fft.fft(self.data[var][:,i,j]-self.data[var][:,i,j].mean()))
	#gain_median = ma.masked_all(nt)
	#for t in range(nt):
	#    gain_median[t] = numpy.median(gain[t,:,:].compressed()[numpy.isfinite(gain[t,:,:].compressed())])
	#freq=numpy.fft.fftfreq(nt)/dt.days
	#import rpy2.robjects as robjects
	#smooth = robjects.r['smooth.spline'](robjects.FloatVector(gain_median[numpy.ceil(nt/2.):]),robjects.FloatVector(-freq[numpy.ceil(nt/2.):]),spar=.4)
	##smooth = robjects.r['smooth.spline'](robjects.FloatVector(-freq[numpy.ceil(nt/2.):]),robjects.FloatVector(gain_median[numpy.ceil(nt/2.):]),spar=.4)
	#s_interp = robjects.r['predict'](smooth,x=0.5)
	#halfpower_period = 1./s_interp.rx2['y'][0]

        # ----
        self.data[var]=output

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
        if 'prefilter' in self.metadata:
            print "Filtering in time"
            if self.metadata['prefilter']['type'] == 'bandpass':
                self.filter(var,l=self.metadata['prefilter']['l'],type=self.metadata['prefilter']['type'], l2=self.metadata['prefilter']['l2'],)
            else:
                self.filter(var,l=self.metadata['prefilter']['l'],type=self.metadata['prefilter']['type'])
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
	self.set_wavelenght()

        #for k in self.data.keys():
        #    print k, self.data[k].shape, type(self.data[k])

        if 'figs' in self.metadata:
	    print "Creating figures for %s modes" % nmodes
            for n in range(nmodes):
                if 'suffix' in self.metadata['figs']:
                    filename="../figs/CEOF_%s_mode%s.eps" % (self.metadata['figs']['suffix'],(n+1))
                else:
                    filename="../figs/CEOF_mode%s.eps" % (n+1)
                limits={'LatIni':-5, 'LatFin':15, 'LonIni':-60, 'LonFin':-25}
                #self.plot(self['eofs'][:,:,n],self['pcs'][:,n],(n+1),self['variancefraction'][n],filename=filename,limits=limits,cumvarfrac=self['variancefraction'][:(n+1)].sum())
                import graphics
                graphics.plot(self['eofs'][:,:,n], self['pcs'][:,n],
                        (n+1),self['variancefraction'][n],
                        filename=filename,
                        data = self.data,
                        limits=limits,
                        cumvarfrac=self['variancefraction'][:(n+1)].sum())
