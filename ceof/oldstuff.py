# Old stuff that I'm not sure yet about what to do with it


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


"""
    This was inside the class CEOF, in the function go().
    The right place for it would be as a function to plot the figures and it would be loaded with the object [C]EOF.


    def filter(self,var,l,type,l2=None):
        #from maud import window_mean
        #from maud import window_1Dmean_grid, get_halfpower_period
        from maud import wmean_1D, wmean_bandpass_1D
        from datetime import timedelta

        if len((set(numpy.diff(self.data['datetime'])))) !=1:
            print "Class incomplete. Can't deal with a non regular time series"
            return

        #d = self.data['datetime'] - self.data['datetime'].min()
        t = np.array([d.days+d.seconds/86400. for d in \
                self.data['datetime'] - self.data['datetime'].min()])
        #dt=self.data['datetime'][1]-self.data['datetime'][0]
        l = l.days + l.seconds/86400.
        if type == 'bandpass':
            #tscale = dt.days+dt.seconds/86400.
            #ll = (l.days+l.seconds/86400.)/tscale
            #ll2 = (l2.days+l2.seconds/86400.)/tscale
            l2 = l2.days + l2.seconds/86400.
            #lowpass = window_mean.window_1Dmean_grid(self.data[var], ll/2., method='hann', axis=0)
            lowpass = wmean_1D(self.data[var], l2, t=t, method='hann', axis=0)
            #lowpass = window_1Dmean_grid(self.data[var], ll/2., method='hann', axis=0)
            #output = window_1Dmean_grid(lowpass, ll2/2., method='hann', axis=0)
            output = wmean_1D(lowpass, l, t=t, method='hann', axis=0)
            #output = lowpass - output
            output = wmean_bandpass_1D(self.data[var], lshorterpass=l2, llongerpass=l, t=t, method='hann', axis=0)

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


    def go():
        if 'prefilter' in self.metadata:
            import pdb; pdb.set_trace()
            print "Filtering in time"
            if self.metadata['prefilter']['type'] == 'bandpass':
                self.filter(var, l=self.metadata['prefilter']['l'],
                        type=self.metadata['prefilter']['type'],
                        l2=self.metadata['prefilter']['l2'])
            else:
                self.filter(var, l=self.metadata['prefilter']['l'],
                        type=self.metadata['prefilter']['type'])

        if 'figs' in self.metadata:
            print "Creating figures for %s modes" % self.nmodes
            for n in range(self.nmodes):
                if 'suffix' in self.metadata['figs']:
                    filename="../figs/CEOF/CEOF_%s_mode%s.png" % (self.metadata['figs']['suffix'],(n+1))
                else:
                    filename="../figs/CEOF/CEOF_mode%s.png" % (n+1)

                limits={'LatIni':-2.5, 'LatFin':15, 'LonIni':-62.5, 'LonFin':-35}
                graphics.plot(self['eofs'][:,:,n], self['pcs'][:,n],
                        (n+1),self['variancefraction'][n],
                        filename=filename,
                        data = self.data,
                        limits=limits,
                        cumvarfrac=self['variancefraction'][:(n+1)].sum())
"""
