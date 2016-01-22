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
