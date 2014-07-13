import numpy as np

def scaleEOF(pcs, eofs, scaletype):
    """ Scale the EOFS and PCS preserving the mode

        If scaled by the pcmedian:
            pc = pc/median(pc)
            eof = eof*median(pc)
        This is the scale that I most use. On this case, the EOF structure
          has that magnitude at least half of the timeseries.
    """

    assert pcs.ndim == 2
    assert eofs.ndim == 2, "Expected a 2D EOF array."
    assert pcs.shape[-1] == eofs.shape[-1]

    nmodes = pcs.shape[-1]

    if scaletype == 'pc_std':
        print "Normalizing by the pc_std"
        for n in range(nmodes):
            fac = (np.absolute(pcs[:,n])).std()
            pcs[:,n] = pcs[:,n]/fac
            eofs[:,n] = eofs[:,n]*fac

    elif scaletype == 'pc_median':
        print "Normalizing by the pc_median"
        for n in range(nmodes):
            fac = np.median((np.absolute(pcs[:,n])))
            pcs[:,n] = pcs[:,n]/fac
            eofs[:,n] = eofs[:,n]*fac

    elif scaletype == 'pc_max':
        print "Normalizing by the pc_max"
        for n in range(nmodes):
            fac = (np.absolute(pcs[:,n])).max()
            pcs[:,n] = pcs[:,n]/fac
            eofs[:,n] = eofs[:,n]*fac

    elif scaletype == 'eof_std':
        print "Normalizing by the eof_std"
        for n in range(nmodes):
            fac = (np.absolute(eofs[:,n])).std()
            eofs[:,n] = eofs[:,n]/fac
            pcs[:,n] = pcs[:,n]*fac

    elif scaletype == 'eof_max':
        print "Normalizing by the eof_max"
        for n in range(nmodes):
            fac = (np.absolute(eofs[:,n])).max()
            eofs[:,n] = eofs[:,n]/fac
            pcs[:,n] = pcs[:,n]*fac

    else:
        print "Don't understand the normalization config: %s" % scaletype
        return

    return pcs, eofs
