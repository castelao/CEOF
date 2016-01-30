#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from numpy import ma
from numpy import pi

t = np.arange(10)
x = np.arange(3)

#y1 = np.sin(2 * pi * t) * np.ones(x.size) + 0.1 * np.sin(2 * pi * t/2) * x.T
#y2=sin(2*pi*t)*ones(size(x’))+0.6*randn(size(y1));

#y1=sin(2*pi*t)*ones(size(x'))+0.1*sin(2*pi*t/2)*x'; >> y2=sin(2*pi*t)*ones(size(x’))+0.6*randn(size(y1));

from ceof.ceof import eof_from_svd, eof_with_gaps, eof_decompose


def check_eof_shape(x):
    """ Check if the outputs of eof have the correct shape

        This first dimension of the input will be the PC dimension
    """
    p, v, e = eof_from_svd(x)
    n, m = x.shape
    nmodes = min(n, m)
    assert p.shape == (n, nmodes)
    assert v.shape == (nmodes,)
    assert e.shape == (m, nmodes)


def compare_eof_approach(X):
    pc1, e1, eof1 = eof_from_svd(X)
    pc2, e2, eof2 = eof_with_gaps(X)

    assert np.allclose(pc1, pc2)
    assert np.allclose(e1, e2)
    assert np.allclose(eof1, eof2)


def test_eof_shapes():
    check_eof_shape(np.random.random((3, 3)))

    check_eof_shape(np.random.random((10, 3)))

    check_eof_shape(np.random.random((3, 4)))


def eof_reconstruct(x):
    """ Decompose and reconstruct the original data
    """
    p, v, e = eof_from_svd(x)
    y = np.zeros_like(x)
    for n in range(len(v)):
        y += (p[:, n][:,None] * e[:, n])

    return y


def test_eof_recons():
    x = np.random.random((4,3))
    assert np.allclose(x, eof_reconstruct(x))

    # FIXME: Create another test for this signal
    # y = np.sin(2 * pi * t)[:, None] + (0.1 * np.sin(2 * pi * t/2))[:, None] * x


def test_sorted_lambda():
    """ Check if eigenvalues are sorted
    """
    for i in range(10):
        X = np.random.random((10,4))
        pcs, exp_var, eofs = eof_with_gaps(X)
        assert (np.diff(exp_var) <= 0).all()

t = np.arange(100)
x = np.arange(-10,11)

def test_linear_sin_2D():

    y = np.sin(2*pi/30 *t)[:, None] * np.ones_like(x)
    noise = y.std()*0.5*(np.random.random(y.shape) - 0.5)
    pcs, exp_var, eofs = eof_decompose(y + noise)
    np.corrcoef(y[:,0], pcs[:,0])
    assert np.abs(np.corrcoef(y[:,0], pcs[:,0])[0,1]) > 0.99


def test_2linear_sin_2D():

    y = np.sin(2*pi/30 *t)[:, None] * np.array([2, 1])
    noise = y.std()*0.3*(np.random.random(y.shape) - 0.5)
    pcs, exp_var, eofs = eof_decompose(y + noise)
    assert np.abs(np.corrcoef(y[:,0], pcs[:,0])[0,1]) > 0.99
    # Bad test below.
    assert np.abs(2 - np.abs(eofs[0,0] / eofs[1,0]))<0.2

    
def test_linear_sin_2D_missingposition():    
    #y = np.sin(2*pi/30 *t)[:, None] * np.ones_like(x)
    y = np.sin(2*pi/30 *t)[:, None] * np.array([2, 1, 4])
    noise = y.std()*0.3*(np.random.random(y.shape) - 0.5)
    mask = np.zeros(y.shape, dtype='bool')
    mask[:,2] = True
    pcs, exp_var, eofs = eof_decompose(
            ma.masked_array(y + noise, mask = mask))
    
    assert np.abs(np.corrcoef(y[:,0], pcs[:,0])[0,1]) > 0.99
    assert eofs.mask[2,:].all()
    assert ~eofs.mask[:2,:].any()
    assert (eofs.mask.any(axis=1) == [False, False, True]).all()

def test_linear_sin_3D():
    X, Y = np.meshgrid(np.arange(-10, 11), np.arange(-10, 11))
    Z = np.sin(2*pi/30 *t)[:, None, None] * np.ones(X.shape)
    #noise = Z.std()*0.3*(np.random.random(Z.shape) - 0.5)
    noise = 0.1*(np.random.random(Z.shape) - 0.5)
    pcs, exp_var, eofs = eof_decompose(Z + noise)
    assert exp_var[0] > 0.99
    assert np.abs(np.corrcoef(Z[:,0,0], pcs[:,0])[0,1]) > 0.99

def test_monopole_sin_3D():
    X, Y = np.meshgrid(np.arange(-10, 11), np.arange(-10, 11))
    R = (X**2 + Y**2) **0.5
    pattern = np.sin(2*pi/50*R)
    Z = np.sin(2*pi/30 *t)[:, None, None] * pattern[None, :, :] 
    #noise = Z.std()*0.3*(np.random.random(Z.shape) - 0.5)
    noise = 0.1*(np.random.random(Z.shape) - 0.5)
    pcs, exp_var, eofs = eof_decompose(Z + noise)
    assert exp_var[0] > 0.99
    assert np.abs(np.corrcoef(Z[:,0,0], pcs[:,0])[0,1]) > 0.99
    assert np.abs(np.corrcoef(np.ravel(pattern), np.ravel(eofs[:,:,0]))[0,1]) > 0.99


def test_maskposition():
    X, Y = np.meshgrid(np.arange(-10, 11), np.arange(-10, 11))
    R = (X**2 + Y**2) **0.5
    pattern = np.sin(2*pi/50*R)
    Z = np.sin(2*pi/30 *t)[:, None, None] * pattern[None, :, :] 
    #noise = Z.std()*0.3*(np.random.random(Z.shape) - 0.5)
    noise = 0.1*(np.random.random(Z.shape) - 0.5)
    data = Z + noise
    pcs, exp_var, eofs = eof_decompose(data)

    mask = np.zeros(data.shape, dtype='bool')
    mask[:,8,8] = True
    data2 = ma.masked_array(data, mask)
    pcs2, exp_var2, eofs2 = eof_decompose(data2)

    assert np.abs(np.corrcoef(pcs[:,0], pcs2[:,0])[0,1]) > 0.99
    assert np.abs(exp_var[0] - exp_var2[0]) < 1e-3
    ind = ~data2.mask.all(axis=0)
    assert np.abs(np.corrcoef(eofs[:,:,0][ind], eofs2[:,:,0][ind])[0,1])>0.99
