# experimental block start
# @PYFAUST_LICENSE_HEADER@
## @package pyfaust.poly @brief This module provides polynomials as Faust objects.

import _FaustCorePy
import scipy.sparse as sp
import numpy as np
import _FaustCorePy
from scipy.sparse import csr_matrix
from pyfaust import (Faust, isFaust, eye as feye, vstack as fvstack, hstack as
                     fhstack)
from scipy.sparse.linalg import eigsh
import threading


def Chebyshev(L, K, ret_gen=False, dev='cpu', T0=None):
    """
    Builds the Faust of the Chebyshev polynomial basis defined on the symmetric matrix L.

    Args:
        L: the symmetric matrix.
        K: the degree of the last polynomial, i.e. the K+1 first polynomials are built.
        dev: the destination device of the polynomial Faust.
        ret_gen: to return a generator of polynomials in addition to the
        polynomial itself (the generator starts from the the
        K+1-degree polynomial, and allows this way to compute the next
        polynomial simply with the instruction: next(generator)).
        T0: to define the 0-degree polynomial as something else than the
        identity.

    Returns:
        The Faust of the K+1 Chebyshev polynomials.
    """
    if not isinstance(L, csr_matrix) and not isFaust(L):
        L = csr_matrix(L)
    twoL = 2*L
    d = L.shape[0]
    # Id = sp.eye(d, format="csr")
    Id = _eyes_like(L, d)
    if isinstance(T0, type(None)):
        T0 = Id
    T1 = _vstack((Id, L))
    rR = _hstack((-1*Id, twoL))
    if ret_gen or isFaust(L):
        g = _chebyshev_gen(L, T0, T1, rR, dev)
        for i in range(0, K):
            next(g)
        if ret_gen:
            return next(g), g
        else:
            return next(g)
    else:
        return _chebyshev(L, K, T0, T1, rR, dev)


def basis(L, K, basis_name, ret_gen=False, dev='cpu', T0=None, impl="native"):
    """
    Builds the Faust of the polynomial basis defined on the symmetric matrix L.

    Args:
        L: the symmetric matrix.
        K: the degree of the last polynomial, i.e. the K+1 first polynomials are built.
        basis_name: 'chebyshev', and others yet to come.
        dev: the destination device of the polynomial Faust.
        ret_gen: to return a generator of polynomials in addition to the
        polynomial itself (the generator starts from the the
        K+1-degree polynomial, and allows this way to compute the next
        polynomial simply with the instruction: next(generator)).
        impl: "native" (by default) for the C++ impl., "py" for the Python
        impl.

    Returns:
        The Faust of the K+1 Chebyshev polynomials.
    """
    if basis_name.lower() == 'chebyshev':
        if impl == "native":
            F = FaustPoly(core_obj=_FaustCorePy.FaustCore.polyBasis(L, K))
            if ret_gen:
                g = F._generator()
                return F, g
            else:
                return F
        elif impl == "py":
            return Chebyshev(L, K, ret_gen=ret_gen, dev=dev, T0=T0)
        else:
            raise ValueError(impl+" is an unknown implementation.")
    else:
        raise ValueError(basis_name+" is not a valid basis name")



def poly(coeffs, basis='chebyshev', L=None, dev='cpu', impl='native'):
    """
        Returns the linear combination of the polynomials defined by basis.

        Args:
            coeffs: the linear combination coefficients (numpy.array).
            basis: either the name of the polynomial basis to build on L or the Faust of
            polynomials if already built externally or an equivalent np.ndarray.
            L: the symmetric matrix on which the polynomials are built, can't be None
            if basis is a function (not a Faust).
            dev: the device to instantiate the returned Faust ('cpu' or 'gpu').

        Returns:
            The linear combination Faust.
    """
    K = coeffs.size-1
    if isinstance(basis, str):
        if L is None:
            raise ValueError('The L matrix must be set to build the'
                             ' polynomials.')
        F = basis(L, K, basis, dev=dev, impl=impl)
    if isFaust(basis):
        F = basis
    elif not isinstance(basis, np.ndarray):
        print("type", type(basis))
        raise TypeError('basis is neither a str neither a Faust nor'
                        ' a numpy.ndarray')
    else:
        F = basis
    if impl == 'py':
        if isFaust(F):
            Id = sp.eye(L.shape[1], format="csr")
            scoeffs = sp.hstack(tuple(Id*coeffs[i] for i in range(0, K+1)),
                                format="csr")
            Fc = Faust(scoeffs, dev=dev) @ F
            return Fc
        else:
           # F is a np.ndarray
           return _poly_arr_py(coeffs, F, L, dev=dev)
    elif impl == 'native':
        if isFaust(F):
            Fc = poly_Faust_cpp(coeffs, F)
            if F.device != dev:
                Fc = Fc.clone(dev=dev)
            return Fc
        else:
            return _poly_arr_cpp(coeffs, F, L, dev='cpu')
    else:
        raise ValueError(impl+" is an unknown implementation.")

def _poly_arr_py(coeffs, basisX, L, dev='cpu'):
    """
    """
    mt = True # multithreading
    n = basisX.shape[1]
    d = L.shape[0]
    K_plus_1 = int(basisX.shape[0]/d)
    Y = np.empty((d, n))
    if n == 1:
        Y[:, 0] = basisX[:, 0].reshape(K_plus_1, d).T @ coeffs
    elif mt:
        nthreads = 4
        threads = []
        def apply_coeffs(i, n):
            for i in range(i,n,nthreads):
                Y[:, i] = basisX[:, i].reshape(K_plus_1, d).T @ coeffs
        for i in range(0,nthreads):
            t = threading.Thread(target=apply_coeffs, args=([i,n]))
            threads.append(t)
            t.start()
        for i in range(0,nthreads):
           threads[i].join()
    else:
         for i in range(n):
                Y[:, i] = basisX[:, i].reshape(K_plus_1, d).T @ coeffs
# other way:
#	Y = coeff[0] * basisX[0:d,:]
#	for i in range(1,K+1):
#		Y += (basisX[d*i:(i+1)*d, :] * coeff[i])
    return Y

def _poly_arr_cpp(coeffs, basisX, L, dev='cpu'):
    d = L.shape[0]
    Y = _FaustCorePy.polyCoeffs(d, basisX, coeffs)
    return Y

def _poly_Faust_cpp(coeffs, basisFaust, dev='cpu'):
    Y = basisFaust.m_faust.polyCoeffs(coeffs)
    return Y


def _chebyshev(L, K, T0, T1, rR, dev='cpu'):
    d = L.shape[0]
    factors = [T0]
    if(K > 0):
        factors.insert(0, T1)
        for i in range(2, K + 1):
            Ti = _chebyshev_Ti_matrix(rR, L, i)
            factors.insert(0, Ti)
    T = Faust(factors, dev=dev)
    return T  # K-th poly is T[K*L.shape[0]:,:]


def _chebyshev_gen(L, T0, T1, rR, dev='cpu'):
    if isFaust(T0):
        T = T0
    else:
        T = Faust(T0)
    yield T
    if isFaust(T1):
        T = T1 @ T
    else:
        T = Faust(T1) @ T
    yield T
    i = 2
    while True:
        Ti = _chebyshev_Ti_matrix(rR, L, i)
        if isFaust(Ti):
            T = Ti @ T
        else:
            T = Faust(Ti) @ T
        yield T
        i += 1


def _chebyshev_Ti_matrix(rR, L, i):
    d = L.shape[0]
    if i <= 2:
        R = rR
    else:
        #zero = csr_matrix((d, (i-2)*d), dtype=float)
        zero = _zeros_like(L, shape=(d, (i-2)*d))
        R = _hstack((zero, rR))
    di = d*i
    Ti = _vstack((_eyes_like(L, shape=di), R))
    return Ti


def _zeros_like(M, shape=None):
    """
    Returns a zero of the same type of M: csr_matrix, pyfaust.Faust.
    """
    if isinstance(shape, type(None)):
        shape = M.shape
    if isFaust(M):
        zero = csr_matrix(([0], ([0], [0])), shape=shape)
        return Faust(zero)
    elif isinstance(M, csr_matrix):
        zero = csr_matrix(shape, dtype=M.dtype)
        return zero
    else:
        raise TypeError('M must be a Faust or a scipy.sparse.csr_matrix.')


def _eyes_like(M, shape=None):
    """
    Returns an identity of the same type of M: csr_matrix, pyfaust.Faust.
    """
    if isinstance(shape, type(None)):
        shape = M.shape[1]
    if isFaust(M):
        return feye(shape)
    elif isinstance(M, csr_matrix):
        return sp.eye(shape, format='csr')
    else:
        raise TypeError('M must be a Faust or a scipy.sparse.csr_matrix.')


def _vstack(arrays):
    _arrays = _build_consistent_tuple(arrays)
    if isFaust(arrays[0]):
        # all arrays are of type Faust
        return fvstack(arrays)
    else:
        # all arrays are of type csr_matrix
        return sp.vstack(arrays, format='csr')


def _hstack(arrays):
    _arrays = _build_consistent_tuple(arrays)
    if isFaust(arrays[0]):
        # all arrays are of type Faust
        return fhstack(arrays)
    else:
        # all arrays are of type csr_matrix
        return sp.hstack(arrays, format='csr')


def _build_consistent_tuple(arrays):
    contains_a_Faust = False
    for a in arrays:
        if isFaust(a):
            contains_a_Faust = True
            break
    if contains_a_Faust:
        _arrays = []
        for a in arrays:
            if not isFaust(a):
                a = Faust(a)
            _arrays.append(a)
        return tuple(_arrays)
    else:
        return arrays

class FaustPoly(Faust):

    def __init__(self, *args, **kwargs):
        super(FaustPoly, self).__init__(*args, **kwargs)

    def _generator(self):
        F = self
        while True:
            F_next = FaustPoly(core_obj=F.m_faust.polyNext())
            F = F_next
            yield F

# experimental block end
