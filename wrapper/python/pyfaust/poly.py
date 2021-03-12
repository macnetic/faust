# experimental block start
# @PYFAUST_LICENSE_HEADER@
## @package pyfaust.poly @brief This module provides polynomials as Faust objects.

import scipy.sparse as sp
import _FaustCorePy
from scipy.sparse import csr_matrix
from pyfaust import (Faust, isFaust, eye as feye, vstack as fvstack, hstack as
                     fhstack)
from scipy.sparse.linalg import eigsh


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


def basis(L, K, basis_name, ret_gen=False, dev='cpu', T0=None):
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

    Returns:
        The Faust of the K+1 Chebyshev polynomials.
    """
    if basis_name.lower() == 'chebyshev':
        return Chebyshev(L, K, ret_gen=ret_gen, dev=dev, T0=T0)


def basis_cpp(L, K, basis_name, ret_gen=False, dev='cpu', T0=None):
    """

    """
    F = Faust(core_obj=_FaustCorePy.FaustCore.polyBasis(L, K))
    return F


def poly(coeffs, L=None, basis=Chebyshev, dev='cpu'):
    """
        Returns the linear combination of the polynomials defined by basis.

        Args:
            coeffs: the linear combination coefficients (numpy.array).
            basis: either the function to build the polynomials basis on L or the Faust of
            polynomials if already built externally.
            L: the symmetric matrix on which the polynomials are built, can't be None
            if basis is a function (not a Faust).
            dev: the device to instantiate the returned Faust ('cpu' or 'gpu').

        Returns:
            The linear combination Faust.
    """
    K = coeffs.size-1
    if isFaust(basis):
        F = basis
        if F.device != dev:
            F = F.clone(dev=dev)
    else:
        if L is None:
            raise ValueError('The L matrix must be set to build the'
                             ' polynomials.')
        F = poly(L, K, dev=dev)
    Id = sp.eye(L.shape[1], format="csr")
    scoeffs = sp.hstack(tuple(Id*coeffs[i] for i in range(0, K+1)),
                        format="csr")
    Fc = Faust(scoeffs, dev=dev) @ F
    return Fc


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
# experimental block end
