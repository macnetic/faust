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


def Chebyshev(L, K, ret_gen=False, dev='cpu', T0=None, impl="native"):
    """
    Builds the Faust of the Chebyshev polynomial basis defined on the sparse matrix L.

    Args:
        L: the sparse scipy matrix in CSR format (scipy.sparse.csr_matrix).
           L can aslo be a Faust if impl is "py".
        K: the degree of the last polynomial, i.e. the K+1 first polynomials are built.
        dev: the device to instantiate the returned Faust ('cpu' or 'gpu').
        'gpu' is not available yet for impl='native'.
        ret_gen: to return a generator of polynomials in addition to the
        polynomial itself (the generator starts from the the
        K+1-degree polynomial, and allows this way to compute the next
        polynomial simply with the instruction: next(generator)).
        T0: to define the 0-degree polynomial as something else than the identity.
        impl: "native" (by default) for the C++ impl., "py" for the Python impl.

    Returns:
        The Faust of the K+1 Chebyshev polynomials.

    Example:
        >>> from pyfaust.poly import Chebyshev
        >>> from scipy.sparse import random
        >>> L = random(50, 50, .02, format='csr')
        >>> L = L@L.T
        >>> K = 3
        >>> F = Chebyshev(L, K)
        >>> F
        Faust size 200x50, density 0.0654, nnz_sum 654, 4 factor(s):
        - FACTOR 0 (real) SPARSE, size 200x150, density 0.00893333, nnz 268
        - FACTOR 1 (real) SPARSE, size 150x100, density 0.0145333, nnz 218
        - FACTOR 2 (real) SPARSE, size 100x50, density 0.0236, nnz 118
        - FACTOR 3 (real) SPARSE, size 50x50, density 0.02, nnz 50

        >>> F, gen = Chebyshev(L, K, ret_gen=True)
        >>> F
        Faust size 200x50, density 0.0684, nnz_sum 684, 4 factor(s):
        - FACTOR 0 (real) SPARSE, size 200x150, density 0.00926667, nnz 278
        - FACTOR 1 (real) SPARSE, size 150x100, density 0.0152, nnz 228
        - FACTOR 2 (real) SPARSE, size 100x50, density 0.0256, nnz 128
        - FACTOR 3 (real) SPARSE, size 50x50, density 0.02, nnz 50

         Generate the next basis (the one with one additional dimension,
         whose the polynomial greatest degree is K+1 = 4)

        >>> G = next(gen)
        >>> G
        Faust size 250x50, density 0.08096, nnz_sum 1012, 5 factor(s):
        - FACTOR 0 (real) SPARSE, size 250x200, density 0.00656, nnz 328
        - FACTOR 1 (real) SPARSE, size 200x150, density 0.00926667, nnz 278
        - FACTOR 2 (real) SPARSE, size 150x100, density 0.0152, nnz 228
        - FACTOR 3 (real) SPARSE, size 100x50, density 0.0256, nnz 128
        - FACTOR 4 (real) SPARSE, size 50x50, density 0.02, nnz 50

        The factors 0 to 3 of G are views of the same factors of F.
        They are not duplicated in memory
    """
    if impl == "py":
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
    elif impl == "native":
        F = FaustPoly(core_obj=_FaustCorePy.FaustCore.polyBasis(L, K))
        if ret_gen:
            g = F._generator()
            return F, g
        else:
            return F
    else:
        raise ValueError(impl+" is an unknown implementation.")



def basis(L, K, basis_name, ret_gen=False, dev='cpu', T0=None, impl="native"):
    """
    Builds the Faust of the polynomial basis defined on the sparse matrix L.

    Args:
        L: the sparse scipy matrix in CSR format (scipy.sparse.csr_matrix).
           L can aslo be a Faust if impl is "py".
        K: the degree of the last polynomial, i.e. the K+1 first polynomials are built.
        basis_name: 'chebyshev', and others yet to come.
        dev: the device to instantiate the returned Faust ('cpu' or 'gpu').
        'gpu' is not available yet for impl='native'.
        ret_gen: to return a generator of polynomials in addition to the
        polynomial itself (the generator starts from the the
        K+1-degree polynomial, and allows this way to compute the next
        polynomial simply with the instruction: next(generator)).
        impl: "native" (by default) for the C++ impl., "py" for the Python impl.

    Returns:
        The Faust of the basis composed of the K+1 orthogonal polynomials.

    Example:
        >>> from pyfaust.poly import basis
        >>> from scipy.sparse import random
        >>> L = random(50, 50, .02, format='csr')
        >>> L = L@L.T
        >>> K = 3
        >>> F = basis(L, K, 'chebyshev')
        >>> F
        Faust size 200x50, density 0.0654, nnz_sum 654, 4 factor(s):
        - FACTOR 0 (real) SPARSE, size 200x150, density 0.00893333, nnz 268
        - FACTOR 1 (real) SPARSE, size 150x100, density 0.0145333, nnz 218
        - FACTOR 2 (real) SPARSE, size 100x50, density 0.0236, nnz 118
        - FACTOR 3 (real) SPARSE, size 50x50, density 0.02, nnz 50

        >>> F, gen = basis(L, K, 'chebyshev', ret_gen=True)
        >>> F
        Faust size 200x50, density 0.0684, nnz_sum 684, 4 factor(s):
        - FACTOR 0 (real) SPARSE, size 200x150, density 0.00926667, nnz 278
        - FACTOR 1 (real) SPARSE, size 150x100, density 0.0152, nnz 228
        - FACTOR 2 (real) SPARSE, size 100x50, density 0.0256, nnz 128
        - FACTOR 3 (real) SPARSE, size 50x50, density 0.02, nnz 50

         Generate the next basis (the one with one additional dimension,
         whose the polynomial greatest degree is K+1 = 4)

        >>> G = next(gen)
        >>> G
        Faust size 250x50, density 0.08096, nnz_sum 1012, 5 factor(s):
        - FACTOR 0 (real) SPARSE, size 250x200, density 0.00656, nnz 328
        - FACTOR 1 (real) SPARSE, size 200x150, density 0.00926667, nnz 278
        - FACTOR 2 (real) SPARSE, size 150x100, density 0.0152, nnz 228
        - FACTOR 3 (real) SPARSE, size 100x50, density 0.0256, nnz 128
        - FACTOR 4 (real) SPARSE, size 50x50, density 0.02, nnz 50

        The factors 0 to 3 of G are views of the same factors of F.
        They are not duplicated in memory.

    """
    if basis_name.lower() == 'chebyshev':
        return Chebyshev(L, K, ret_gen=ret_gen, dev=dev, T0=T0, impl=impl)
    else:
        raise ValueError(basis_name+" is not a valid basis name")


def poly(coeffs, basis='chebyshev', L=None, dev='cpu', impl='native'):
    """
        Computes the linear combination of the polynomials defined by basis.

        Args:
            coeffs: the linear combination coefficients (vector as a numpy.ndarray).
            basis: either the name of the polynomial basis to build on L or the
            basis if already built externally (as a FaustPoly or an equivalent
            np.ndarray).
            L: the sparse scipy matrix in CSR format (scipy.sparse.csr_matrix).
            L can aslo be a Faust if impl is "py". It can't be None if basis is not a FaustPoly or a numpy.ndarray.
            dev: the device to instantiate the returned Faust ('cpu' or 'gpu').
            'gpu' is not available yet for impl='native'.

        Returns:
            The linear combination Faust or np.ndarray depending on if basis is itself a Faust or a np.ndarray.

        Example:
            >>> import numpy as np
            >>> from pyfaust.poly import basis, poly
            >>> from scipy.sparse import random
            >>> K = 3
            >>> L = random(50, 50, .02, format='csr')
            >>> L = L@L.T
            >>> F = basis(L, K, 'chebyshev')
            >>> coeffs = np.array([.5, 1, 2, 3])
            >>> G = poly(coeffs, F)
            >>> G
            Faust size 50x50, density 0.3608, nnz_sum 902, 5 factor(s):
            - FACTOR 0 (real) SPARSE, size 50x200, density 0.02, nnz 200
            - FACTOR 1 (real) SPARSE, size 200x150, density 0.00946667, nnz 284
            - FACTOR 2 (real) SPARSE, size 150x100, density 0.0156, nnz 234
            - FACTOR 3 (real) SPARSE, size 100x50, density 0.0268, nnz 134
            - FACTOR 4 (real) SPARSE, size 50x50, density 0.02, nnz 50

            Above G is a Faust because F is too.
            Below the full array of the Faust F is passed, so an array is returned into GA.
            >>> GA = poly(coeffs, F.toarray())
            >>> type(GA)
            numpy.ndarray

            But of course they are equal:

            >>> np.allclose(GA, G.toarray())
            True

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
        raise TypeError('basis is neither a str neither a Faust nor'
                        ' a numpy.ndarray')
    else:
        F = basis
    if L == None:
        d = F.shape[1]
    else:
        d = L.shape[0]
    if impl == 'py':
        if isFaust(F):
            Id = sp.eye(d, format="csr")
            scoeffs = sp.hstack(tuple(Id*coeffs[i] for i in range(0, K+1)),
                                format="csr")
            Fc = Faust(scoeffs, dev=dev) @ F
            return Fc
        else:
           # F is a np.ndarray
           return _poly_arr_py(coeffs, F, d, dev=dev)
    elif impl == 'native':
        if isFaust(F):
            Fc = _poly_Faust_cpp(coeffs, F)
            if F.device != dev:
                Fc = Fc.clone(dev=dev)
            return Fc
        else:
            return _poly_arr_cpp(coeffs, F, d, dev='cpu')
    else:
        raise ValueError(impl+" is an unknown implementation.")

def _poly_arr_py(coeffs, basisX, d, dev='cpu'):
    """
    """
    mt = True # multithreading
    n = basisX.shape[1]
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

def _poly_arr_cpp(coeffs, basisX, d, dev='cpu'):
    Y = _FaustCorePy.polyCoeffs(d, basisX, coeffs)
    return Y

def _poly_Faust_cpp(coeffs, basisFaust, dev='cpu'):
    Y = Faust(core_obj=basisFaust.m_faust.polyCoeffs(coeffs))
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
    """
    Subclass of Faust specialized for orthogonal polynomial basis.

    This class is used only for the native implementation of the poly functions.

    NOTE: it is not advisable to use this class directly.

    """
    def __init__(self, *args, **kwargs):
        super(FaustPoly, self).__init__(*args, **kwargs)

    def _generator(self):
        F = self
        while True:
            F_next = FaustPoly(core_obj=F.m_faust.polyNext())
            F = F_next
            yield F

# experimental block end
