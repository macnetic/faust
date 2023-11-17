# -*- coding: utf-8 -*-
# @PYFAUST_LICENSE_HEADER@

## @package pyfaust.tools @brief The pyfaust tools module

import numpy as np
from scipy.sparse.linalg import spsolve_triangular
from numpy.linalg import solve, lstsq
from scipy.sparse import hstack, vstack, csr_matrix
from numpy import concatenate as cat
from numpy import zeros, argmax, empty, ndarray


def omp(y, D, maxiter=None, tol=0, relerr=True, verbose=False):
    """
    Runs the greedy OMP algorithm optimized by Cholesky decomposition.

    Args:
        y: (np.ndarray)
            the vector to approximate by D@x.
        D: (np.ndarray or a Faust)
            the dictionary as a numpy array or a Faust.
        maxiter: (int or NoneType)
            the maximum number of iterations of the algorithm.
            By default (None) it's y's dimension: max(y.shape).
        tol: (float)
            the tolerance error to reach for the algorithm to stop. By default,
            it's zero for not stopping on error criterion.
        relerr: (bool)
            the type of error stopping criterion. Default to
            True to use relative error, otherwise (False) the absolute error is used.
        verbose: (bool)
            to enable the verbosity (True).

    Returns:
        x: the solution of y = D@x (according to the error).

    Example:
        >>> from pyfaust.tools import omp
        >>> from scipy.sparse import random
        >>> from pyfaust import rand, seed
        >>> import numpy as np
        >>> np.random.seed(42) # just for reproducibility
        >>> seed(42)
        >>> D = rand(1024, 1024)
        >>> x0 = random(1024, 1, .01)
        >>> y = D @ x0
        >>> maxiter = 17
        >>> x = omp(y, D, maxiter, tol=10**-16)
        Stopping. Exact signal representation found!
        >>> # omp() runs at most maxiter iterations until the error tolerance is
        >>> # reached

    """
    from pyfaust import Faust
    # check y is a numpy.ndarray (or a matrix_csr ?)
    m = D.shape[1]
    sp = y.shape
    if(sp[1] == 1):
        n = sp[0]
    elif sp[0] == 1:
        y = y.T.conj()
        n = sp[1]
    else:
        raise Exception("y must be a vector")
    if(Faust.isFaust(D) or isinstance(D, ndarray)):
        P = lambda z: D@z
        Pt = lambda z: D.T.conj()@z
    else:
        raise Exception("D must be a Faust or a numpy.ndarray. Here D is "
                        "a:"+str(type(D)))

    # pre-calculate Pt(y) because it's used repeatedly
    Ptx = Pt(y)

    if(not maxiter):
        maxiter = int(n)
    else:
        if(not (type(maxiter) in [int, float])):
            raise ValueError("maxiter must be a number.")

    if(relerr):
        tolerr = tol * (y.T.conj()@y)[0,0]
    else:
        tolerr = tol**2

    # try if enough memory
    # TODO: maybe it should be done directly on empty()
    # below
    try:
        R = zeros(maxiter)
    except:
        raise Exception("Variable size is too large.")

    r_count = 0
    # initialize
    s_initial = zeros((m,1)).astype(np.complex128)
    residual = y
    s = s_initial
    R = empty((maxiter+1,maxiter+1)).astype(np.complex128)
    # TODO: why do we use complex, we should do in function of D, y
    # likewise H should be T in real case
    oldErr = y.T.conj()@y
    # err_mse = []

    t = 0
    IN = []
    DR = Pt(residual)
    done = False
    it = 1

    MAXITER = n

    while not done:

        # select new element
        DR[IN] = 0
        I = argmax(abs(DR))
        IN += [I]
        # update R
        R[0:r_count+1, 0:r_count+1] = UpdateCholesky(R[0:r_count,
                                                       0:r_count], P, Pt,
                                                     IN, m)

        r_count+=1

        Rs = lstsq(R[0:r_count, 0:r_count].T.conj(), Ptx[IN], rcond=-1)[0]
        s[IN] = lstsq(R[0:r_count, 0:r_count], Rs, rcond=-1)[0]

        residual = y - P(s)
        DR = Pt(residual)

        err = residual.T.conj()@residual #/n

        done = it >= maxiter
        if(not done and verbose):
            print("Iteration",  it, "---", maxiter-it, "iterations to go")

        if(err < tolerr and r_count > 1):
            print("Stopping. Exact signal representation found!")
            done=True

        if(it >= MAXITER):
            print("Stopping. Maximum number of iterations reached!")
            done = True
        else:
            it += 1
            oldErr = err

    if((s.imag == 0).all()):
        return s.real
    return s


def UpdateCholesky(R,P,Pt,index,m):
    """
    Args:
        R: (scipy.sparse.csr_matrix or a np.ndarray)
            A triangular matrix
        P: (callable)
            A function to implement the matrix product.
        Pt: (callable)
            A function to implement the matrix transpose product.
        index: (list[int])
            a list of all non-zero indices of length R.shape[1]+1
        m: (int)
            dimension of space from which P maps.

    Returns:
        R is a triangular matrix such that R'*R = Pt(index)*P(index)
    """

    if(isinstance(R,csr_matrix)):
        return UpdateCholeskySparse(R,P,Pt,index,m)
    else:
        return UpdateCholeskyFull(R,P,Pt,index,m)

def UpdateCholeskyFull(R,P,Pt,index,m):
    """
    Same as :py:func:`.UpdateCholesky` but R is a np.ndarray.
    """
    li = len(index)

    if(li != R.shape[1]+1):
        raise Exception("Incorrect index length or size of R.")


    mask = zeros((m,1))
    mask[index[-1]] = 1
    new_vector = P(mask)


    if(li == 1):
        R = np.sqrt(new_vector.T.conj()@new_vector.astype(np.complex128))
    else:
        Pt_new_vector = Pt(new_vector)
        # linsolve_options_transpose.UT = true;
        # linsolve_options_transpose.TRANSA = true;
        # matlab opts for linsolve() are respected here
        new_col = lstsq(R.T.conj(), Pt_new_vector[index[0:-1]], rcond=-1)[0] # solve() only works
        # for full rank square matrices, that's why we use ltsq
        R_ii = np.sqrt(new_vector.T.conj()@new_vector -
                       new_col.T.conj()@new_col.astype(np.complex128))
        R = cat((
                cat((R,new_col),axis=1),
                cat((zeros((1, R.shape[1])), R_ii),axis=1)
                ),axis=0)
    #assert(np.allclose((R.T.conj()@R).T.conj(), R.T.conj()@R))
    #D = P(np.eye(m,m))
    #assert(np.allclose(D[:,index].T.conj()@D[:,index], R.T.conj()@R))
    return R

def UpdateCholeskySparse(R,P,Pt,index,m):
    """
    Same as :py:func:`.UpdateCholesky` but R is a scipy.sparse.csr_matrix.
    """

    li = len(index)

    if(li != R.shape[1]+1):
        raise Exception("Incorrect index length or size of R.")


    mask = csr_matrix((m,1)) #np.zeros((m,1))
    mask[index[-1]] = 1
    new_vector = P(mask)

    if(li == 1):
        R = np.sqrt(new_vector.T.conj()@new_vector.astype(np.complex128))
    else:
        Pt_new_vector = Pt(new_vector)
        # linsolve_options_transpose.UT = true;
        # linsolve_options_transpose.TRANSA = true;
        # matlab opts for linsolve() are respected here
        new_col = spsolve_triangular(R.T.conj(), Pt_new_vector[index[0:-1]])
        R_ii = np.sqrt(new_vector.T.conj()@new_vector -
                       new_col.T.conj()@new_col.astype(np.complex128))
        R = vstack((hstack((R, new_col)), hstack((csr_matrix((1, R.shape[1])),
                                                 R_ii))))

    return R

greed_omp_chol = omp

def _sanitize_dtype(dtype):
    """
    Verifies the dtype is pyfaust-compatible and returns it as str.

    Returns:
        one of the str of the dtype (float32, float64, complex).
    """
    if dtype in [np.float64, np.float64, np.double, float, 'float',
                 'double', 'float64']:
        return 'float64'
    elif dtype in [np.float32, 'float32']:
        return 'float32'
    elif dtype in [np.complex128, np.complex128, 'complex', 'complex128']:
        return 'complex'
    else:
        raise TypeError(str(dtype)+' is not a dtype compatible with pyfaust'
                        ' (float32, float64/double/float, complex128/complex)')
