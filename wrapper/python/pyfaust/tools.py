# -*- coding: utf-8 -*-
# @PYFAUST_LICENSE_HEADER@

## @package pyfaust.tools @brief The pyfaust tools module

import numpy as np
from scipy.sparse.linalg import spsolve_triangular
from numpy.linalg import solve, lstsq, norm
from scipy.sparse import hstack, vstack, csr_matrix
from numpy import concatenate as cat
from numpy import matrix, zeros, argmax, empty
from pyfaust import Faust

def omp(y, D, maxiter=None, tol=0, relerr=True, verbose=False):
    """
    Runs the greedy OMP algorithm optimized by Cholesky decomposition.

    Args:
        D: the dictionary as a numpy matrix or a Faust.
        y: the vector to approximate by D*x.
        maxiter: the maximum number of iterations of the algorithm.
        By default (None) it's y's dimension: max(y.shape).
        tol: the tolerance error under what the algorithm stops. By default,
        it's zero for not stopping on error criterion.
        relerr: the type of error stopping criterion. Default to
        True to use relative error, otherwise (False) the absolute error is used.
        verbose: to enable the verbosity (value to True).

    Returns:
        x: the solution of y = D*x (according to the error).
    """
    # check y is a numpy.matrix (or a matrix_csr ?)
    if(isinstance(y, np.ndarray)): y = matrix(y, copy=False)
    m = D.shape[1]
    sp = y.shape
    if(sp[1] == 1):
        n = sp[0]
    elif sp[0] == 1:
        y = y.H
        n = sp[1]
    else:
        raise Exception("y must be a vector")
    if(Faust.isFaust(D) or isinstance(D, matrix)):
        P = lambda z : matrix(D*z, copy=False)
        Pt = lambda z : matrix(D.H*z, copy=False)
    else:
        raise Exception("D must be a Faust or a numpy.matrix. Here D is "
                        "a:"+str(type(D)))

    # pre-calculate Pt(y) because it's used repeatedly
    Ptx = Pt(y)

    if(not maxiter):
        maxiter = int(n)
    else:
        if(not (type(maxiter) in [int, float])):
            raise ValueError("maxiter must be a number.")

    if(relerr):
        tolerr = tol * (y.H*y)[0,0]
    else:
        tolerr = tol**2

    # try if enough memory
    try:
        R = matrix(zeros(maxiter))
    except:
        raise Exception("Variable size is too large.")

    r_count = 0
    # initialize
    s_initial = matrix(zeros((m,1))).astype(np.complex)
    residual = y
    s = s_initial
    R = matrix(empty((maxiter+1,maxiter+1))).astype(np.complex)
    oldErr = y.H*y
    err_mse = []

    t=0
    IN = []
    DR = Pt(residual)
    done = False
    it = 1

    MAXITER=n

    while not done:

        # select new element
        DR[IN] = 0
        I = argmax(abs(DR))
        IN += [I]
        # update R
        R[0:r_count+1, 0:r_count+1] = UpdateCholeskyFull(R[0:r_count,
                                                           0:r_count], P, Pt,
                                                         IN, m)

        r_count+=1

        Rs = lstsq(R[0:r_count, 0:r_count].H, Ptx[IN], rcond=-1)[0]
        s[IN] = lstsq(R[0:r_count, 0:r_count], Rs, rcond=-1)[0]

        residual = y - P(s)
        DR = Pt(residual)

        err = residual.H*residual #/n

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
        R: must be a csr_matrix or a numpy.matrix
        P: a function
        Pt: a function
        index: a list of all non-zero indices of length R.shape[1]+1
        m: dimension of space from which P maps.

    Returns:
        R is a triangular matrix such that R'*R = Pt(index)*P(index)
    """

    if(isinstance(R,csr_matrix)):
        return UpdateCholeskySparse(R,P,Pt,index,m)
    else:
        return UpdateCholeskyFull(R,P,Pt,index,m)

def UpdateCholeskyFull(R,P,Pt,index,m):
    """
    Args:
        R: must be a numpy.matrix
        P: a function returning a numpy.matrix
        Pt: a function returning a numpy.matrix
        index: a list of all non-zero indices of length R.shape[1]+1
        m: dimension of space from which P maps.

    Returns:
        R is a triangular matrix such that R'*R = Pt(index)*P(index)
    """
    li = len(index)

    if(li != R.shape[1]+1):
        raise Exception("Incorrect index length or size of R.")


    mask = matrix(np.zeros((m,1)))
    mask[index[-1]] = 1
    new_vector = P(mask)


    if(li == 1):
        R = np.sqrt(new_vector.H*new_vector.astype(np.complex))
    else:
        Pt_new_vector = Pt(new_vector)
        # linsolve_options_transpose.UT = true;
        # linsolve_options_transpose.TRANSA = true;
        # matlab opts for linsolve() are respected here
        new_col = lstsq(R.H, Pt_new_vector[index[0:-1]], rcond=-1)[0] # solve() only works
        # for full rank square matrices, that's why we use ltsq
        R_ii = np.sqrt(new_vector.H*new_vector -
                       new_col.H*new_col.astype(np.complex))
        R = cat((
                cat((R,new_col),axis=1),
                cat((np.zeros((1, R.shape[1])), R_ii),axis=1)
                ),axis=0)
    #assert(np.allclose((R.H*R).H, R.H*R))
    #D = P(np.eye(m,m))
    #assert(np.allclose(D[:,index].H*D[:,index], R.H*R))
    return R

def UpdateCholeskySparse(R,P,Pt,index,m):
    """
    Args:
        R: must be a csr_matrix.
        P: a function
        Pt: a function
        index: a list of all non-zero indices of length R.shape[1]+1
        m: dimension of space from which P maps.

    Returns:
        R is a triangular matrix such that R'*R = Pt(index)*P(index)
    """

    li = len(index)

    if(li != R.shape[1]+1):
        raise Exception("Incorrect index length or size of R.")


    mask = csr_matrix((m,1)) #np.zeros((m,1))
    mask[index[-1]] = 1
    new_vector = P(mask)

    if(li == 1):
        R = np.sqrt(new_vector.H*new_vector.astype(np.complex))
    else:
        Pt_new_vector = Pt(new_vector)
        # linsolve_options_transpose.UT = true;
        # linsolve_options_transpose.TRANSA = true;
        # matlab opts for linsolve() are respected here
        new_col = spsolve_triangular(R.H, Pt_new_vector[index[0:-1]])
        R_ii = np.sqrt(new_vector.H*new_vector -
                       new_col.H*new_col.astype(np.complex))
        R = vstack((hstack((R, new_col)), hstack((csr_matrix((1, R.shape[1])),
                                                 R_ii))))

    return R

greed_omp_chol = omp
