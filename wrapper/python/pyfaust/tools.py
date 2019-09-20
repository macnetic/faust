# -*- coding: utf-8 -*-
# @PYFAUST_LICENSE_HEADER@

import numpy as np
from scipy.sparse.linalg import spsolve_triangular
from numpy.linalg import solve, lstsq, norm
from scipy.sparse import hstack, vstack, csr_matrix
from numpy import concatenate as cat
from numpy import matrix, zeros, argmax, empty
from pyfaust import Faust

def greed_omp_chol(x, A, m, stopTol=None, verbose=False):
    """
    Runs greedy OMP algorithm optimized by Cholesky decomposition.

    Params:
        A: the dictionary as a numpy matrix or a Faust.
        x: the vector to approximate by A*y.
        m: the dimension from which A maps (the number of columns of A).
        stopTol:  stopping criterion based on the number of iterations (by default it's a
        quarter of the size of x) and it can't be greater than the size of x. This
        criterion doesn't disable the other criterion based on the precision of
        the solution (the algorithm will stop anyway if a precision lower or equal to
        1^-16 is reached).
        verbose: to enable the verbosity (value to True).

    Returns:
        y: the solution of x = A*y.
    """
    #TODO: check x is a numpy.matrix (or a matrix_csr ?)
    sp = x.shape
    if(sp[1] == 1):
        n = sp[0]
    elif sp[0] == 1:
        x = x.H
        n = sp[1]
    else:
        raise Exception("x must be a vector")
    if(Faust.isFaust(A) or isinstance(A, matrix)):
        P = lambda z : matrix(A*z, copy=False)
        Pt = lambda z : matrix(A.H*z, copy=False)
    else:
        raise Exception("A must be a Faust or a numpy.matrix. Here A is "
                        "a:"+str(type(A)))
    #TODO: raise exception for x also if not a matrix-vector (type matrix)

    # pre-calculate Pt(x) because it's used repeatedly
    Ptx = Pt(x)

    maxM = stopTol
    # try if enough memory
    try:
        R = matrix(zeros(maxM))
    except:
        raise Exception("Variable size is too large.")

    if(not stopTol):
        stopTol = int(n/4)

    r_count = 0
    sigsize = x.H*x/n
    # initialize
    s_initial = matrix(zeros((m,1))).astype(np.complex)
    residual = x
    s = s_initial
    R = matrix(empty((stopTol+1,stopTol+1))).astype(np.complex)
    oldErr = sigsize
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

        Rs = lstsq(R[0:r_count, 0:r_count].H, Ptx[IN])[0]
        s[IN] = lstsq(R[0:r_count, 0:r_count], Rs)[0]

        residual = x - P(s)
        DR = Pt(residual)

        err = residual.H*residual/n

        done = it >= stopTol
        if(not done and verbose):
            print("Iteration",  it, "---", stopTol-it, "iterations to go")

        if(err<1**-16 and r_count > 1):
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
        new_col = lstsq(R.H, Pt_new_vector[index[0:-1]])[0] # solve() only works
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


