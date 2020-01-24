# -*- coding: utf-8 -*-
# @PYFAUST_LICENSE_HEADER@
## @package pyfaust.proj @brief This module provides matrix projectors.

from pyfaust.factparams import *
if sys.version_info > (3,0):
    from abc import ABC, abstractmethod
else:
    from abc import abstractmethod
    ABC = object # trick to handle py2 missing ABC
                 # but not using abstract class in py2.7

class proj_gen(ABC):
    """
    The parent abstract class to represent projectors (as functors).
    """
    @abstractmethod
    def __init__(self):
        pass

    def __call__(self, M):
        return self.constraint.project(M)

class toeplitz(proj_gen):
    """
    Functor for the TOEPLITZ projector.
    """
    def __init__(self, shape, normalized=True, pos=False):
        self.constraint = ConstraintMat('toeplitz', np.empty(shape), normalized, pos)

#    def __call__(self, M):
#        P = np.empty(M.shape)
#        for i in range(M.shape[0]):
#            m = np.mean(np.diag(M,i))
#            m_ = np.mean(np.diag(M,-i))
#            dl = len(np.diag(M,i))
#            I =  list(range(0,dl))
#            J =  [ i+j for j in range(0,dl) ]
#            P[I,J] = m
#            P[J,I] = m_
#        return P

class circ(proj_gen):
    """
    Functor for the CIRC(ulant) projector.
    """
    def __init__(self, shape, normalized=True, pos=False):
        self.constraint = ConstraintMat('circ', np.empty(shape), normalized, pos)

#    def __call_(self, M):
#        P = np.empty(M.shape)
#        for i in range(M.shape[0]):
#            j = i - M.shape[0]
#            m = np.mean(np.hstack((np.diag(M,i), np.diag(M,j))))
#            print("m=", m, "mi=", np.mean(np.diag(M,i)), "mj=",
#                                             np.mean(np.diag(M,j)), "i=", i,
#                  "j=", j)
#            dli = M.shape[0]-i
#            dlj = M.shape[0]+j
#            I =  list(range(dli))
#            J = [ k+i for k in range(dli) ]
#            P[I,J] = m
#            I = [ k+-j for k in range(dlj) ]
#            J = [ k for k in range(dlj) ]
#            P[I,J] = m
#        return P

class hankel(proj_gen):
    """
    Functor for the HANKEL projector.
    """
    def __init__(self, shape, normalized=True, pos=False):
        self.constraint = ConstraintMat('hankel', np.empty(shape), normalized, pos)


class sp(proj_gen):
    """
    Functor for the SP projector. A, the image matrix, is such that \f$ \| A \|_0 = k,  \| A\|_F = 1\f$.
    """

    def __init__(self, shape, k, normalized=True, pos=False):
        """

        Args:
            k: the number of nonzeros of the projection image.
            normalized: True to normalize the projection image according to its 2-norm.
            pos: True to skip negative values (replaced by zero) of the matrix to project.

        """
        self.constraint = ConstraintInt('sp', shape[0], shape[1], k, normalized, pos)

class splin(proj_gen):
    """
    Functor for the SPLIN projector. A, the image matrix, is defined by \f$ \forall i \in \{0,...,shape[0]-1\} \| \f$ the i-th row \f$ A_{i,*}\f$ is such that \f$ \| A_{i,*}\|_0 = k,  \| A\|_F = 1\f$.
    """
    def __init__(self, shape, k, normalized=True, pos=False):
        """
        Args:
            k: the number of nonzeros of the projection image.
            normalized: True to normalize the projection image according to its 2-norm.
            pos: True to skip negative values (replaced by zero) of the matrix to project.

        """
        self.constraint = ConstraintInt('splin', shape[0], shape[1], k, normalized, pos)

class spcol(proj_gen):
    """
    Functor for the SPCOL projector. A, the image matrix, is defined by \f$ \forall j \in \{0,...,shape[1]-1\} \f$ the j-th column \f$ A_{*,j}\f$ is such that \f$ \| A_{*,j}\|_0 = k,  \| A\|_F = 1\f$.
    """
    def __init__(self, shape, k, normalized=True, pos=False):
        """

        Args:
            S: the support matrix.
            normalized: True to normalize the projection image according to its 2-norm.
            pos: True to skip negative values (replaced by zero) of the matrix to project.

        """
        self.constraint = ConstraintInt('spcol', shape[0], shape[1], k, normalized, pos)

class splincol(proj_gen):
    """
    Functor for the SPLINCOL projector.

    It's the union of SPLIN and SPCOL projectors.
    """
    def __init__(self, shape, k, normalized=True, pos=False):
        """

        Args:
            S: the support matrix.
            normalized: True to normalize the projection image according to its 2-norm.
            pos: True to skip negative values (replaced by zero) of the matrix to project.

        """
        self.constraint = ConstraintInt('splincol', shape[0], shape[1], k, normalized, pos)

class supp(proj_gen):
    """
        Functor for the SUPP projector. A, the image matrix, is such that np.nonzero(A) == np.nonzero(S).
    """
    def __init__(self, S, normalized=True, pos=False):
        """

        Args:
            S: the support matrix.
            normalized: True to normalize the projection image according to its 2-norm.
            pos: True to skip negative values (replaced by zero) of the matrix to project.
        """
        self.constraint = ConstraintMat('supp', S, normalized, pos)

class const(proj_gen):
    """
        Functor for the CONST projector. A, the image matrix, is such that A == C.
    """
    def __init__(self, C, normalized=False):
        """

        Args:
            C: the constant matrix.
            normalized: True to normalize the projection image according to its 2-norm.
            pos: True to skip negative values (replaced by zero) of the matrix to project.

        """
        self.constraint = ConstraintMat('const', C, normalized, pos=False)

class normcol(proj_gen):
    """
        Functor for the NORMCOL projector. A, the image matrix, is defined by \f$ \forall j \in \{0,...,shape[1]-1\} \f$ the j-th column \f$ A_{*,j} \f$ is such that \f$\| A_{*,j}\|_2 = s  \f$.

    """
    def __init__(self, shape, s, normalized=False, pos=False):
        """
        Args:
            s: the column 2-norm.
            normalized: True to normalize the projection image according to its 2-norm.
            pos: True to skip negative values (replaced by zero) of the matrix to project.

        """
        self.constraint = ConstraintReal('normcol', shape[0], shape[1], s, normalized, pos)

class normlin(proj_gen):
    """
        Functor for the NORMLIN projector. A, the image matrix, is defined by \f$ \forall i \in \{0,...,shape[0]-1\}\f$ the i-th row \f$ A_{i,*} \f$ is such that \f$ \| A_{i,*} \|_2 = s \f$.

    """

    def __init__(self, shape, s, normalized=False, pos=False):
        """
        Args:
            s: the row 2-norm.
            normalized: True to normalize the projection image according to its 2-norm.
            pos: True to skip negative values (replaced by zero) of the matrix to project.

        """
        self.constraint = ConstraintReal('normlin', shape[0], shape[1], s, normalized, pos)


