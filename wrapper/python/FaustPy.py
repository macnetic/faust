# -*- coding: utf-8 -*-
##############################################################################
##                              Description:                                ##
##                                                                          ##
##          FaustPy is a python module  which delivered                     ##
##          a class named Faust which represents a dense matrix             ##
##          by a product of 'sparse' factors (i.e Faust)                    ##
##          Python wrapper class implemented in C++                         ##
##                                                                          ##
##  For more information on the FAuST Project, please visit the website     ##
##  of the project : <http://faust.inria.fr>                         ##
##                                                                          ##
##                              License:                                    ##
##  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      ##
##                      Luc Le Magoarou, Remi Gribonval                     ##
##                      INRIA Rennes, FRANCE                                ##
##                      http://www.inria.fr/                                ##
##                                                                          ##
##  The FAuST Toolbox is distributed under the terms of the GNU Affero      ##
##  General Public License.                                                 ##
##  This program is free software: you can redistribute it and/or modify    ##
##  it under the terms of the GNU Affero General Public License as          ##
##  published by the Free Software Foundation.                              ##
##                                                                          ##
##  This program is distributed in the hope that it will be useful, but     ##
##  WITHOUT ANY WARRANTY; without even the implied warranty of              ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    ##
##  See the GNU Affero General Public License for more details.             ##
##                                                                          ##
##  You should have received a copy of the GNU Affero General Public        ##
##  License along with this program.                                        ##
##  If not, see <http://www.gnu.org/licenses/>.                             ##
##                                                                          ##
##                             Contacts:                                    ##
##      Nicolas Bellot  : nicolas.bellot@inria.fr                           ##
##      Adrien Leman    : adrien.leman@inria.fr                             ##
##      Thomas Gautrais : thomas.gautrais@inria.fr                          ##
##      Luc Le Magoarou : luc.le-magoarou@inria.fr                          ##
##      Remi Gribonval  : remi.gribonval@inria.fr                           ##
##                                                                          ##
##############################################################################

import copy

import numpy as np
from scipy.io import savemat, loadmat
import FaustCorePy


class Faust:
    """FAµST class.
     This class represents a given dense matrix by a product of sparse matrices
     (i.e Faust).
     The main goal of Faust representation is to speed up operations on that
     matrix, especially the multiplication. Besides the time optimization, a Faust
     can reduce the memory space size needed both for storage and loading.

     Although the sparse matrices are more interesting for optimization it's not
     forbidden to define a Faust as a product of dense matrices or a mix up of dense
     and sparse matrices.

     The matrices composing the Faust product, also called the factors, are
     defined on complex or real fields. Hence a Faust can be a complex Faust or a
     real Faust.
    """

    def  __init__(F, factors=None, alpha=1.0, filepath=None, core_obj=None):
        """ Creates a Faust from a list of factors or alternatively from a file.

            Another easy way to create a Faust is to call the static method Faust.randFaust().

        Args:
            factors: list/tuple of numpy matrices.<br/>
                     The factors must respect the dimensions needed for
                     the product to be defined (for i=0 to len(factors)-1,
                     factors[i].shape[1] == factors[i+1].shape[0]).<br/>
                     The factors can be sparse or dense matrices
                     (either scipy.sparse.csr.csr_matrix or numpy.ndarray).
            filepath: the file from which a Faust is created.<br/>
                      The format is Matlab version 5 (.mat extension).<br/>
                      The file must have been saved before with Faust.save().
            alpha: multiplicative scalar applied to the factor product before
            to set the Faust with.
            core_obj: for internal purpose only. Please don't fill this argument.

        WARNING: filepath and factors arguments are multually exclusive. If you
        use filepath you must explicitely set argument with the keyword.

        Examples:
            >>> from FaustPy import Faust
            >>> import numpy as np
            >>> from scipy import sparse
            >>> factors = []
            >>> is_sparse = False
            >>> for i in range(0,5):
            >>>     if(is_sparse):
            >>>         factors += [ sparse.random(100,100, dtype=np.float64, format='csr',
            >>>                                   density=0.1)]
            >>>     else:
            >>>         factors += [ np.random.rand(100, 100).astype(np.float64) ]
            >>>     is_sparse = not is_sparse

            >>> # define a Faust with those factors
            >>> F = Faust(factors)

            >>> alpha = 2
            >>> G = Faust(factors, alpha) # G == alpha*F

            >>> F.save("F.mat")
            >>> # define a Faust from file
            >>> H = Faust(filepath="F.mat")
            >>> I = Faust(filepath="F.mat", alpha) # I == alpha*H

        <b/> See also Faust.save, Faust.randFaust

        """
        if(core_obj):
            F.m_faust = core_obj
        else:
            if(filepath and isinstance(filepath, str)):
                    contents = loadmat(filepath)
                    factors = contents['faust_factors'][0]
            if(factors is not None):
                F.m_faust = FaustCorePy.FaustCore(factors, alpha);
            #else:
                #TODO: manage empty Faust

    @staticmethod
    def randFaust(faust_type, field, min_num_factors, max_num_factors, min_dim_size,
                   max_dim_size, density=0.1):
        """
        Generates a random Faust.

            Args:
                faust_type: must be one of RandFaustType.DENSE, RandFaustType.SPARSE or RandFaustType.MIXTE.
                field: must be one of RandFaustType.REAL or RandFaustType.COMPLEX.
                min_num_factors: the minimal number of factors.
                max_num_factors: the maximal number of factors.
                min_dim_size: the minimal size of column and row dimensions.
                max_dim_size: the maximal size of column and row dimensions.
                density: the approximate density of factors.

        Returns:
            the random generated Faust.

        Examples:
            >>> from FaustPy import Faust, RandFaustType
            >>> F = Faust.randFaust(RandFaustType.MIXTE, RandFaustType.COMPLEX, 2, 3, 10, 20,.5)
            """
        if(not isinstance(faust_type, int) or faust_type not in [RandFaustType.SPARSE,
                                                       RandFaustType.DENSE,
                                                       RandFaustType.MIXTE]):
            raise ValueError("Faust.randFaust(): first arg must be an"
                             " integer equal to RandFaustType.SPARSE,"
                             " RandFaustType.DENSE or RandRandFaustType.MIXTE.")

        if(not isinstance(field, int) or field not in
           [RandFaustType.REAL, RandFaustType.COMPLEX]):
            raise ValueError("Faust.randFaust(): second arg must be an"
                             " integer equal to RandFaustType.REAL,"
                             " RandFaustType.COMPLEX")
        rF = Faust(core_obj=FaustCorePy.FaustCore.randFaust(faust_type, field, min_num_factors, max_num_factors,
                                                            min_dim_size, max_dim_size, density))
        return rF

    def get_nb_rows(F):
        """
        Gives the Faust's number of rows.

        Returns:
                the number of rows.
        """
        return F.m_faust.shape()[0]

    def get_nb_cols(F):
        """
        Gives the Faust's number of columns.

        Returns:
                the number of columns.
        """
        return F.m_faust.shape()[1]

    def size(F):
        """
        Gives the size of the Faust.

        Returns:
            the Faust size tuple: get_nb_rows(), get_nb_cols().
        """
        #return F.shape[0], F.shape[1]
        return F.m_faust.shape()

    def transpose(F):
        """
        Transposes the current Faust.

        Examples:
            >>> F.transpose()
        """
        F_trans = Faust(core_obj=F.m_faust.transpose())
        return F_trans

    def conj(F):
        """
        Returns the conjugate of the current Faust.

        Args:
            F: the Faust object.

        Returns:
            the conjugate as Faust object.
            <br/> if F is a complex Faust, the value Fc returned is as
            Fc.todense().imag == - F.todense().imag and
            Fc.todense().real == F.todense().real.


        Examples:
            >>> from FaustPy import Faust, RandFaustType
            >>> F = Faust.randFaust(RandFaustType.MIXTE, RandFaustType.COMPLEX, 2,
            >>>                     5, 50, 100, .5)
            >>> Fc = F.conj()

        <b/> See also Faust.transpose
        """
        F_conj = Faust(core_obj=F.m_faust.conj())
        return F_conj

    def getH(F):
        """
        Returns the conjugate transpose of the current Faust.
        """
        F_ctrans = Faust(core_obj=F.m_faust.getH())
        return F_ctrans

    def display(F):
        """
        Displays information about F.


        WARNING: Currently a bug is affecting this function. When the Faust is
        transposed the dimensions are inverted in the display (like the Faust
        hasn't been transposed).


        Args:
            F: the Faust object.

        Examples:
            >>> from FaustPy import Faust, RandFaustType
            >>> F = Faust.randFaust(RandFaustType.MIXTE, RandFaustType.REAL,
            >>>                     1, 2, 50, 100, .5)
            >>> F.display()

        <b/>See also Faust.nnz, Faust.RCG, Faust.size, Faust.get_factor,
        Faust.get_nb_factors

        """
        F.m_faust.display()

    def __mul__(F, A):
        """
        Multiplies F by the numpy matrix A.

        This method overloads the Python operator *.

        Args:
            F: the Faust object.
            A: is a 2D numpy matrix (ndarray).
            <br/> A must be Fortran contiguous (i.e. Column-major order;
                `order' argument passed to np.ndararray() must be equal to str
                'F').

        Returns:
            the result of the multiplication as a numpy matrix.

        Raises:
            ValueError


        Examples:
            >>> from FaustPy import Faust, RandFaustType
            >>> import numpy as np

            >>> F = Faust.randFaust(RandFaustType.MIXTE, RandFaustType.REAL, 2,
            >>>                     5, 50, 100, .5)
            >>> A = np.random.rand(F.size()[1], 50)
            >>> B = F*A
            >>> # is equivalent to B = F.__mul__(A)

        """
        return F.m_faust.multiply(A)

    def todense(F):
        """
        Converts the current Faust into a numpy matrix.

        WARNING: this function costs F.get_nb_factors()-1 matrix multiplications.

        Returns:
            A numpy matrix.
        """
        identity = np.eye(F.get_nb_cols(), F.get_nb_cols())
        F_dense = F*identity
        return F_dense

    def __getitem__(F, indices):
        """
        Gets a submatrix of the full matrix of F.

        This function is a Python built-in overload.

        WARNING: this function costs as much as Faust.__mul__.

        Args:
            F: the Faust object.
            indices: array of length 2 which elements must be slice, integer or
            Ellipsis (...) (see examples below).

        Returns:
            the numpy subarray requested.

        Examples:
            >>> from FaustPy import Faust, RandFaustType
            >>> import numpy as np
            >>> from random import randint

            >>> F = Faust.randFaust(RandFaustType.MIXTE, RandFaustType.REAL, 2,
            >>>                        5, 50, 100, .5)
            >>> i1 = randint(0, min(F.size())-1)
            >>> i2 = randint(0, min(F.size())-1)

            >>> F[i1,i2] # element at line i1, column i2

            >>> F[:, i2] # full column i2

            >>> F[2:4, 1:4] # submatrix from line 2 to 3, each line containing
                            # only elements from column 1 to 3

            >>> F[::, 4:-1] # submatrix from line 0 to end line, each line
                            # containing only elements from column 4 to
                            # column before the last one.

            >>> F[0:i1, ...] # equivalent to F[0:i1, ::]
        """
        # check if indices has a 2 index (row and column)
        if (len(indices) != 2):
            raise ValueError('indices must contains 2 elements, '
                             'the row index and the col index')

        # check if the index are slice or integer or Ellipsis
        for id in indices:
            if (not isinstance(id, slice) and
                not isinstance(id, int) and
                id != Ellipsis):
                raise ValueError('indices must contains slice (1:n)'
                                 ' or Ellipsis (...) or int (2)')

        keyCol = indices[1]
        keyRow = indices[0]
        identity = np.eye(F.get_nb_cols(), F.get_nb_cols())
        if(keyCol != Ellipsis):
            identity = identity[..., keyCol]
        submatrix = F*identity
        submatrix = submatrix[keyRow, :]
        return submatrix

    def nnz(F):
        """
        Gives the number of non-zero elements in the Faust.

        The function sums together the number of non-zeros elements of
        each factor and returns the result. Note that in fact the sum is
        computed at Faust creation time and kept in cache.

        Returns:
            the number of non-zero elements as an integer.

        Examples:
            >>> F.nnz()
        """
        return F.m_faust.nnz()

    def density(F):
        """ Calculates the density of F.

        The density of F is equal to the number of non-zeros in factors over
        the total number of elements in dense matrix of F (which is equal to
        F.size()[0]*F.size()[1]).

        NOTE: this definition of density allows the value to be greater than 1.

        Args:
            F: the Faust object.

        Returns:
            the density value (float).

        Examples:
        >>> from FaustPy import Faust, RandFaustType
        >>> F = Faust.randFaust(RandFaustType.MIXTE, RandFaustType.COMPLEX, 2,
        >>>                     5, 50, 100, .5)
        >>> F.density()

        <b/> See also Faust.nnz, Faust.RCG
        """
        return float(F.nnz())/(F.get_nb_cols()*F.get_nb_rows())

    def RCG(F):
        """
        Computes the Relative Complexity Gain (inverse of the density).

        Returns:
            the RCG value (float).
            If the density is zero it will be float("inf").
            If the density is negative it will be -1.

        Examples:
            >>> F.RCG()
        """
        d = F.density()
        if(d > 0):
            return 1/d
        elif(d == 0):
            return float("inf")
        else:
            return float(-1)

    def norm(F, ord=2):
        """
        Computes the norm of the Faust. Several types of norm are available: 1-norm, 2-norm and Frobenius norm.

        Args:
            ord: (Optional) the norm order (1 or 2) or "fro" for
            Frobenius norm (by default the 2-norm is computed).

        Returns:
            the norm (float).

        Raises:
            ValueError.


        Examples:
            >>> F.norm()
            >>> F.norm(2)
            >>> F.norm('fro')
        """
        if(ord not in [1, 2, "fro"]):
            raise ValueError("ord must have the value 1, 2 or 'fro'.")
        return F.m_faust.norm(ord)

    def get_nb_factors(F):
        """
        Gives the Faust's number of factors.

        Returns:
            the number of factors.

        Examples:
            >>> F.get_nb_factors()
        """
        return F.m_faust.get_nb_factors()

    def get_factor(F, i):
        """
        Returns the Faust's i-th factor as a numpy.ndarray.

        Args:
            i: the integer index of the factor to get back.

        Returns:
            the i-th factor as a dense matrix (of type numpy.ndarray).

        Raises:
            ValueError.


        Examples:
            >>> F.get_factor(0)
        """
        fact = F.m_faust.get_fact(i)
        return fact

    def save(F, filepath, format="Matlab"):
        """
        Saves the Faust into file.

        Args:
            filepath: The path for saving the Faust (should ends with .mat
            if Matlab format used). It can be either an absolute or relative path.
            format: (Optional) It designates the format to use for
            writing. By default, it's "Matlab" to save the Faust in a .mat
            file (currently only that format is available).

        Raises:
            ValueError.


        Examples:
            >>> F.save("myFaust.mat")
        """
        if(format not in ["Matlab"]):
            raise ValueError("Only Matlab or Matlab_core format is supported.")
        if(format == "Matlab"):
            F.m_faust.save_mat_file(filepath)

    def isReal(F):
        """
        Returns True if F is a real Faust and False if it's a complex Faust.
        """
        return F.m_faust.isReal()

class RandFaustType:
    DENSE=0
    SPARSE=1
    MIXTE=2
    REAL=3
    COMPLEX=4
