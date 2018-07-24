# -*- coding: utf-8 -*-
##############################################################################
##                              Description:                                ##
##                                                                          ##
##          pyfaust is a python module  which delivered                     ##
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
    """FAµST Python class.
     This class represents a given dense matrix by a product of sparse matrices
     (i.e. Faust).
     The main goal of Faust representation is to speed up operations on that
     matrix, especially the multiplication. Besides the time optimization, a Faust
     can reduce the memory space size needed both for storage and loading.

     Although the sparse matrices are more interesting for optimization it's not
     forbidden to define a Faust as a product of dense matrices or a mix of dense
     and sparse matrices.

     The matrices composing the Faust product, also called the factors, are
     defined on complex or real fields. Hence a Faust can be a complex Faust or a
     real Faust.

     A Faust has ideally to be seen and used as a numpy dense matrix/array, but
     this matrix exists only virtually and is actually represented by its factors.
     In order to use a Faust like a numpy matrix, a certain number of Python
     built-ins available for numpy matrices are implemented in this class but take
     note that not all are.

     You have the capability to retrieve the dense matrix with the
     method Faust.toarray but it will cost the multiplication of the Faust's factors.
     It's noteworthy that in this documentation the expression 'dense matrix'
     designates the numpy dense matrix corresponding to a Faust, that is the
     matrix obtained by the multiplication of the
     previously mentioned Faust's factors.

     Likewise, other Faust's methods need to calculate the factor product, and
     they will be indicated with a warning in this documentation. You should avoid
     to use them if it's not really necessary (for example you might limit their use
     to test purposes).

     Other important limitation is that contrary to a numpy dense matrix
     a Faust is immutable. It means that you cannot affect elements of a Faust using
     the affectation operator `=' like you do with a numpy matrix (e.g. `M[i,j] =
     2').
     That limitation is the reason why the Python built-in `__setitem__()' is not
     implemented in this class.

    """

    def  __init__(F, factors=None, scale=1.0, filepath=None, core_obj=None):
        """ Creates a Faust from a list of factors or alternatively from a file.

            Another easy way to create a Faust is to call the static method Faust.rand().

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
            scale: a multiplicative scalar (see examples below).
            core_obj: for internal purpose only. Please don't fill this argument.

        WARNING: filepath and factors arguments are multually exclusive. If you
        use filepath you must explicitely set argument with the keyword.

        Examples:
            >>> from pyfaust import Faust
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

            >>> scale = 2
            >>> G = Faust(factors, scale) # G == scale*F

            >>> F.save("F.mat")
            >>> # define a Faust from file
            >>> H = Faust(filepath="F.mat")
            >>> I = Faust(filepath="F.mat", scale) # I == scale*H

        <b/> See also Faust.save, Faust.rand

        """
        if(core_obj):
            F.m_faust = core_obj
        else:
            if(filepath and isinstance(filepath, str)):
                    contents = loadmat(filepath)
                    factors = contents['faust_factors'][0]
            if(factors is not None):
                F.m_faust = FaustCorePy.FaustCore(factors, scale);
            #else:
                #TODO: manage empty Faust

    @staticmethod
    def rand(num_factors, dim_sizes, density=0.1, fac_type="mixed",
                  is_real=True):
        """
        Generates a random Faust.

            Args:
                num_factors: If it's an integer it will be the number of random factors to set in the Faust.
                            If num_factors is a list or tuple of 2 integers then the
                            number of factors will be set randomly between
                            num_factors[0] and num_factors[1] (inclusively).
                dim_sizes: if it's an integer it will be the order of the square
                matrix factors (of size size_dims**2).
                            If it's a list or tuple of 2 integers then the
                            number of rows and columns will
                            be a random number between size_dims[0] and
                            size_dims[1] (inclusively).
                density: (optional) the approximate density of factors. The
                default value is 0.1. It should be a floating point number
                between 0 and 1.
                fac_type: (optional) the type of density of factors. Must be
                "sparse", "dense" or "mixed" if you want a mix of dense and
                sparse matrices in the generated Faust (choice's done according
                to an uniform distribution).
                            The default value is "mixed".
                is_real: (optional) a boolean set to True to generate a real Faust,
                        set to False to generate a complex Faust. The default value
                        is True.

        Returns:
            the random Faust.

        Examples:
            >>> from pyfaust import Faust
            >>> F = Faust.rand(2, 10, .5, is_real=False)
            >>> G = Faust.rand([2, 5], [10, 20], .5, fac_type="dense")
            >>> F
            Faust of size : 10x10, nb factor 2, RCG 1.05263,nnz 95<br/>
             FACTOR 0 type : SPARSE size 10x10, density 0.52, nnz 52<br/>
             FACTOR 1 type : DENSE size 10x10, density 0.43, nnz 43<br/>
            >>> G
            Faust of size : 18x20, nb factor 5, RCG 0.590164,nnz 610<br/>
             FACTOR 0 type : DENSE size 18x17, density 0.496732, nnz 152<br/>
             FACTOR 1 type : DENSE size 17x17, density 0.49481, nnz 143<br/>
             FACTOR 2 type : DENSE size 17x10, density 0.470588, nnz 80<br/>
             FACTOR 3 type : DENSE size 10x16, density 0.48125, nnz 77<br/>
             FACTOR 4 type : DENSE size 16x20, density 0.49375, nnz 158<br/>

        <b/> See also Faust.__init__
        """
        DENSE=0
        SPARSE=1
        MIXED=2
        REAL=3
        COMPLEX=4
        # set density type of factors
        if(not isinstance(fac_type, str) or fac_type not in ["sparse",
                                                             "dense",
                                                             "mixed"]):
            raise ValueError('Faust.rand(): argument fac_type must be a'
                             ' str equal to "sparse",'
                             ' "dense" or "mixed".')

        fac_type_map = {"sparse" : SPARSE, "dense" : DENSE, "mixed" : MIXED}
        # set field of matrices/factors
        if(not isinstance(is_real, bool)):
            raise ValueError('Faust.rand(): argument is_real must be a'
                             'boolean.')
        if(is_real):
            field = REAL
        else:
            field = COMPLEX
        if((isinstance(num_factors, list) or isinstance(num_factors, tuple)) and
        len(num_factors) == 2):
            min_num_factors = num_factors[0]
            max_num_factors = num_factors[1]
        elif(isinstance(num_factors, int)):
            min_num_factors = max_num_factors = num_factors
        else:
            raise ValueError("Faust.rand(): num_factors argument must be an "
                             "integer or a list/tuple of 2 integers.")
        if((isinstance(dim_sizes, list) or isinstance(dim_sizes, tuple)) and
        len(dim_sizes) == 2):
            min_dim_size = dim_sizes[0]
            max_dim_size = dim_sizes[1]
        elif(isinstance(num_factors, int)):
            min_dim_size = max_dim_size = dim_sizes
        else:
            raise ValueError("Faust.rand(): dim_sizes argument must be an "
                             "integer or a list/tuple of 2 integers.")
        rF = Faust(core_obj=FaustCorePy.FaustCore.randFaust(
            fac_type_map[fac_type], field, min_num_factors, max_num_factors,
            min_dim_size, max_dim_size, density))
        return rF

    @property
    def shape(F):
        """
        Gives the size of the Faust F.

        This function is intended to be used as a property (see the examples).

        The size is a pair of numbers: the number of rows and the number of
        columns of the equivalent dense matrix of F.

        Args:
            F: the Faust object.

        Returns:
            the Faust shape tuple, with at first index the number of rows, and
            at second index the number of columns.

        Examples:
            >>> from pyfaust import Faust
            >>> F = Faust.rand(2, 50, is_real=False)
            >>> nrows, ncols = F.shape
            >>> nrows = F.shape[0]
            >>> ncols = F.shape[1]

        <b/> See also Faust.display
        """
        return F.m_faust.shape()

    def transpose(F):
        """
        Returns the transpose of the Faust F.

        Args:
            F: the Faust object.

        Returns:
            F transpose as a Faust object.

        Examples:
            >>> tF = F.transpose()

        <b/> Faust.getH, Faust.H, Faust.T
        """
        F_trans = Faust(core_obj=F.m_faust.transpose())
        return F_trans

    @property
    def T(F):
        """
        Returns the transpose of the Faust F.

        This function is intended to be used as a property (see the examples).

        Args:
            F: the Faust object.

        Returns:
            F transpose as a Faust object.

        Examples:
            >>> tF = F.T

        <b/> Faust.getH, Faust.H, Faust.T
        """
        return F.transpose()

    def conj(F):
        """
        Returns the conjugate of the current Faust.

        Args:
            F: the Faust object.

        Returns:
            the conjugate as a Faust object.
            <br/> if F is a real Faust then F_conj == F.
            <br/> if F is a complex Faust, the value Fc returned verifies the
            next assertion for i in range(0,F.get_num_factors()):

            <code>F.get_factor(i).conj() == Fc.get_factor(i)</code>


        Examples:
            >>> from pyfaust import Faust
            >>> F = Faust.rand(5, 50, is_real=False)
            >>> Fc = F.conj()

        <b/> See also Faust.transpose, Faust.get_factor, Faust.get_num_factors,
        Faust.getH, Faust.H
        """
        F_conj = Faust(core_obj=F.m_faust.conj())
        return F_conj

    def getH(F):
        """
        Returns the conjugate transpose of F.

        Args:
            F: the Faust object.

        Returns:
            the conjugate transpose of F as a Faust object.


        Examples:
            >>> from pyfaust import Faust
            >>> F = Faust.rand(5, 50, is_real=False)
            >>> H1 = F.getH()
            >>> H2 = F.transpose()
            >>> H2 = H2.conj()
            >>> (H1.toarray() == H2.toarray()).all()
            True

        <b/> See also Faust.transpose, Faust.conj, Faust.H
        """
        F_ctrans = Faust(core_obj=F.m_faust.getH())
        return F_ctrans

    @property
    def H(F):
        """
        Returns the conjugate transpose of F.

        This function is intended to be used as a property (see the examples).

        Returns:
            the conjugate transpose of F as a Faust object.


        Examples:
            >>> from pyfaust import Faust
            >>> F = Faust.rand(5, 50, is_real=False)
            >>> H1 = F.H
            >>> H2 = F.transpose()
            >>> H2 = H2.conj()
            >>> (H1.toarray() == H2.toarray()).all()
            True

        <b/> See also Faust.transpose, Faust.conj, Faust.getH
        """
        return F.getH()

    def __repr__(F):
        """
        Returns a str object representing the Faust object.

        This method overloads a Python function.

        NOTE: Ideally this function is intended to return a valid Python
        expression but here this is not the case. Only information is
        displayed.

        Args:
            F: the Faust object.

        Examples:
            >>> from pyfaust import Faust
            >>> F = Faust.rand(2, 50)
            >>> F.__repr__()
            >>> # the same function is called when typing F in a terminal:
            >>> F

        <b/>See also Faust.nnz_sum, Faust.rcg, Faust.shape, Faust.get_factor,
        <b/>Faust.get_num_factors, Faust.display

        """
        #_str = super(object, F).__repr__()
        _str = str(F.m_faust.to_string())
        return _str

    def display(F):
        """
        Displays information about F.

        Args:
            F: the Faust object.

        Examples:
            >>> from pyfaust import Faust
            >>> F = Faust.rand([1, 2], [50, 100], .5)
            >>> F.display()
            >>> # equivalent to:
            >>> F

        <b/>See also Faust.nnz_sum, Faust.rcg, Faust.shape, Faust.get_factor,
        <b/>Faust.get_num_factors, Faust.__repr__

        """
        print(F.__repr__())
        #F.m_faust.display()

    def __mul__(F, A):
        """
        Multiplies F by the numpy matrix A.

        This method overloads a Python function/operator.

        WARNING: this function costs F.get_num_factors() matrix
        multiplications. its use is discouraged except for test purpose.

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
            >>> from pyfaust import Faust
            >>> import numpy as np

            >>> F = Faust.rand(5, [50, 100])
            >>> A = np.random.rand(F.shape[1], 50)
            >>> B = F*A
            >>> # is equivalent to B = F.__mul__(A)

        """
        return F.m_faust.multiply(A)

    def toarray(F):
        """
        Converts the current Faust into a numpy array.

        WARNING: this function costs F.get_num_factors()-1 matrix multiplications.

        Returns:
            A numpy ndarray.
        """
        identity = np.eye(F.shape[1], F.shape[1])
        F_dense = F*identity
        return F_dense

    def todense(F):
        """
        Converts the current Faust into a numpy matrix.

        WARNING: this function costs F.get_num_factors()-1 matrix multiplications.
        Besides it implies the loss of space compression allowed by the
        Faust representation.
        Using this function is discouraged except for test purpose.

        Returns:
            A numpy matrix M which is such that M*x == F*x for any vector x.
        """
        return np.matrix(F.toarray())


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
            >>> from pyfaust import Faust
            >>> import numpy as np
            >>> from random import randint

            >>> F = Faust.rand(5, [50, 100])
            >>> i1 = randint(0, min(F.shape)-1)
            >>> i2 = randint(0, min(F.shape)-1)

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
        identity = np.eye(F.shape[1], F.shape[1])
        if(keyCol != Ellipsis):
            identity = identity[..., keyCol]
        submatrix = F*identity
        submatrix = submatrix[keyRow, :]
        return submatrix

    def nnz_sum(F):
        """
        Gives the total number of non-zero elements in the factors of F.

        The function sums together the number of non-zeros elements of
        each factor and returns the result. Note that in fact the sum is
        computed at Faust creation time and kept in cache.

        Returns:
            the number of non-zeros.

        <b/> See also Faust.rcg, Faust.density.
        """
        return F.m_faust.nnz()

    def density(F):
        """ Calculates the density of F.

        The density of F is equal to the number of non-zeros in factors over
        the total number of elements in dense matrix of F (which is equal to
        F.shape[0]*F.shape[1]).

        NOTE: this definition of density allows the value to be greater than 1.

        Args:
            F: the Faust object.

        Returns:
            the density value (float).

        Examples:
        >>> from pyfaust import Faust
        >>> F = Faust.rand(5, 50, .5)
        >>> dens = F.density()

        <b/> See also Faust.nnz_sum, Faust.rcg
        """
        return float(F.nnz_sum())/(F.shape[1]*F.shape[0])

    def rcg(F):
        """
        Computes the Relative Complexity Gain (inverse of Faust.density).

        RCG is the theoretical gain brought by Faust representation relatively to its dense
        matrix equivalent. The higher is the RCG, the more computational
        savings will be made. That gain applies both for storage and computation time.

        Args:
            F: the Faust object.

        Returns:
            the RCG value (float).
            If the density is zero it will be float("inf").
            If the density is negative it will be -1.

        Examples:
            >>> F.rcg()

        <b/> See also: Faust.density, Faust.nnz_sum.
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
        Computes the norm of F. Several types of norm are available: 1-norm, 2-norm and Frobenius norm.

        NOTE: The norm of F is equal to the norm of its dense matrix.

        WARNING: this function costs at least as much as Faust.__mul__.

        Args:
            F: the Faust object.
            ord: (optional) the norm order (1 or 2) or "fro" for
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

    def get_num_factors(F):
        """
        Gives the number of factors of F.

        Returns:
            the number of factors.

        Examples:
        >>> from pyfaust import Faust
        >>> F = Faust.rand(2, 100, .5)
        >>> nf = F.get_num_factors()

        <b/> See also Faust.get_factor
        """
        return F.m_faust.get_nb_factors()

    def get_factor(F, i):
        """
        Returns the i-th factor of F.

        Args:
            F: the Faust object.
            i: the factor index.

        Returns:
            the i-th factor as a dense matrix (of type numpy.ndarray).

        Raises:
            ValueError.


        Examples:
            >>> from pyfaust import Faust
            >>> F = Faust.rand(5, [50, 100], .5)
            >>> f0 = F.get_factor(0)

        <b/> See also Faust.get_num_factors
        """
        fact = F.m_faust.get_fact(i)
        return fact

    def save(F, filepath, format="Matlab"):
        """
        Saves the Faust F into file.

        Args:
            F: the Faust object.
            filepath: the path for saving the Faust (should ends with .mat
            if Matlab format is used).
            format: (optional) it designates the format to use for
            writing. By default, it's "Matlab" to save the Faust in a .mat
            file (currently only that format is available).

        Raises:
            ValueError.


        Examples:
            >>> from pyfaust import Faust
            >>> F = Faust.rand([2, 3], [10, 20],.5, is_real=False)
            >>> F.save("F.mat")
            >>> G = Faust(filepath="F.mat")
            >>> H = Faust(filepath="F.mat", scale=2)
            >>> (H.toarray()/G.toarray() != 2).any()
            False

            <b/> See also Faust.__init__.
        """
        if(format not in ["Matlab"]):
            raise ValueError("Only Matlab or Matlab_core format is supported.")
        if(format == "Matlab"):
            F.m_faust.save_mat_file(filepath)

    @property
    def dtype(F):
        """
        Returns the dtype of the Faust.

        This function is intended to be used as a property (see the examples).

        Args:
            F: the Faust object.

        Returns:
            A numpy.complex128 if the Faust is complex otherwise (if it's
            a real Faust) it returns a numpy.float64.


        Examples:
            >>> from pyfaust import Faust
            >>> F = Faust.rand([2, 3], [10, 20],.5, is_real=False)
            dtype('complex128')
            >>> F = Faust.rand([2, 3], [10, 20],.5)
            dtype('float64')


        """
        if(F.m_faust.isReal()):
            return np.dtype(np.float64)
        else:
            return np.dtype(np.complex)
