# -*- coding: utf-8 -*-
##############################################################################
##                              Description:                                ##
##                                                                          ##
##          pyfaust is a python module  which delivers                      ##
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
    """FAµST Python wrapper main class.

    This class provides a numpy-like interface for operations
    with FAµST data structures, which correspond to matrices that can be
    written exactly as the product of sparse matrices.

    A FAµST data structure is designed to allow fast matrix-vector multiplications
    together with reduced memory storage compared to what would be obtained by
    manipulating directly the corresponding (dense) numpy array/matrix.

    A particular example is the matrix associated to the discrete Fourier
    transform, which can be represented exactly as a MuST, leading to a fast and
    compact implementation.

    Although the sparse matrices are more interesting for optimization it's not
    forbidden to define a Faust as a product of dense matrices or a mix of dense
    and sparse matrices.

    The matrices composing the Faust product, also called the factors, are
    defined on complex or real fields. Hence a Faust can be a complex Faust or a
    real Faust.

    A Faust has ideally to be seen and used as a numpy dense matrix/array, but
    this matrix exists only virtually and is actually represented by its factors.
    In order to use a Faust like a numpy matrix, a certain number of Python
    functions available for numpy matrices are implemented in this class but take
    note that not all are.

    It's possible to retrieve the dense matrix with the method Faust.toarray but
    it will cost the multiplication of the Faust's factors.
    It's noteworthy that in this documentation the expression 'dense matrix'
    designates the numpy dense matrix corresponding to a Faust, that is the
    matrix obtained by the multiplication of the previously mentioned Faust's
    factors.

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

    For more information about FAµST take a look at http://faust.inria.fr.
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
            Faust size 10x10, density 0.99, nnz_sum 99, 2 factors:<br/>
            FACTOR 0 (complex) SPARSE, size 10x10, density 0.4, nnz 40<br/>
            FACTOR 1 (complex) DENSE, size 10x10, density 0.59, nnz 59<br/>
            >>> G
            Faust size 19x16, density 1.37171, nnz_sum 417, 4 factors:<br/>
            FACTOR 0 (real) DENSE, size 19x17, density 0.49226, nnz 159<br/>
            FACTOR 1 (real) DENSE, size 17x10, density 0.517647, nnz 88<br/>
            FACTOR 2 (real) DENSE, size 10x13, density 0.515385, nnz 67<br/>
            FACTOR 3 (real) DENSE, size 13x16, density 0.495192, nnz 103<br/>

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
        elif(isinstance(dim_sizes, int)):
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

        <b/> See also Faust.nnz_sum, Faust.size
        """
        return F.m_faust.shape()

    @property
    def size(F):
        """
        Gives the number of elements in the Faust F.

        It's equivalent to numpy.prod(F.shape)).

        This function is intended to be used as a property (see the examples).

        Args:
            F: the Faust object.

        Returns:
            The number of elements in the Faust F.

        Examples:
            >>> from pyfaust import Faust
            >>> F = Faust.rand(2, 50, is_real=False)
            >>> size = F.size

        <b/> See also Faust.shape

        """
        return np.prod(F.shape)

    def transpose(F):
        """
        Returns the transpose of the Faust F.todense().

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

        <b/> See also Faust.transpose, Faust.get_factor, Faust.get_num_factors, Faust.getH, Faust.H
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
            Faust size 98x82, density 0.686909, nnz_sum 5520, 2 factors:<br/>
            FACTOR 0 (real) SPARSE, size 98x78, density 0.395081, nnz 3020<br/>
            FACTOR 1 (real) SPARSE, size 78x82, density 0.390869, nnz 2500<br/>

            >>> F
            Faust size 98x82, density 0.686909, nnz_sum 5520, 2 factors:<br/>
            FACTOR 0 (real) SPARSE, size 98x78, density 0.395081, nnz 3020<br/>
            FACTOR 1 (real) SPARSE, size 78x82, density 0.390869, nnz 2500<br/>
            <!-- >>> -->

        <b/>See also Faust.nnz_sum, Faust.density, Faust.shape, Faust.get_factor,
        <b/>Faust.get_num_factors, Faust.__repr__

        """
        print(F.__repr__())
        #F.m_faust.display()

    def __mul__(F, A):
        """
        Multiplies F by A which is a dense matrix or a Faust object.

        This method overloads a Python function/operator.

        NOTE: The primary goal of this function is to implement “fast” multiplication by a
        MuST, which is the operation performed when A is a standard matrix (dense or
        sparse).
        In the best cases, F*A is F.rcg() times faster than performing the
        equivalent F.toarray()*A.
        <br/>When A is a MuST, F*A is itself a MuST.

        Args:
            F: the Faust object.
            A: is a Faust object or a 2D dense matrix (numpy.ndarray, numpy.matrix).
            <br/> In the latter case, A must be Fortran contiguous (i.e. Column-major order;
                `order' argument passed to np.ndararray() must be equal to str
                'F').

        Returns:
            the result of the multiplication as a numpy.ndarray if A is a
            ndarray.<br/>
            The result of the multiplication as a Faust object if A is a Faust.

        Raises:
            ValueError


        Examples:
            >>> from pyfaust import Faust
            >>> import numpy as np

            >>> F = Faust.rand(5, [50, 100])
            >>> A = np.random.rand(F.shape[1], 50)
            >>> B = F*A
            >>> # is equivalent to B = F.__mul__(A)
            >>> G = Faust.rand(5, F.shape[1])
            >>> H = F*G
            >>> # H is a Faust because F and G are

        <b/>See also Faust.__init__, Faust.rcg
        """
        if(isinstance(A, Faust)):
            if(F.shape[1] != A.shape[0]): raise ValueError("The dimensions of "
                                                          "the two Fausts must "
                                                          "agree.")
            return Faust(core_obj=F.m_faust.multiply(A.m_faust))
        else:
            return F.m_faust.multiply(A)

    def toarray(F):
        """
        Converts the current Faust into a numpy array.

        WARNING: Using this function is discouraged except for test purposes,
        as it loses the main potential interests of the MuST structure:
        compressed memory storage and faster matrix-vector multiplication
        compared to its equivalent dense matrix representation.

        Returns:
            A numpy ndarray.

        Raises:
            MemoryError


        Examples:
            >>> from pyfaust import Faust
            >>> F = Faust.rand(5, 10**6, .00001, 'sparse')
               Faust size 1000000x1000000, density 5e-05, nnz_sum 49999995, 5 factors:<br/>
                FACTOR 0 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999<br/>
                FACTOR 1 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999<br/>
                FACTOR 2 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999<br/>
                FACTOR 3 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999<br/>
                FACTOR 4 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999<br/>
            >>> # an attempt to convert F to a dense matrix is most likely to raise a memory error
            >>> # the sparse format is the only way to handle a so big Faust
            >>> F.toarray()
            ...
            MemoryError

        <b/>See also Faust.todense
        """
        identity = np.eye(F.shape[1], F.shape[1])
        F_dense = F*identity
        return F_dense

    def todense(F):
        """
        Converts the current Faust into a numpy matrix.

        WARNING: Using this function is discouraged except for test purposes,
        as it loses the main potential interests of the MuST structure:
        compressed memory storage and faster matrix-vector multiplication
        compared to its equivalent dense matrix representation.

        Returns:
            A numpy matrix M which is such that M*x == F*x for any vector x.

        Raises:
            MemoryError


        Examples:
            >>> from pyfaust import Faust
            >>> F = Faust.rand(5, 10**6, .00001, 'sparse')
            >>> F
            Faust size 1000000x1000000, density 5e-05, nnz_sum 49999995, 5 factors:<br/>
            FACTOR 0 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999<br/>
            FACTOR 1 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999<br/>
            FACTOR 2 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999<br/>
            FACTOR 3 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999<br/>
            FACTOR 4 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999<br/>
            >>> # an attempt to convert F to a dense matrix is most likely to raise a memory error
            >>> # the sparse format is the only way to handle a so big Faust
            >>> F.todense()
            (...)
            MemoryError

        <b/>See also Faust.toarray
        """
        return np.matrix(F.toarray())


    def __getitem__(F, indices):
        """
        Returns a Faust representing a submatrix of the dense matrix of F.

        This function is a Python built-in overload.

        WARNING: this function doesn't handle a slice step different to 1 (e.g. F[i:j:2,:]
        where slice step is 2.)

        Args:
            F: the Faust object.
            indices: array of length 1 or 2 which elements must be slice, integer or
            Ellipsis (...) (see examples below). Note that using Ellipsis for
            more than two indices is forbidden.

        Returns:
            the Faust object requested.

        Raises:
            IndexError


        Examples:
            >>> from pyfaust import Faust
            >>> import numpy as np
            >>> from random import randint

            >>> F = Faust.rand(5, [50, 100])
            >>> i1 = randint(0, min(F.shape)-1)
            >>> i2 = randint(0, min(F.shape)-1)

            >>> F[i1,i2] # a Faust representing a matrix with only one element
                         # at row i1, column i2 of F's dense matrix

            >>> F[:, i2] # full column i2

            >>> F[2:4, 1:4] # from row 2 to 3, each row containing
                            # only elements from column 1 to 3

            >>> F[::, 4:-1] # from row 0 to end row, each row
                            # containing only elements from column 4 to
                            # column before the last one.

            >>> F[0:i1, ...] # equivalent to F[0:i1, ::]
        """
        #TODO: refactor (by index when indices == tuple(2), error message,
        #      out_indices checking on the end)
        #TODO: check that step == 1 on slices and stop on failure if it is or
        # implement negative steps
        if(indices == Ellipsis): # F[...]
            out_indices = [slice(0,F.shape[0]), slice(0, F.shape[1])]
        elif(isinstance(indices,int)): # F[i] # a line
            out_indices = [slice(indices, indices+1), slice(0, F.shape[1])]
        elif(isinstance(indices,slice)):
            #F[i:j] a group of contiguous lines
            out_indices = [indices, slice(0, F.shape[1])]
        elif(isinstance(indices, tuple)):
            if(len(indices) == 2):
                out_indices = [0,0]
                if(indices[0] == Ellipsis):
                    if(indices[1] == Ellipsis):
                        raise IndexError('an index can only have a single ellipsis '
                                         '(\'...\')')
                    else:
                        # all lines
                        out_indices[0] = slice(0, F.shape[0])
                elif(isinstance(indices[0], int)):
                    # line F[i]
                    out_indices[0] = slice(indices[0], indices[0]+1)
                elif(isinstance(indices[0], slice)):
                    out_indices[0] = indices[0]
                else:
                     raise IndexError("only integers, slices (`:`), ellipsis"
                                     " (`...`), and integer are valid indices")
                if(indices[1] == Ellipsis):
                    # all lines
                    out_indices[1] = slice(0, F.shape[1])
                elif(isinstance(indices[1], int)):
                    # line F[i]
                    out_indices[1] = slice(indices[1], indices[1]+1)
                elif(isinstance(indices[1], slice)):
                    out_indices[1] = indices[1]
                else:
                     raise IndexError("only integers, slices (`:`), ellipsis"
                                    " (`...`), and integer are valid indices")
            else:
                raise IndexError('Too many indices.')
        else:
            print("type indices:", type(indices))
            raise IndexError("only integers, slices (`:`), ellipsis"
                             " (`...`), and integer are valid indices")
        if(out_indices[0].start == None and out_indices[0].stop == None): #F[::] or F[::,any]
            out_indices[0] = slice(0,F.shape[0])
        if(out_indices[0].stop < 0):
            out_indices[0] = slice(out_indices[0].start,
                                   F.shape[0]+out_indices[0].stop)
        if(out_indices[1].start == None and out_indices[1].stop == None): #F[any, ::]
            out_indices[1] = slice(0,F.shape[1])
        if(out_indices[1].stop < 0):
            out_indices[1] = slice (out_indices[1].start,
                                    F.shape[1]+out_indices[1].stop)
        for i in range(0,2):
            if(out_indices[i].start >= F.shape[i] or out_indices[i].stop > F.shape[i]):
                raise IndexError("index "+
                                 str(max(out_indices[i].start,out_indices[i].stop-1))+
                                 " is out of bounds for axis "+str(i)+" with size "+
                                 str(F.shape[i]))
        slice_F = Faust(core_obj=F.m_faust.slice(out_indices[0].start,
                                         out_indices[0].stop,
                                         out_indices[1].start,
                                         out_indices[1].stop))
        return slice_F

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
        """ Calculates the density of F such that F.nnz_sum() == F.density()/F.size.

        NOTE: A value of density below one indicates potential memory savings
        compared to storing the corresponding dense matrix F.todense(), as well
        as potentially faster matrix-vector multiplication when applying F*x
        instead of F.todense()*x.

        NOTE: A density above one is possible but prevents any saving.

        Args:
            F: the Faust object.

        Returns:
            the density value (float).

        Examples:
        >>> from pyfaust import Faust
        >>> F = Faust.rand(5, 50, .5)
        >>> dens = F.density()

        <b/> See also Faust.nnz_sum, Faust.rcg, Faust.size, Faust.todense
        """
        return float(F.nnz_sum())/F.size

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
            F.rcg() == 1/F.density().

        <b/> See also: Faust.density, Faust.nnz_sum, Faust.shape.
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
        Computes the norm of F.todense(). Several types of norm are available: 1-norm, 2-norm and Frobenius norm.

        NOTE: The norm of F is equal to the norm of its dense matrix
        F.toarray().

        WARNING: This function costs at least as much as Faust.__mul__.

        Args:
            F: the Faust object.
            ord: (optional) the norm order (1 or 2) or "fro" for
            Frobenius norm (by default the 2-norm is computed).

        Returns:
            the norm (float).
            <br/>If ord == 1, then the norm is the maximum absolute column sum of the
            matrix F.todense().
            <br/>If ord == 2, then the norm is approximately
            max(scipy.linalg.svd(F.toarray())[1]) (the greatest singular value). This is
            equivalent to norm(F).
            <br/>If ord == 'fro', then the norm is the Frobenius norm of
            F.todense().

        Raises:
            ValueError.


        Examples:
            >>> from pyfaust import Faust
            >>> F = Faust.rand([1, 2], [50, 100], .5)
            >>> F.norm()
            80.0914231648143
            >>> F.norm(2)
            80.0914231648143
            >>> F.norm('fro')
            133.96231832406144
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

        NOTE: storing F should typically use rcg(F) times less storage than
        storing F.todense().

        Args:
            F: the Faust object.
            filepath: the path for saving the Faust (should end with .mat
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

            <b/> See also Faust.__init__, Faust.rcg.
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

        NOTE: if F is a real Faust, all factors are real. If F is a complex
        Faust, at least one factor is complex.

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
