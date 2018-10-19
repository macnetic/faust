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
    transform, which can be represented exactly as a Faust, leading to a fast and
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

            Another easy way to create a Faust is to call the static method FaustFactory.rand().

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

        WARNING: filepath and factors arguments are mutually exclusive. If you
        use filepath you must explicitely set argument with the keyword.

        Examples:
            >>> from pyfaust import FaustFactory
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

        <b/> See also Faust.save, FaustFactory.rand

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
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand(2, 50, is_real=False)
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
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand(2, 50, is_real=False)
            >>> size = F.size

        <b/> See also Faust.shape

        """
        return np.prod(F.shape)

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
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand(5, 50, is_real=False)
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
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand(5, 50, is_real=False)
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
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand(5, 50, is_real=False)
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
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand(2, 50)
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
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand([1, 2], [50, 100], .5)
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
        Multiplies F by A which is a dense matrix, a Faust object or a scalar number.

        This method overloads a Python function/operator.

        NOTE: The primary goal of this function is to implement “fast” multiplication by a
        Faust, which is the operation performed when A is a standard matrix (dense or
        sparse).
        In the best cases, F*A is F.rcg() times faster than performing the
        equivalent F.toarray()*A.
        <br/>When A is a Faust, F*A is itself a Faust.

        Args:
            F: the Faust object.
            A: is a scalar number or a Faust object or a 2D dense matrix (numpy.ndarray,
            numpy.matrix).
            <br/> In the latter case, A must be Fortran contiguous (i.e. Column-major order;
                `order' argument passed to np.ndararray() must be equal to str
                'F').

        Returns:
            the result of the multiplication as a numpy.ndarray if A is a
            ndarray.<br/>
            The result of the multiplication as a Faust object if A is a Faust
            or a scalar number.

        Raises:
            ValueError


        Examples:
            >>> from pyfaust import FaustFactory
            >>> import numpy as np

            >>> F = FaustFactory.rand(5, [50, 100])
            >>> A = np.random.rand(F.shape[1], 50)
            >>> B = F*A
            >>> # is equivalent to B = F.__mul__(A)
            >>> G = FaustFactory.rand(5, F.shape[1])
            >>> H = F*G
            >>> # H is a Faust because F and G are

        <b/>See also Faust.__init__, Faust.rcg
        """
        if(isinstance(A, Faust)):
            if(F.shape[1] != A.shape[0]): raise ValueError("The dimensions of "
                                                          "the two Fausts must "
                                                          "agree.")
            return Faust(core_obj=F.m_faust.multiply(A.m_faust))
        elif(isinstance(A, float) or isinstance(A, int) or isinstance(A, np.complex)):
            return Faust(core_obj=F.m_faust.multiply_scal(A))
        else:
            return F.m_faust.multiply(A)

    def toarray(F):
        """
        Converts the current Faust into a numpy array.

        WARNING: Using this function is discouraged except for test purposes,
        as it loses the main potential interests of the Faust structure:
        compressed memory storage and faster matrix-vector multiplication
        compared to its equivalent dense matrix representation.

        Returns:
            A numpy ndarray.

        Raises:
            MemoryError


        Examples:
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand(5, 10**6, .00001, 'sparse')
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
        as it loses the main potential interests of the Faust structure:
        compressed memory storage and faster matrix-vector multiplication
        compared to its equivalent dense matrix representation.

        Returns:
            A numpy matrix M which is such that M*x == F*x for any vector x.

        Raises:
            MemoryError


        Examples:
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand(5, 10**6, .00001, 'sparse')
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
        Returns a Faust representing a submatrix of the dense matrix of F or a scalar element if that Faust can be reduced to a single element.

        This function is a Python built-in overload.

        WARNING:
            - This function doesn't handle a slice step different to 1 (e.g. F[i:j:2,:]
        where slice step is 2.).
            - It is not advised to use this function as an element accessor
        (e.g. F(0,0)) because such an use induces to convert the Faust to its
        dense matrix representation and that is a very expensive computation if used
        repetitively.

        Args:
            F: the Faust object.
            indices: array of length 1 or 2 which elements must be slice, integer or
            Ellipsis (...) (see examples below). Note that using Ellipsis for
            more than two indices is forbidden.

        Returns:
            the Faust object requested or just the corresponding scalar if that Faust has
            a shape equal to (1,1).

        Raises:
            IndexError


        Examples:
            >>> from pyfaust import FaustFactory
            >>> import numpy as np
            >>> from random import randint

            >>> F = FaustFactory.rand(5, [50, 100])
            >>> i1 = randint(0, min(F.shape)-1)
            >>> i2 = randint(0, min(F.shape)-1)

            >>> F[i1,i2] # is the scalar element located at
                         # at row i1, column i2 of the F's dense matrix

            >>> F[:, i2] # full column i2

            >>> F[2:4, 1:4] # from row 2 to 3, each row containing
                            # only elements from column 1 to 3

            >>> F[::, 4:-1] # from row 0 to end row, each row
                            # containing only elements from column 4 to
                            # column before the last one.

            >>> F[0:i1, ...] # equivalent to F[0:i1, ::]
            >>> F[2::, :3:] # equivalent to F[2:F.shape[0],0:3]
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
        elif(out_indices[0].start == None): # out_indices[0].stop != None
            out_indices[0] = slice(0, out_indices[0].stop)
        elif(out_indices[0].stop == None): # out_indices[0].start != None
            out_indices[0] = slice(out_indices[0].start, F.shape[0])
        if(out_indices[0].stop < 0):
            out_indices[0] = slice(out_indices[0].start,
                                   F.shape[0]+out_indices[0].stop)
        if(out_indices[1].start == None and out_indices[1].stop == None): #F[any, ::]
            out_indices[1] = slice(0,F.shape[1])
        elif(out_indices[1].start == None): # out_indices[1].stop != None
            out_indices[1] = slice(1, out_indices[1].stop)
        elif(out_indices[1].stop == None): # out_indices[1].start != None
            out_indices[1] = slice(out_indices[1].start, F.shape[1])
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
        >>> from pyfaust import FaustFactory
        >>> F = FaustFactory.rand(5, 50, .5)
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

    def norm(F, ord='fro'):
        """
        Computes the norm of F.todense(). Several types of norm are available: 1-norm, 2-norm and Frobenius norm.

        NOTE: The norm of F is equal to the norm of its dense matrix
        F.toarray().

        WARNING: This function costs at least as much as Faust.__mul__.

        Args:
            F: the Faust object.
            ord: (optional) the norm order (1 or 2) or "fro" for
            Frobenius norm (by default the Frobenius norm is computed).

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
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand([1, 2], [50, 100], .5)
            >>> F.norm()
            133.96231832406144
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
        >>> from pyfaust import FaustFactory
        >>> F = FaustFactory.rand(2, 100, .5)
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
            a copy of the i-th factor, the copy type is:
                - numpy.ndarray if it is a full storage matrix or,
                - scipy.sparse.csc.matrix_csc if it's a sparse matrix of a
                transposed Faust,
                - scipy.sparse.csr.csr_matrix if it's a sparse matrix of a
                non-transposed Faust.

        Raises:
            ValueError.


        Examples:
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand(5, [50, 100], .5)
            >>> f0 = F.get_factor(0)

        <b/> See also Faust.get_num_factors, Faust.transpose
        """
        fact = F.m_faust.get_fact_opt(i)
        return fact

    def get_factor_nonopt(F, i):
        """
        DEPRECATED: use Faust.get_factor
        Returns the i-th factor of F.

        Args:
            F: the Faust object.
            i: the factor index.

        Returns:
            a copy of the i-th factor as a dense matrix (of type numpy.ndarray).

        Raises:
            ValueError.


        Examples:
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand(5, [50, 100], .5)
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
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand([2, 3], [10, 20],.5, is_real=False)
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
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand([2, 3], [10, 20],.5, is_real=False)
            dtype('complex128')
            >>> F = FaustFactory.rand([2, 3], [10, 20],.5)
            dtype('float64')


        """
        if(F.m_faust.isReal()):
            return np.dtype(np.float64)
        else:
            return np.dtype(np.complex)

class FaustFactory:
    """
    This factory class provides methods for generating a Faust especially by factorization of a dense matrix.

    This class gives access to the main factorization algorithms of
    FAµST. Those algorithms can factorize a dense matrix to a sparse product
    (i.e. a Faust object).

    There are two algorithms for factorization.

    The first one is Palm4MSA :
    which stands for Proximal Alternating Linearized Minimization for
    Multi-layer Sparse Approximation. Note that Palm4MSA is not
    intended to be used directly. You should rather rely on the second algorithm.

    The second one is the Hierarchical Factorization algorithm:
    this is the central algorithm to factorize a dense matrix to a Faust.
    It makes iterative use of Palm4MSA to proceed with the factorization of a given
    dense matrix.

    A more secondary functionality of this class is the pseudo-random generation of a
    Faust with FaustFactory.rand().

    """

    @staticmethod
    def fact_palm4msa(M, p):
        """
        Factorizes the matrix M with Palm4MSA algorithm using the parameters set in p.

        Args:
            M: the numpy matrix to factorize.
            p: the ParamsPalm4MSA instance to define the algorithm parameters.

        Returns:
            The Faust object result of the factorization.

        Examples:
        >>> from pyfaust import FaustFactory, ParamsPalm4MSA, ConstraintReal,\
        >>> ConstraintInt, ConstraintName, StoppingCriterion
        >>> import numpy as np
        >>> num_facts = 2
        >>> is_update_way_R2L = False
        >>> init_lambda = 1.0
        >>> M = np.random.rand(500, 32)
        >>> cons1 = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5)
        >>> cons2 = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32, 32, 1.0)
        >>> stop_crit = StoppingCriterion(num_its=200)
        >>> # default step_size is 1e-16
        >>> param = ParamsPalm4MSA(num_facts, is_update_way_R2L, init_lambda,
        >>>                        [cons1, cons2], stop_crit,
        >>>                        is_verbose=False)
        >>> F = FaustFactory.fact_palm4msa(M, param)
        >>> F.display()
        Faust size 500x32, density 0.22025, nnz_sum 3524, 2 factor(s):<br/>
        FACTOR 0 (real) SPARSE, size 500x32, density 0.15625, nnz 2500<br/>
        FACTOR 1 (real) SPARSE, size 32x32, density 1, nnz 1024<br/>
        >>> # verify if NORMCOL constraint is respected by the 2nd factor
        >>> from numpy.linalg import norm
        >>> for i in range(0, F.get_factor(1).shape[1]):
        >>>     if(i<F.get_factor(1).shape[1]-1):
        >>>             end=' '
        >>>     else:
        >>>             end='\n'
        >>>     print(norm(F.get_factor(1)[:,i]), end=end)
        1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
        """
        if(not isinstance(p, ParamsPalm4MSA)):
            raise ValueError("p must be a ParamsPalm4MSA object.")
        FaustFactory._check_fact_mat('FaustFactory.fact_palm4msa()', M)
        return Faust(core_obj=FaustCorePy.FaustFact.fact_palm4MSA(M, p))

    @staticmethod
    def fact_hierarchical(M, p):
        """
        Factorizes the matrix M with Hierarchical Factorization using the parameters set in p.

        Args:
            M: the numpy matrix to factorize.
            p: the ParamsHierarchicalFact instance to define the algorithm parameters.

        Returns:
            The Faust object result of the factorization.

        Examples:
            >>> from pyfaust import FaustFactory, ParamsHierarchicalFact, ConstraintReal,\
            >>> ConstraintInt, ConstraintName, StoppingCriterion
            >>> import numpy as np
            >>> num_facts = 4
            >>> is_update_way_R2L = False
            >>> init_lambda = 1.0
            >>> M = np.random.rand(500, 32)
            >>> fact0_cons = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5)
            >>> fact1_cons = ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 96)
            >>> fact2_cons = ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 96)
            >>> res0_cons = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32, 32, 1)
            >>> res1_cons =  ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 666)
            >>> res2_cons =  ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 333)
            >>> stop_crit1 = StoppingCriterion(num_its=200)
            >>> stop_crit2 = StoppingCriterion(num_its=200)
            >>> param = ParamsHierarchicalFact(num_facts, is_update_way_R2L, init_lambda,
            >>>                                [fact0_cons, fact1_cons, fact2_cons],
            >>>                                [res0_cons, res1_cons, res2_cons],
            >>>                                M.shape[0],M.shape[1],[stop_crit1,
            >>>                                                       stop_crit2],
            >>>                                is_verbose=False)
            >>> F = FaustFactory.fact_hierarchical(M, param)
            Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorisation
            1/3<br/>
            Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorisation
            2/3<br/>
            Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorisation
            3/3<br/>
            >>> F.display()
            Faust size 500x32, density 0.189063, nnz_sum 3025, 4 factor(s):
                - FACTOR 0 (real) SPARSE, size 500x32, density 0.15625, nnz 2500
                - FACTOR 1 (real) SPARSE, size 32x32, density 0.09375, nnz 96
                - FACTOR 2 (real) SPARSE, size 32x32, density 0.09375, nnz 96
                - FACTOR 3 (real) SPARSE, size 32x32, density 0.325195, nnz 333
                """
        if(not isinstance(p, ParamsHierarchicalFact)):
            raise ValueError("p must be a ParamsHierarchicalFact object.")
        FaustFactory._check_fact_mat('FaustFactory.fact_hierarchical()', M)
        return Faust(core_obj=FaustCorePy.FaustFact.fact_hierarchical(M, p))

    @staticmethod
    def hadamard(n):
        """
           Constructs a Faust whose the full matrix is the Hadamard matrix of size 2^n*2^n.
        """
        H = Faust(core_obj=FaustCorePy.FaustCore.hadamardFaust(n))
        return H

    @staticmethod
    def fourier(n):
        """
            Constructs a Faust whose the full matrix is the FFT matrix of size 2^n*2^n.
        """
        F = Faust(core_obj=FaustCorePy.FaustCore.fourierFaust(n))
        return F

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
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand(2, 10, .5, is_real=False)
            >>> G = FaustFactory.rand([2, 5], [10, 20], .5, fac_type="dense")
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
            raise ValueError('FaustFactory.rand(): argument fac_type must be a'
                             ' str equal to "sparse",'
                             ' "dense" or "mixed".')

        fac_type_map = {"sparse" : SPARSE, "dense" : DENSE, "mixed" : MIXED}
        # set field of matrices/factors
        if(not isinstance(is_real, bool)):
            raise ValueError('FaustFactory.rand(): argument is_real must be a'
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
            raise ValueError("FaustFactory.rand(): num_factors argument must be an "
                             "integer or a list/tuple of 2 integers.")
        if((isinstance(dim_sizes, list) or isinstance(dim_sizes, tuple)) and
        len(dim_sizes) == 2):
            min_dim_size = dim_sizes[0]
            max_dim_size = dim_sizes[1]
        elif(isinstance(dim_sizes, int)):
            min_dim_size = max_dim_size = dim_sizes
        else:
            raise ValueError("FaustFactory.rand(): dim_sizes argument must be an "
                             "integer or a list/tuple of 2 integers.")
        rF = Faust(core_obj=FaustCorePy.FaustCore.randFaust(
            fac_type_map[fac_type], field, min_num_factors, max_num_factors,
            min_dim_size, max_dim_size, density))
        return rF

    @staticmethod
    def _check_fact_mat(funcname, M):
        if(not isinstance(M, np.ndarray)):
            raise Exception(funcname+" 1st argument must be a numpy ndarray.")
        if(not isinstance(M[0,0], np.complex) and not isinstance(M[0,0],
                                                                np.float)):
            raise Exception(funcname+" 1st argument must be a float or complex "
                            "ndarray.")
        #if(isinstance(M[0,0], np.complex)):
        #   raise Exception(funcname+" doesn't yet support complex matrix "
        #                   "factorization.")



class ParamsFact(object):

    def __init__(self, num_facts, is_update_way_R2L, init_lambda,
                 constraints, step_size, constant_step_size=False, is_verbose=False,
                 ):
        self.num_facts = num_facts
        self.is_update_way_R2L = is_update_way_R2L
        self.init_lambda = init_lambda
        self.step_size = step_size
        self.constraints = constraints
        self.is_verbose = is_verbose
        self.constant_step_size = constant_step_size

class ParamsPalm4MSA(ParamsFact):

    def __init__(self, num_facts, is_update_way_R2L, init_lambda,
                 constraints, stop_crit, init_facts=None, step_size=10.0**-16,
                 constant_step_size=False,
                 is_verbose=False):
        super(ParamsPalm4MSA, self).__init__(num_facts, is_update_way_R2L,
                                             init_lambda,
                                             constraints, step_size,
                                             constant_step_size,
                                             is_verbose)
        if(init_facts != None and (not isinstance(init_facts, list) and not isinstance(init_facts,
                                                               tuple) or
           len(init_facts) != num_facts)):
            raise ValueError('ParamsPalm4MSA init_facts argument must be a '
                             'list/tuple of '+str(num_facts)+" (num_facts) arguments.")
        else:
            self.init_facts = init_facts
        if(not isinstance(stop_crit, StoppingCriterion)):
           raise TypeError('ParamsPalm4MSA stop_crit argument must be a StoppingCriterion '
                           'object')
        self.stop_crit = stop_crit
        #TODO: verify number of constraints is consistent with num_facts

class ParamsHierarchicalFact(ParamsFact):

    def __init__(self, num_facts, is_update_way_R2L, init_lambda,
                 fact_constraints, res_constraints, data_num_rows,
                 data_num_cols, stop_crits,
                 step_size=10.0**-16, constant_step_size=False,
                 is_verbose=False,
                 is_fact_side_left = False):
        constraints = fact_constraints + res_constraints
        super(ParamsHierarchicalFact, self).__init__(num_facts,
                                                     is_update_way_R2L,
                                                     init_lambda,
                                                     constraints, step_size,
                                                     constant_step_size,
                                                     is_verbose)
        self.data_num_rows = data_num_rows
        self.data_num_cols = data_num_cols
        self.stop_crits = stop_crits
        self.is_fact_side_left = is_fact_side_left
        #TODO: verify number of constraints is consistent with num_facts in
        if((not isinstance(stop_crits, list) and not isinstance(stop_crits,
                                                                tuple)) or
           len(stop_crits) != 2 or
           not isinstance(stop_crits[0],StoppingCriterion) or not
           isinstance(stop_crits[1],StoppingCriterion)):
            raise TypeError('ParamsHierarchicalFact stop_crits argument must be a list/tuple of two '
                            'StoppingCriterion objects')
        if((not isinstance(constraints, list) and not isinstance(constraints,
                                                                tuple)) or
           np.array([not isinstance(constraints[i],ConstraintGeneric) for i in
                    range(0,len(constraints))]).any()):
            raise TypeError('constraints argument must be a list/tuple of '
                            'ConstraintGeneric (or subclasses) objects')

class StoppingCriterion(object):

    def __init__(self, is_criterion_error = False , error_treshold = 0.3,
                 num_its = 500,
                 max_num_its = 1000):
        self.is_criterion_error = is_criterion_error
        self.error_treshold = error_treshold
        self.num_its = num_its
        self.max_num_its = max_num_its
        #TODO: check_validity() like C++ code does
        if(is_criterion_error and num_its != 500):
            raise ValueError("It's forbidden to set a number of iterations as stopping"
                             " criterion when is_criterion_error == True.")
        elif(not is_criterion_error and (max_num_its != 1000 or error_treshold
                                         != 0.3)):
            raise ValueError("When is_criterion_error == True it's forbidden to use"
                             " other arguments than num_its argument to define "
                             "the stopping criterion.")

class ConstraintName:
    """
    Attributes:
        SP: Designates a constraint on the sparsity/0-norm of a matrix.
        SPCOL: Designates a sparsity/0-norm constraint on the columns of a
        matrix.
        SPLIN: Designates a sparsity/0-norm constraint on the lines of a
        matrix.
        SPLINCOL: Designates a constraint that imposes both SPLIN and SPCOL
        constraints.
        SP_POS: Designates a constraint that imposes a SP constraints and
        besides erase the negative coefficients (it doesn't apply to complex
        matrices).
        NORMCOL: Designates a 2-norm constraint on the columns of a matrix.
        NORMLIN: Designates a 2-norm constraint on the lines of a matrix.
        CONST: Designates a constraint imposing to a matrix to be constant.

        
    """
    SP = 0 # Int Constraint
    SPCOL = 1 # Int Constraint
    SPLIN=2 # Int Constraint
    NORMCOL = 3 # Real Constraint
    SPLINCOL = 4 # Int Constraint
    CONST = 5 # Mat Constraint
    SP_POS = 6 # Int Constraint
    #BLKDIAG = 7 # ?? Constraint #TODO
    SUPP = 8 # Mat Constraint
    NORMLIN = 9 # Real Constraint

    def __init__(self, name):
        if(not isinstance(name, np.int) or name < ConstraintName.SP or name > ConstraintName.NORMLIN):
            raise ValueError("name must be an integer among ConstraintName.SP,"
                             "ConstraintName.SPCOL, ConstraintName.NORMCOL,"
                             "ConstraintName.SPLINCOL, ConstraintName.CONST,"
                             "ConstraintName.SP_POS," # ConstraintName.BLKDIAG,
                             "ConstraintName.SUPP, ConstraintName.NORMLIN")
        self.name = name

    def is_int_constraint(self):
        return self.name in [ ConstraintName.SP, ConstraintName.SPCOL,
                             ConstraintName.SPLIN, ConstraintName.SPLINCOL,
                             ConstraintName.SP_POS ]

    def is_real_constraint(self):
        return self.name in [ ConstraintName.NORMCOL, ConstraintName.NORMLIN ]

    def is_mat_constraint(self):
        return self.name in [ConstraintName.SUPP, ConstraintName.CONST ]

class ConstraintGeneric(object):
    """
    This is the parent class for representing a factor constraint in FAµST factorization algorithms.

    This class shouldn't be instantiated, rather rely on sub-classes.

    <b/> See also: ConstraintInt, ConstraintReal, ConstraintMat, FaustFactory.fact_palm4msa, FaustFactory.fact_hierarchical.

    Attributes:
        _name: The name of the constraint applied to the factor (ConstraintName instance).
        _num_rows: The number of rows to constrain the factor to.
        _num_cols: The number of columns to constrain the factor to.
        _cons_value: The parameter value of the constraint.

    """

    def __init__(self, name, num_rows, num_cols, cons_value):
        """
        Constructs a generic constraint.

        Args:
            name: the name of the constraint applied to the factor (ConstraintName instance).
            num_rows: the number of rows to constrain the factor to.
            num_cols: the number of columns to constrain the factor to.
            cons_value: the parameter value of the constraint.

        """
        self._name = name
        self._num_rows = num_rows
        self._num_cols = num_cols
        self._cons_value = cons_value

    @property
    def name(self):
        """
            Property to access the ConstraintName of the constraint.
        """
        return self._name.name

    def is_int_constraint(self):
        """
            Returns True if this constraint is a ConstraintInt, False otherwise.
        """
        return self._name.is_int_constraint()

    def is_real_constraint(self):
        """
            Returns True if this constraint is a ConstraintReal, False otherwise.
        """
        return self._name.is_real_constraint()

    def is_mat_constraint(self):
        """
            Returns True if this constraint is a ConstraintMat, False otherwise.
        """
        return self._name.is_mat_constraint()

class ConstraintReal(ConstraintGeneric):
    """
        This class represents a real constraint on a matrix-factor.

        It constrains a matrix by a column/row-vector 2-norm
        (ConstraintName.NORMCOL, ConstraintName.NORMLIN).

    """
    def __init__(self, name, num_rows, num_cols, cons_value):
        """
        Constructs a real type constraint.

        Args:
            name: the name of the constraint applied to the factor. It must be
            a ConstraintName instance.
            num_rows: the number of rows to constrain the factor to.
            num_cols: the number of columns to constrain the factor to.
            cons_value: the parameter value of the constraint, it must be a
            float number that designates the 2-norm imposed to all columns (if
            name.name == ConstraintName.NORMCOL) or rows (if
            name.name == ConstraintName.NORMLIN).
        """
        super(ConstraintReal, self).__init__(name, num_rows, num_cols, cons_value)
        if(not isinstance(cons_value, np.float) and not isinstance(cons_value, np.int)):
            raise TypeError('ConstraintReal must receive a float as cons_value '
                            'argument.')
        self._cons_value = float(self._cons_value)
        if(not isinstance(name, ConstraintName) or not name.is_real_constraint()):
            raise TypeError('ConstraintReal first argument must be a '
                            'ConstraintName with a real type name '
                            '(name.is_real_constraint() must return True).')

class ConstraintInt(ConstraintGeneric):
    """
        This class represents an integer constraint on a matrix-factor.

        It constrains a matrix by its column/row-vectors sparsity or 0-norm
        (ConstraintName.SPLIN, ConstraintName.SPCOL, ConstraintName.SPLINCOL).

    """
    def __init__(self, name, num_rows, num_cols, cons_value):
        super(ConstraintInt, self).__init__(name, num_rows, num_cols, cons_value)
        if(not isinstance(cons_value, np.int)):
            raise TypeError('ConstraintInt must receive a int as cons_value '
                            'argument.')
        if(not isinstance(name, ConstraintName) or not name.is_int_constraint()):
            raise TypeError('ConstraintInt first argument must be a '
                            'ConstraintName with a int type name '
                            '(name.is_int_constraint() must return True).')

class ConstraintMat(ConstraintGeneric):

    def __init__(self, name, num_rows, num_cols, cons_value):
        super(ConstraintMat, self).__init__(name, num_rows, num_cols, cons_value)
        if(not isinstance(cons_value, np.matrix) and not isinstance(cons_value,
                                                               np.ndarray)):
            raise TypeError('ConstraintMat must receive a numpy matrix as cons_value '
                            'argument.')
        self.cons_value = float(self.cons_value)
        if(not isinstance(name, ConstraintName) or not name.is_mat_constraint()):
            raise TypeError('ConstraintMat first argument must be a '
                            'ConstraintName with a matrix type name'
                            '(name.is_mat_constraint() must return True).')

