# -*- coding: utf-8 -*-
#                              Description:                                ##
#                                                                          ##
#          pyfaust is a python module  which delivers                      ##
#          a class named Faust which represents a dense matrix             ##
#          by a product of 'sparse' factors (i.e Faust)                    ##
#          Python wrapper class implemented in C++                         ##
#                                                                          ##
#  For more information on the FAuST Project, please visit the website     ##
#  of the project : <http://faust.inria.fr>                                ##
#                                                                          ##
#                              License:                                    ##
#  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      ##
#                      Luc Le Magoarou, Remi Gribonval                     ##
#                      INRIA Rennes, FRANCE                                ##
#                      http://www.inria.fr/                                ##
#                                                                          ##
#  The FAuST Toolbox is distributed under the terms of the GNU Affero      ##
#  General Public License.                                                 ##
#  This program is free software: you can redistribute it and/or modify    ##
#  it under the terms of the GNU Affero General Public License as          ##
#  published by the Free Software Foundation.                              ##
#                                                                          ##
#  This program is distributed in the hope that it will be useful, but     ##
#  WITHOUT ANY WARRANTY; without even the implied warranty of              ##
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    ##
#  See the GNU Affero General Public License for more details.             ##
#                                                                          ##
#  You should have received a copy of the GNU Affero General Public        ##
#  License along with this program.                                        ##
#  If not, see <http://www.gnu.org/licenses/>.                             ##
#                                                                          ##
#                             Contacts:                                    ##
#      Nicolas Bellot  : nicolas.bellot@inria.fr                           ##
#      Adrien Leman    : adrien.leman@inria.fr                             ##
#      Thomas Gautrais : thomas.gautrais@inria.fr                          ##
#      Luc Le Magoarou : luc.le-magoarou@inria.fr                          ##
#      Remi Gribonval  : remi.gribonval@inria.fr                           ##
#                                                                          ##
#############################################################################

## @package pyfaust @brief <b> The FAµST Python Wrapper</b>
import copy

import numpy as np, scipy
from scipy.io import savemat, loadmat
from scipy.sparse import csr_matrix, csc_matrix
import FaustCorePy
import pyfaust

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

    Although sparse matrices are more interesting for optimization it's not
    forbidden to define a Faust as a product of dense matrices or a mix of dense
    and sparse matrices.

    The matrices composing the Faust product, also called the factors, are
    defined on complex or real fields. Hence a Faust can be a complex Faust or a
    real Faust.

    Several Python builtins have been overloaded to ensure that a Faust is
    almost handled as a native numpy array/matrix.

    The main exception is that contrary to a numpy matrix a Faust is immutable.
    It means that you cannot affect elements of a Faust using
    the affectation operator `=' like you do with a numpy matrix (e.g. `M[i,j] =
    2').
    That limitation is the reason why the Python built-in `__setitem__()' is not
    implemented in this class.

    Other notable limitations are that one cannot:
        - compute the real and imaginary parts of a Faust,
        - perform elementwise operations between two Fausts (e.g. elementwise
        multiplication),
        - reshape.

    Primarily for convenience and test purposes, a Faust can be converted into
    the corresponding full matrix using the function Faust.todense or
    Faust.toarray.

    Warning: using Faust.todense or Faust.toarray is discouraged except for test purposes, as it
    loses the main potential interests of the FAuST structure: compressed
    memory storage and faster matrix-vector multiplication compared to its
    equivalent full matrix representation.

    In this documentation, the expression 'full matrix' designates the array
    Faust.toarray() obtained by the multiplication of the Faust factors.

   List of functions that are memory costly: Faust.toarray(), Faust.todense(),
   Faust.pinv().

    For more information about FAµST take a look at http://faust.inria.fr.
    """

    def  __init__(F, factors=None, filepath=None, **kwargs):
        """ Creates a Faust from a list of factors or alternatively from a file.

            Other easy ways to create a Faust is to call one of the
            FaustFactory static methods: FaustFactory.rand(),
            FaustFactory.dft() or FaustFactory.wht().

        Args:
            factors: list of numpy/scipy array/matrices or a single array/matrix.<br/>
                     The factors must respect the dimensions needed for
                     the product to be defined (for i=0 to len(factors)-1,
                     factors[i].shape[1] == factors[i+1].shape[0]).<br/>
                     The factors can be sparse or dense matrices
                     (either scipy.sparse.csr.csr_matrix or
                     numpy.ndarray).<br/>
                     Passing only an array or or sparse matrix to the
                     constructor is equivalent to
                     passing a list of a single factor.
            filepath: the file from which a Faust is created.<br/>
                      The format is Matlab version 5 (.mat extension).<br/>
                      The file must have been saved before with Faust.save().
            *kwargs: internal purpose arguments.

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

            >>> F.save("F.mat")
            >>> # define a Faust from file
            >>> H = Faust(filepath="F.mat")

            >>> Faust(np.random.rand(10,10)) # creating a Faust with only one
                                             # factor

        <b/> See also Faust.save, FaustFactory.rand

        """
        if("scale" in kwargs.keys()):
            # scale hidden argument
            scale = kwargs['scale']
            if(not (isinstance(scale, float) or isinstance(scale, int) or
               isinstance(scale, np.complex))):
                raise Exception("Scale must be a number.")
        else:
            scale = 1.0
        if("core_obj" in kwargs.keys()):
            core_obj = kwargs['core_obj']
            if(core_obj):
                F.m_faust = core_obj
        else:
            if(filepath and isinstance(filepath, str)):
                    contents = loadmat(filepath)
                    factors = contents['faust_factors'][0]
            if((isinstance(factors, np.ndarray) and factors.ndim == 2)
               or isinstance(factors,
                             scipy.sparse.csc.csc_matrix)
               or isinstance(factors, scipy.sparse.csr.csr_matrix)):
                factors = [ factors ]
            if(not isinstance(factors, list) and not
               isinstance(factors, np.ndarray)):
                raise Exception("factors must be a non-empty list of/or a numpy.ndarray, "
                                "scipy.sparse.csr.csr_matrix/csc.csc_matrix.")
            if(factors is not None and len(factors) > 0):
                F.m_faust = FaustCorePy.FaustCore(factors, scale);
            else:
                raise Exception("Cannot create an empty Faust.")

    @property
    def shape(F):
        """
        Gives the shape of the Faust F.

        This function is intended to be used as a property (see the examples).

        The shape is a pair of numbers: the number of rows and the number of
        columns of F.todense().

        Args:
            F: the Faust object.

        Returns:
            the Faust shape tuple, with at first index the number of rows, and
            at second index the number of columns.

        Examples:
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand(2, 50, field='complex')
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
            >>> F = FaustFactory.rand(2, 50, field='complex')
            >>> size = F.size

        <b/> See also Faust.shape

        """
        return np.prod(F.shape)

    def transpose(F):
        """
        Returns the transpose of F.

        Args:
            F: the Faust object.

        Returns:
            a Faust object implementing the transpose of F.todense(), such
            that:
            <code>F.transpose().todense() == F.todense().transpose()</code>

        Examples:
            >>> tF = F.transpose()

        <b/> See also Faust.conj, Faust.getH, Faust.H, Faust.T
        """
        F_trans = Faust(core_obj=F.m_faust.transpose())
        return F_trans

    @property
    def T(F):
        """
        Returns the transpose of F.

        This function is intended to be used as a property (see the examples).

        Args:
            F: the Faust object.

        Returns:
            a Faust object implementing the transpose of F.todense(), such
            that:
            <code>F.T.todense() == F.todense().T</code>

        Examples:
            >>> tF = F.T

        <b/> See also Faust.conj, Faust.getH, Faust.H, Faust.T
        """
        return F.transpose()

    def conj(F):
        """
        Returns the complex conjugate of F.

        Args:
            F: the Faust object.

        Returns:
            a Faust object Fc implementing the conjugate of F.todense() such
            that:
            <code>Fc.todense() == F.todense().conj()</code>


        Examples:
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand(5, 50, field='complex')
            >>> Fc = F.conj()

        <b/> See also Faust.transpose, Faust.getH, Faust.H
        """
        F_conj = Faust(core_obj=F.m_faust.conj())
        return F_conj

    def getH(F):
        """
        Returns the conjugate transpose of F.

        Args:
            F: the Faust object.

        Returns:
            a Faust object H implementing the conjugate transpose of
            F.todense() such that:
            <code>H.todense() == F.todense().getH()</code>

        Examples:
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand(5, 50, field='complex')
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
            a Faust object H implementing the conjugate transpose of
            F.todense() such that:
            <code>H.todense() == F.todense().H</code>

        Examples:
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand(5, 50, field='complex')
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

    def __add__(F, *args):
        """
        Sums F to one or a sequence of variables. Faust objects, arrays or scalars.

        This method overloads the Python function/operator +.

        Args:
            *args: the list of variables to sum all together with F. These can
            be: Faust objects, arrays (and scipy.sparse.csr_matrix) or scalars.
            Faust and arrays/matrices must be of compatible size.

        Returns: the sum as a Faust object.

        <b/>See also Faust.__sub__
        """
        for i in range(0,len(args)):
            G = args[i]
            if(isinstance(G,Faust)):
                if(F.shape != G.shape):
                    raise Exception('Dimensions must agree')
                C = F.concatenate(G, axis=1)
                # Id = np.eye(int(C.shape[1]/2))
                Fid = FaustFactory.eye(int(C.shape[1]/2))
                #F = C*Faust(np.concatenate((Id,Id)),axis=0)
                F = C*Fid.concatenate(Fid,axis=0)
            elif(isinstance(G,np.ndarray) or
                 isinstance(G,scipy.sparse.csr_matrix) or
                 isinstance(G,scipy.sparse.csc_matrix)):
                F = F+Faust(G)
            elif(isinstance(G,int) or isinstance(G,float) or isinstance(G, np.complex)):
                if(F.shape[0] <= F.shape[1]):
                    F = F+Faust([np.eye(F.shape[0], F.shape[1]),
                                 np.ones((F.shape[1], 1))*G,
                                 np.ones((1, F.shape[1]))])
                else:
                    F = F+Faust([np.eye(F.shape[1], F.shape[0]),
                                 np.ones((F.shape[0], 1))*G,
                                 np.ones((1, F.shape[0]))]).T
            else:
                raise Exception("Cannot add a Faust to something that is not a"
                                " Faust, a matrix/array or a scalar.")
        return F

    def __sub__(F,*args):
        """
        Substracts from F one or a sequence of variables. Faust objects, arrays or scalars.

        This method overloads the Python function/operator -.

        Args:
            *args: the list of variables to compute the difference with F. These can
            be: Faust objects, arrays (and scipy.sparse.csr_matrix) or scalars.
            Faust and arrays/matrices must be of compatible size.

        Returns: the difference as a Faust object.

        <b/>See also Faust.__add__
        """

        nargs = []
        for arg in args:
            nargs += [ arg*-1 ]
        return F.__add__(*nargs)

    def __truediv__(F,s):
        """
        Divides F by the scalar s.

        This method overloads the Python function/operator `/' (whether s is a
        float or an integer).

        Args:
        F: the Faust object.
        s: the scalar to divide the Faust object with.

        Returns:
            the division result as a Faust object.

        <b/> See also Faust.__mul__
        """
        if(isinstance(s, float) or isinstance(s, np.complex) or isinstance(s,
                                                                           int)):
            return F*(1./s)
        else:
            raise Exception("unsupported operand type(s) for /: a Faust can only be "
                  "divided by a scalar.")


    def __mul__(F, A):
        """
        Multiplies F by A which is a full matrix, a Faust object or a scalar number.

        This method overloads a Python function/operator.

        <b>The primary goal</b> of this function is to implement “fast” multiplication by a
        Faust, which is the operation performed when A is a dense matrix.<br/>
        In the best case, F*A is F.rcg() times faster than equivalent F.toarray()*A.

        <b>Other use cases</b> are available for this function:
        - If A is a Faust, no actual multiplication is performed, instead a
        new Faust is built to implement the multiplication.<br/>
        This Faust verifies that:<br/>
            <code>
            (F*A).todense() == F.todense()*A.todense()
            </code>
            <br/>N.B.: you could have an elementwise non-significant absolute
            difference between the two members.
        - If A is a scalar, F*A is also a Faust such that:<br/>
        <code>
        (F*A).get_factor(0) ==  F.get_factor(0)*A
        </code>

        Args:
            F: the Faust object.
            A: is a scalar number, a Faust object or a 2D full matrix (numpy.ndarray,
            numpy.matrix).
            <br/> In the latter case, A must be Fortran contiguous (i.e. Column-major order;
                `order' argument passed to np.ndararray() must be equal to str
                'F').

        Returns: The result of the multiplication
            - as a numpy.ndarray if A is a ndarray,<br/>
            - as a Faust object if A is a Faust
            or a scalar number.
            - When either F or A is complex, G=F*A is also complex.

        Raises:
            ValueError
            The multiplying operand A is a sparse matrix:
            <code>
            >>> from scipy.sparse import csr_matrix
            >>> from scipy.sparse import random as srand
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand(5,50)
            >>> S = srand(50,5000,.1)
            >>> F*S
            ValueError input M must a numpy.ndarray or a numpy.matrix.
            </code>
            ValueError
            F is real but A is a complex scalar:
            <code>
            >>> import numpy
            >>> F*numpy.complex(0,1)
            ValueError You cannot multiply a real Faust by a complex scalar (not yet implemented).
            </code>

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
            >>> F_times_two = F*2

        <b/>See also Faust.__init__, Faust.rcg, Faust.__mul__
        """
        if(isinstance(A, Faust)):
            if(F.shape[1] != A.shape[0]): raise ValueError("The dimensions of "
                                                          "the two Fausts must "
                                                          "agree.")
            return Faust(core_obj=F.m_faust.multiply(A.m_faust))
        elif(isinstance(A, float) or isinstance(A, int) or isinstance(A, np.complex)):
            return Faust(core_obj=F.m_faust.multiply_scal(A))
        elif(isinstance(A, np.ndarray) and isinstance(A[0,0], np.complex)):
            j = np.complex(0,1)
            return F.m_faust.multiply(A.real).astype(np.complex) + j*F.m_faust.multiply(A.imag)
        else:
            return F.m_faust.multiply(A)

    #def concatenate(F, *args, axis=0): # py. 2 doesn't handle this signature
    def concatenate(F, *args, **kwargs):
        """
            Concatenates F with len(args) Faust objects, numpy arrays or scipy sparse matrices.

            The resulting Faust <code>C = F.concatenate(G, H, ... ,axis)</code> verifies that:
            <code>
            C.toarray() == numpy.concatenate((F.toarray(), G.toarray(),
            H.toarray(),...), axis)
            </code>
            <br/>N.B.: you could have an elementwise non-significant absolute
            difference between the two members.


           Args:
               F: the Faust to concatenate to.
               args: the Fausts or matrices (numpy array or
               scipy.sparse.csr/csc_matrix) to be concatenated to F. If args[i] is a
               matrix it will be Faust-converted on the fly.
               axis (optional): the dimension index (0 or 1) along to concatenate the
               Faust objects. By default, the axis is 0 (for vertical
               concatenation).

            Returns:
                The concatenation result as a Faust.

            Raises:
                ValueError
                The dimensions of the two Fausts must agree.
                <code>
                >>> F = FaustFactory.rand(2,51);
                >>> G = FaustFactory.rand(2,25);
                >>> F.concatenate(G, 0)
                </code>
                ValueError
                Axis must be 0 or 1.
                <code>
                >>> F = FaustFactory.rand(2,51);
                >>> G = FaustFactory.rand(2,25);
                >>> F.concatenate(G, 5)
                </code>
                ValueError
                You can't concatenate a Faust with something that is not a
                Faust, a numpy array or a scipy sparse matrix.
                <code>
                >>> F = FaustFactory.rand(2,51);
                >>> F.concatenate(['a', 'b', 'c'], 5)
                </code>

            Examples:
                >>> import pyfaust import FaustFactory
                >>> F = FaustFactory.rand(5, 50)
                >>> G = FaustFactory.rand(6, 50)
                >>> F.concatenate(G) # equivalent to F.concatenate(G, 0)
                Faust size 100x50, density 0.5634, nnz_sum 2817, 7 factor(s)<br/>
                FACTOR 0 (real) SPARSE, size 100x100, density 0.0473, nnz 47<br/>
                FACTOR 1 (real) SPARSE, size 100x100, density 0.0469, nnz 46<br/>
                FACTOR 2 (real) SPARSE, size 100x100, density 0.0504, nnz 50<br/>
                FACTOR 3 (real) SPARSE, size 100x100, density 0.0502, nnz 50<br/>
                FACTOR 4 (real) SPARSE, size 100x100, density 0.0482, nnz 48<br/>
                FACTOR 5 (real) SPARSE, size 100x100, density 0.0287, nnz 28<br/>
                FACTOR 6 (real) SPARSE, size 100x50, density 0.02, nnz 10<br/>
                >>> F.concatenate(G, 1)
                Faust size 50x100, density 0.5634, nnz_sum 2817, 7 factor(s)<br/>
                FACTOR 0 (real) SPARSE, size 50x100, density 0.02, nnz 10<br/>
                FACTOR 1 (real) SPARSE, size 100x100, density 0.0286, nnz 28<br/>
                FACTOR 2 (real) SPARSE, size 100x100, density 0.0469, nnz 46<br/>
                FACTOR 3 (real) SPARSE, size 100x100, density 0.0504, nnz 50<br/>
                FACTOR 4 (real) SPARSE, size 100x100, density 0.0504, nnz 50<br/>
                FACTOR 5 (real) SPARSE, size 100x100, density 0.0479, nnz 47<br/>
                FACTOR 6 (real) SPARSE, size 100x100, density 0.0475, nnz 47<br/>
                >>> from numpy.random import rand
                >>> F.concatenate(rand(34, 50), axis=0) # The random array is auto-converted to a Faust before the vertical concatenation
                Faust size 84x50, density 0.765476, nnz_sum 3215, 6 factor(s):<br/>
                FACTOR 0 (real) SPARSE, size 84x100, density 0.23119, nnz 1942<br/>
                FACTOR 1 (real) SPARSE, size 100x100, density 0.0301, nnz 301<br/>
                FACTOR 2 (real) SPARSE, size 100x100, density 0.0286, nnz 286<br/>
                FACTOR 3 (real) SPARSE, size 100x100, density 0.0285, nnz 285<br/>
                FACTOR 4 (real) SPARSE, size 100x100, density 0.0301, nnz 301<br/>
                FACTOR 5 (real) SPARSE, size 100x50, density 0.02, nnz 100<br/>
                >>> from scipy.sparse import rand as sprand
                >>> F.concatenate(sprand(50, 24, format='csr'), axis=1) # The sparse random matrix is auto-converted to a Faust before the horizontal concatenation
                Faust size 50x74, density 0.412703, nnz_sum 1527, 6 factor(s):<br/>
                FACTOR 0 (real) SPARSE, size 50x100, density 0.02, nnz 100<br/>
                FACTOR 1 (real) SPARSE, size 100x100, density 0.0292, nnz 292<br/>
                FACTOR 2 (real) SPARSE, size 100x100, density 0.0301, nnz 301<br/>
                FACTOR 3 (real) SPARSE, size 100x100, density 0.0286, nnz 286<br/>
                FACTOR 4 (real) SPARSE, size 100x100, density 0.0285, nnz 285<br/>
                FACTOR 5 (real) SPARSE, size 100x74, density 0.0355405, nnz 263<br/>
                >>> F.concatenate(F, G, F, G, rand(34,50), F, G) # it's allowed to concatenate an arbitrary number of Fausts
        """
        if("axis" in kwargs.keys()):
            axis = kwargs['axis']
        else:
            axis = 0
        if(axis not in [0,1]):
            raise ValueError("Axis must be 0 or 1.")

        for G in args:
            if(axis == 0 and F.shape[1] != G.shape[1] or axis == 1 and F.shape[0]
               != G.shape[0]): raise ValueError("The dimensions of "
                                                "the two Fausts must "
                                                "agree.")

        C=F
        for G in args:
           if(isinstance(G, np.ndarray) or isinstance(G, csr_matrix) \
               or isinstance(G, csc_matrix)):
               G = Faust([G])
           if(not isinstance(G, Faust)): raise ValueError("You can't concatenate a "
                                                           "Faust with something "
                                                           "that is not a Faust, "
                                                           "a numpy array or scipy "
                                                           "sparse matrix.")
           if(axis==0):
                C = Faust(core_obj=C.m_faust.vertcat(G.m_faust))
           elif(axis==1):
                C = Faust(core_obj=C.m_faust.horzcat(G.m_faust))
        return C

    def toarray(F):
        """
        Returns a numpy array for the full matrix implemented by F.

        WARNING: Using this function is discouraged except for test purposes,
        as it loses the main potential interests of the Faust structure:
        compressed memory storage and faster matrix-vector multiplication
        compared to its equivalent full matrix representation.

        Returns:
            A numpy ndarray.

        Raises:
            MemoryError

        WARNING: running the example below is likely to raise a memory
        error or freeze your computer for a certain amount of time.

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
            >>> # the sparse format is the only way to handle such a large Faust
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
        Returns a numpy matrix for the full matrix implemented by F.

        WARNING: Using this function is discouraged except for test purposes,
        as it loses the main potential interests of the Faust structure:
        compressed memory storage and faster matrix-vector multiplication
        compared to its equivalent full matrix representation.

        Returns:
            A numpy matrix M such that M*x == F*x for any vector x.

        Raises:
            MemoryError

        WARNING: running the example below is likely to raise a memory
        error or freeze your computer for a certain amount of time.

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
            >>> # the sparse format is the only way to handle such a large Faust
            >>> F.todense()
            (...)
            MemoryError

        <b/>See also Faust.toarray
        """
        return np.matrix(F.toarray())


    def __getitem__(F, indices):
        """
        Indexes or Slices a Faust.

        Returns a Faust representing a submatrix of F.toarray() or a scalar element if that Faust can be reduced to a single element.

        This function overloads a Python built-in.

        WARNING:
            - This function doesn't implement F[l1,l2] where l1 and l2 are
            integer list objects, rather use F[l1][:,l2].
            - It is not advised to use this function as an element accessor
        (e.g. F(0,0)) because such a use induces to convert the Faust to its
        dense matrix representation and that is a very expensive computation if used
        repetitively.
            - Subindexing a Faust which would create an empty Faust will raise
            an error.
            - fancy indexing must be done with a list not a numpy array.

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
            >>> i1 = randint(1, min(F.shape)-1)
            >>> i2 = randint(1, min(F.shape)-1)

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
            >>> F[0:i2:2,:] # takes every row of even index until the row i2 (excluded)
            >>> F[-1:-3:-1,:] # takes the two last rows in reverse order
            >>> F[i2:0:-2,:] # starts from row i2 and goes backward to take one in two rows until the first one (reversing order of F)
            >>> F[[1,18,2],:] # takes in this order the rows 1, 18 and 2
            >>> F[:, [1,18,2]] # takes in this order the columns 1, 18 and 2
            >>> F[[1,18,2], [1,2]] # takes the rows 1, 18 and 2 but keeps only columns 1 and 2 in these rows
        """
        #TODO: refactor (by index when indices == tuple(2), error message,
        #      out_indices checking on the end)
        empty_faust_except = Exception("Cannot create empty Faust.")
        idx_error_exception = IndexError("only integers, slices (`:`), ellipsis"
                                     " (`...`), and integer are valid indices")
        if(indices == Ellipsis): # F[...]
            out_indices = [slice(0,F.shape[0]), slice(0, F.shape[1])]
        elif(isinstance(indices,int)): # F[i] # a line
            out_indices = [slice(indices, indices+1), slice(0, F.shape[1])]
        elif(isinstance(indices,slice)):
            #F[i:j] a group of contiguous lines
            out_indices = [indices, slice(0, F.shape[1])]
        elif(isinstance(indices, list)):
            out_indices = [indices, list(range(0,F.shape[1]))]
            #TODO: check indices are all integers lying into F shape
        elif(isinstance(indices, tuple)):
            if(len(indices) == 2):
                out_indices = [0,0]
                if(isinstance(indices[0], int) and isinstance(indices[1],int)):
                    return F.todense()[indices[0],indices[1]]
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
                elif(isinstance(indices[0], list)):
                    if(len(indices[0]) == 0): raise empty_faust_except
                    out_indices[0] = indices[0]
                else:
                     raise idx_error_exception
                if(indices[1] == Ellipsis):
                    # all lines
                    out_indices[1] = slice(0, F.shape[1])
                elif(isinstance(indices[1], int)):
                    # line F[i]
                    out_indices[1] = slice(indices[1], indices[1]+1)
                elif(isinstance(indices[1], slice)):
                    out_indices[1] = indices[1]
                elif(isinstance(indices[1], list)):
                    if(isinstance(indices[0],list)): raise \
                    Exception("F[list1,list2] error: fancy indexing "
                              "on both dimensions is not implemented "
                              "rather use F[list1][:,list2].")
                    if(len(indices[1]) == 0): raise empty_faust_except
                    out_indices[1] = indices[1]
                else:
                     raise idx_error_exception
            else:
                raise IndexError('Too many indices.')
        else:
            raise idx_error_exception

        for i in range(0,2):
            if(isinstance(out_indices[i], slice)):
                if(out_indices[i].start == None and out_indices[i].stop == None):
                    #F[::] or F[::,any] or F[any, ::]
                    out_indices[i] = slice(0,F.shape[i],out_indices[i].step)
                elif(out_indices[i].start == None): # out_indices[i].stop != None
                    out_indices[i] = slice(0, out_indices[i].stop,
                                           out_indices[i].step)
                elif(out_indices[i].stop == None): # out_indices[i].start != None
                    out_indices[i] = slice(out_indices[i].start,
                                           F.shape[i],out_indices[i].step)
                if(out_indices[i].stop < 0):
                    out_indices[i] = slice(out_indices[i].start,
                                           F.shape[i]+out_indices[i].stop,
                                           out_indices[i].step)

                if(out_indices[i].step == None):
                    out_indices[i] = slice(out_indices[i].start,
                                                out_indices[i].stop,
                                                1)
                if(out_indices[i].start < 0):
                    out_indices[i] = \
                    slice(F.shape[i]+out_indices[i].start,out_indices[i].stop,
                         out_indices[i].step)

                if(out_indices[i].start >= F.shape[i] or out_indices[i].stop > F.shape[i]):
                    raise IndexError("index "+
                                     str(max(out_indices[i].start,out_indices[i].stop-1))+
                                     " is out of bounds for axis "+str(i)+" with size "+
                                     str(F.shape[i]))
                
                # transform slice with neg. step to a list for using fancy
                # indexing
                # likewise for step > 1
                if(out_indices[i].step < 0 or out_indices[i].step > 1):
                    out_indices[i] = \
                    list(range(out_indices[i].start,out_indices[i].stop,out_indices[i].step))
                    if(len(out_indices[i]) == 0): raise empty_faust_except
                elif(out_indices[i].start >= out_indices[i].stop):
                    raise empty_faust_except
                elif(out_indices[i].step == 0):
                    raise ValueError("slice step cannot be zero")

        if(isinstance(out_indices[0],list) or
                      isinstance(out_indices[1], list)):
            sub_F = Faust(core_obj=F.m_faust.fancy_idx(out_indices))
        else:
            sub_F = Faust(core_obj=F.m_faust.slice(out_indices))

        return sub_F

    def nnz_sum(F):
        """
        Gives the total number of non-zero elements in the factors of F.

        The function sums together the number of non-zero elements of
        each factor and returns the result. Note that for efficiency the sum is
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
        Computes the Relative Complexity Gain.

        RCG is the theoretical gain brought by the Faust representation relatively to its dense
        matrix equivalent. <br/>The higher is the RCG, the more computational
        savings will be made.
        This gain applies both for storage space and computation time.

        NOTE: F.rcg() == 1/F.density()

        Args:
            F: the Faust object.

        Returns:
            the RCG value (float).

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
        Computes the norm of F.

        Several types of norm are available: 1-norm, 2-norm, inf-norm and Frobenius norm.

        The norm of F is equal to the numpy.linalg.norm of F.toarray().

        WARNING: if [n,m] == size(F), the computation time can be expected to
        be of the order of min(n,m) times that of multipliying F by a vector.
        Nevertheless, the implementation ensures that memory usage remains
        controlled by avoiding to explicitly compute F.toarray().

        Args:
            F: the Faust object.
            ord: (optional) the norm order (1, 2, numpy.inf) or "fro" for
            Frobenius norm (by default the Frobenius norm is computed).

        Returns:
            the norm (float).
            <br/>If ord == 1, the norm is <code>norm(F.toarray(),1) ==
            max(sum(abs(F.toarray())))</code>.
            <br/>If ord == 2, the norm is the maximum singular value of F
            or approximately <code>norm(F.toarray(), 2) ==
            max(scipy.linalg.svd(F.toarray())[1])</code>.
            <br/> &nbsp;&nbsp;&nbsp; This is the default norm calculated when
            calling to norm(F).
            <br/>If ord = numpy.inf the norm is
            <code>norm(F.toarray(),numpy.inf) == max(sum(abs(F.T.toarray())))</code>
            <br/>If ord == 'fro', the norm is <code>norm(F.toarray(),
            'fro')</code>.

        Raises:
            ValueError.


        Examples:
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand([1, 2], [50, 100], .5)
            >>> F.norm()
            23.55588891399667
            >>> F.norm(2)
            5.929720822717308
            >>> F.norm('fro')
            23.55588891399667
            >>> F.norm(numpy.inf)
            18.509101197826254
        """
        if(ord not in [1, 2, "fro", np.inf]):
            raise ValueError("ord must have the value 1, 2, 'fro' or numpy.inf.")
        return F.m_faust.norm(ord)


    def normalize(F, ord='fro', axis=1):
        """
        Normalize F along the axis dimension using the ord-norm.

        Args:
            ord: the norm order to use (see Faust.norm).
            axis: if 1 the columns are normalized, if 0 the rows.

        Returns:
            the normalized Faust

        <b/> See also: Faust.norm
        """
        if(ord not in [1, 2, np.inf, "fro"]):
            raise ValueError("ord must have the value 1, 2, 'fro' or "
                             "numpy.inf.")
        if(axis not in [0,1]):
            raise ValueError("Invalid axis.")
        if(ord == float('Inf') or ord == np.inf):
            ord = -1
        elif(ord == "fro"):
            ord = -2
        if(axis == 0):
            tF = F.T
            NF = Faust(core_obj=tF.m_faust.normalize(ord))
            NF = NF.T
        else:
            NF = Faust(core_obj=F.m_faust.normalize(ord))
        return NF

    def get_num_factors(F):
        """
        Returns the number of factors of F.

        Returns:
            the number of factors.

        Examples:
        >>> from pyfaust import FaustFactory
        >>> F = FaustFactory.rand(2, 100, .5)
        >>> nf = F.get_num_factors()

        <b/> See also Faust.get_factor
        """
        return F.m_faust.get_nb_factors()

    def get_factor(F, indices):
        """
        Returns the i-th factor of F.

        Args:
            F: the Faust object.
            indices: the factor contiguous indices.

        Returns:
            if indices is a single index: a copy of the i-th factor.
            Otherwise: a new Faust composed of copies of the contigous factors of F
            pointed by indices.

            Each copy type is:
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
            >>> G = F.get_factor(range(3:5)) # a new Faust composed of the two last factors of F

        <b/> See also Faust.get_num_factors, Faust.transpose
        """
        if(hasattr(indices, '__iter__')):
            indices = list(indices)
        else:
            indices = list([indices])
        factors = []
        oi = None
        for i in indices:
            if(not isinstance(i, int)):
                raise not_int_e
            if(oi != None and i-oi != 1):
                raise Exception("Index must be contiguous.")
            factors += [F.m_faust.get_fact_opt(i)]
            oi = i
        if(len(factors) == 1):
            return factors[0]
        else:
            return pyfaust.Faust(factors)

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
        Saves the Faust F into a file.

        The file is saved in Matlab format version 5 (.mat extension).

        NOTE: storing F should typically use rcg(F) times less disk space than
        storing F.todense().

        Args:
            F: the Faust object.
            filepath: the path for saving the Faust (should end with .mat
            if Matlab format is used).
            format: (optional) The format to use for
            writing. By default, it's "Matlab" to save the Faust in a .mat
            file (currently only that format is available).

        Raises:
            ValueError.


        Examples:
            >>> from pyfaust import FaustFactory, Faust
            >>> F = FaustFactory.rand([2, 3], [10, 20],.5, field='complex')
            >>> F.save("F.mat")
            >>> G = Faust(filepath="F.mat")

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
            >>> F = FaustFactory.rand([2, 3], [10, 20],.5, field='complex')
            dtype('complex128')
            >>> F = FaustFactory.rand([2, 3], [10, 20],.5)
            dtype('float64')


        """
        if(F.m_faust.isReal()):
            return np.dtype(np.float64)
        else:
            return np.dtype(np.complex)

    def imshow(F):
        """
        Displays image of F's full matrix and its factors.

        <b> See also Faust.display. </b>

        Examples:
        >>> from pyfaust import FaustFactory
        >>> import matplotlib.pyplot as plt
        >>> F = FaustFactory.rand([2, 3], [10, 20],.5, field='complex')
        >>> F.imshow()
        >>> plt.show()


        """
        import matplotlib.pyplot as plt
        nf = F.get_num_factors()
        plt.subplot(1,nf+1,1)
        plt.imshow(abs(F.toarray()),extent=[0,100,0,1], aspect='auto')
        plt.xticks([]); plt.yticks([])
        for i in range(0,nf):
            plt.subplot(1,F.get_num_factors()+1,i+2)
            fac = F.get_factor(i)
            if(not isinstance(fac, np.ndarray)):
                fac = fac.toarray()
            plt.xticks([]); plt.yticks([])
            plt.imshow(abs(fac),extent=[0,100,0,1], aspect='auto')

    def pinv(F):
        """
            Computes the (Moore-Penrose) pseudo-inverse of a Faust full matrix.

        """
        from numpy.linalg.linalg import pinv
        return pinv(F.toarray())

    @staticmethod
    def isFaust(obj):
        """
        Returns True if obj is a Faust object, False otherwise.
        """
        return isinstance(obj, Faust)


import pyfaust.factparams

class FaustFactory:
    """
    This factory class provides methods for generating a Faust especially by factorization of a dense matrix.

    This class gives access to the main factorization algorithms of
    FAµST. Those algorithms can factorize a dense matrix to a sparse product
    (i.e. a Faust object).

    There are two algorithms for factorization.

    - The first one is Palm4MSA :
    which stands for Proximal Alternating Linearized Minimization for
    Multi-layer Sparse Approximation. Note that Palm4MSA is not
    intended to be used directly. You should rather rely on the second algorithm.

    - The second one is the Hierarchical Factorization algorithm:
    this is the central algorithm to factorize a dense matrix to a Faust.
    It makes iterative use of Palm4MSA to proceed with the factorization of a given
    dense matrix.

    A more secondary functionality of this class is the Faust generation.
    Several methods are available:

    - The pseudo-random generation of a Faust with FaustFactory.rand(),
    - the discrete Fourier transform with FaustFactory.dft(),
    - the Hadamard transform with FaustFactory.wht(),
    - and the identity Faust with FaustFactory.eye().

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
        >>> from pyfaust import FaustFactory
        >>> from pyfaust.factparams import ParamsPalm4MSA, ConstraintReal,\
        >>> ConstraintInt, ConstraintName, StoppingCriterion
        >>> import numpy as np
        >>> M = np.random.rand(500, 32)
        >>> cons = [ ConstraintInt('splin', 500, 32, 5),
        >>>          ConstraintReal('normcol', 32, 32, 1.0)]
        >>> stop_crit = StoppingCriterion(num_its=200)
        >>> param = ParamsPalm4MSA(cons, stop_crit)
        >>> F = FaustFactory.fact_palm4msa(M, param)
        >>> F.display()
        Faust size 500x32, density 0.22025, nnz_sum 3524, 2 factor(s):<br/>
        FACTOR 0 (real) SPARSE, size 500x32, density 0.15625, nnz 2500<br/>
        FACTOR 1 (real) SPARSE, size 32x32, density 1, nnz 1024<br/>
        """
        if(not isinstance(p, pyfaust.factparams.ParamsPalm4MSA)):
            raise TypeError("p must be a ParamsPalm4MSA object.")
        FaustFactory._check_fact_mat('FaustFactory.fact_palm4msa()', M)
        if(not p.is_mat_consistent(M)):
            raise ValueError("M's number of columns must be consistent with "
                             "the last residuum constraint defined in p. "
                             "Likewise its number of rows must be consistent "
                             "with the first factor constraint defined in p.")
        return Faust(core_obj=FaustCorePy.FaustFact.fact_palm4msa(M, p))

    @staticmethod
    def fact_hierarchical(M, p, ret_lambda=False, ret_params=False):
        """
        Factorizes the matrix M with Hierarchical Factorization using the parameters set in p.

        Args:
            M: the numpy matrix to factorize.
			p: is a set of factorization parameters. It might be a fully defined instance of parameters (pyfaust.factparams.ParamsHierarchicalFact) or a simplified expression which designates a pre-defined parametrization:
				- 'squaremat' to use pre-defined parameters typically used to factorize a Hadamard square matrix of order a power of two (see pyfaust.demo.hadamard).
				- ['rectmat', j, k, s] to use pre-defined parameters used for instance in factorization of the MEG matrix which is a rectangular matrix of size m*n such that m < n (see pyfaust.demo.bsl); j is the number of factors, k the sparsity of the main factor's columns, and s the sparsity of rows for all other factors except the residuum (that is the first factor here because the factorization is made toward the left -- is_side_fact_left == true, cf. pyfaust.factparams.ParamsHierarchicalFact).
				</br>The residuum has a sparsity of P*rho^(num_facts-1). <br/> By default, rho == .8 and P = 1.4. It's possible to set custom values with for example p == ( ['rectmat', j, k, s], {'rho':.4, 'P':.7 }). <br/>The sparsity is here the number of non-zero elements.
				@note - The fully defined parameters (ParamsHierarchicalFact
                instance) used/generated by the function are available in the
                return result (so one can consult what precisely mean the
            simplified parametrizations and possibly adjust the attributes to factorize again).
				@note - This function has its shorthand pyfaust.faust_fact(). For convenience you might use it like this:
				<code>
				from pyfaust import *
				F = faust_fact(M, p) % equiv. to FaustFactory.fact_hierarchical(M, p)
				</code>

            ret_lambda: set to True to ask the function to return the scale factor (False by default).
            ret_params: set to True to ask the function to return the
            ParamsHierarchicalFact instance used (False by default).

        Returns:
            F the Faust object result of the factorization.<br/>
            if ret_lambda == True (and ret_params == False), then the function
            returns the tuple (F,_lambda) (_lambda is the scale factor at the
            end of factorization).<br/>
            if ret_params == True (and ret_lambda == False), then the function
            returns the tuple (F, p) (p being the ParamsHierarchicalFact
            instance really used by the algorithm).<br/>
            if ret_lambda == True and ret_params == True, then the function
            returns the tuple (F, _lambda, p).

        Examples:
            <b> 1. Fully Defined Parameters for a Random Matrix Factorization </b>
            >>> from pyfaust import FaustFactory
            >>> from pyfaust.factparams import ParamsHierarchicalFact, ConstraintReal,\
            >>> ConstraintInt, ConstraintName, StoppingCriterion
            >>> import numpy as np
            >>> M = np.random.rand(500, 32)
            >>> fact_cons = [ ConstraintInt('splin', 500, 32, 5), ConstraintInt('sp', 32, 32, 96),
            >>>               ConstraintInt('sp', 32, 32, 96)]
            >>> res_cons = [ ConstraintReal('normcol', 32, 32, 1), ConstraintInt('sp', 32, 32, 666),
            >>>              ConstraintInt('sp', 32, 32, 333)]
            >>> stop_crit1 = StoppingCriterion(num_its=200)
            >>> stop_crit2 = StoppingCriterion(num_its=200)
            >>> param = ParamsHierarchicalFact(fact_cons, res_cons, stop_crit1, stop_crit2)
            >>> F = FaustFactory.fact_hierarchical(M, param)
            Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorisation
            1/3<br/>
            Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorisation
            2/3<br/>
            Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorisation
            3/3<br/>
            >>> F
            Faust size 500x32, density 0.189063, nnz_sum 3025, 4 factor(s):
                - FACTOR 0 (real) SPARSE, size 500x32, density 0.15625, nnz 2500
                - FACTOR 1 (real) SPARSE, size 32x32, density 0.09375, nnz 96
                - FACTOR 2 (real) SPARSE, size 32x32, density 0.09375, nnz 96
                - FACTOR 3 (real) SPARSE, size 32x32, density 0.325195, nnz 333

            <b>2. Simplified Parameters for Hadamard Factorization</b>

            >>> from pyfaust import FaustFactory as Facto
            >>> from numpy.linalg import norm
            >>> # generate a Hadamard Faust of size 32x32
            >>> FH = Facto.wht(5);
            >>> H = FH.toarray(); # the full matrix version
            >>> # factorize it
            >>> FH2 = Facto.fact_hierarchical(H, 'squaremat');
            >>> # test the relative error
            >>> (FH-FH2).norm('fro')/FH.norm('fro') # the result is 1.1015e-16, the factorization is accurate
            >>> FH
            Faust size 32x32, density 0.3125, nnz_sum 320, 5 factor(s):
            - FACTOR 0 (real) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 1 (real) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 2 (real) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 3 (real) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 4 (real) SPARSE, size 32x32, density 0.0625, nnz 64

            >>> FH2
            Faust size 32x32, density 0.3125, nnz_sum 320, 5 factor(s):
                - FACTOR 0 (real) SPARSE, size 32x32, density 0.0625, nnz 64
                - FACTOR 1 (real) SPARSE, size 32x32, density 0.0625, nnz 64
                - FACTOR 2 (real) SPARSE, size 32x32, density 0.0625, nnz 64
                - FACTOR 3 (real) SPARSE, size 32x32, density 0.0625, nnz 64
                - FACTOR 4 (real) SPARSE, size 32x32, density 0.0625, nnz 64

           <b> 3. Simplified Parameters for a Rectangular Matrix Factorization
           (the BSL demo MEG matrix)</b>
           >>> from pyfaust import *
           >>> from scipy.io import loadmat
           >>> from pyfaust.demo import get_data_dirpath
           >>> d = loadmat(get_data_dirpath()+'/matrix_MEG.mat')
           >>> MEG = d['matrix']
           >>> num_facts = 9
           >>> k = 10
           >>> s = 8
           >>> MEG16 = FaustFactory.fact_hierarchical(MEG, ['rectmat', num_facts, k, s])
           >>> MEG16
           Faust size 204x8193, density 0.0631655, nnz_sum 105573, 9 factor(s):
               - FACTOR 0 (real) SPARSE, size 204x204, density 0.293613, nnz
               12219
               - FACTOR 1 (real) SPARSE, size 204x204, density 0.0392157, nnz
               1632
               - FACTOR 2 (real) SPARSE, size 204x204, density 0.0392157, nnz
               1632
               - FACTOR 3 (real) SPARSE, size 204x204, density 0.0392157, nnz
               1632
               - FACTOR 4 (real) SPARSE, size 204x204, density 0.0392157, nnz
               1632
               - FACTOR 5 (real) SPARSE, size 204x204, density 0.0392157, nnz
               1632
               - FACTOR 6 (real) SPARSE, size 204x204, density 0.0392157, nnz
               1632
               - FACTOR 7 (real) SPARSE, size 204x204, density 0.0392157, nnz
               1632
               - FACTOR 8 (real) SPARSE, size 204x8193, density 0.0490196, nnz
               81930

           >>> # verify the constraint k == 10, on column 4
           >>> count_nonzero(MEG16.get_factor(8)[:,4].toarray())
           10
           >>> # now verify the s constraint is respected on MEG16 factor 1
           >>> count_nonzero(MEG16.get_factor(1).toarray())/MEG16.shape[0]
           8.0

        """
        from pyfaust.factparams import (ParamsHierarchicalFact,
        ParamsFactFactory)
        if(not isinstance(p, ParamsHierarchicalFact) and
           ParamsFactFactory.is_a_valid_simplification(p)):
            p = ParamsFactFactory.createParams(M, p)
        if(not isinstance(p, ParamsHierarchicalFact)):
            raise TypeError("p must be a ParamsHierarchicalFact object.")
        FaustFactory._check_fact_mat('FaustFactory.fact_hierarchical()', M)
        if(not isinstance(ret_lambda, bool)):
            raise TypeError("ret_lambda must be a bool.")
        if(not p.is_mat_consistent(M)):
            raise ValueError("M's number of columns must be consistent with "
                             "the last residuum constraint defined in p. "
                             "Likewise its number of rows must be consistent "
                             "with the first factor constraint defined in p.")
        core_obj,_lambda = FaustCorePy.FaustFact.fact_hierarchical(M, p)
        F = Faust(core_obj=core_obj)
        ret_list = [ F ]
        if(ret_lambda):
            ret_list += [ _lambda ]
        if(ret_params):
            ret_list += [ p ]
        if(ret_lambda or ret_params):
            return ret_list
        else:
            return F

    @staticmethod
    def wht(n):
        """
           Constructs a Faust implementing the Hadamard transform of dimension 2^n.

           The resulting Faust has n sparse factors of order 2^n, each one having 2 non-zero
           elements per row and per column.

           Args:
               n: the power of two exponent for a Hadamard matrix of order 2^n
               and a factorization in n factors.

           Returns:
               The Faust implementing the Hadamard transform of dimension 2^n.

          Examples:
              >>> from pyfaust import FaustFactory
              >>>  FaustFactory.wht(10)
              Faust size 1024x1024, density 0.0195312, nnz_sum 20480, 10
              factor(s):
                  - FACTOR 0 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
                  - FACTOR 1 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
                  - FACTOR 2 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
                  - FACTOR 3 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
                  - FACTOR 4 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
                  - FACTOR 5 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
                  - FACTOR 6 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
                  - FACTOR 7 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
                  - FACTOR 8 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
                  - FACTOR 9 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048

        """
        H = Faust(core_obj=FaustCorePy.FaustCore.hadamardFaust(n))
        return H

    @staticmethod
    def dft(n):
        """
            Constructs a Faust whose the full matrix is the Discrete Fourier Transform square matrix of order 2^n.

            The factorization algorithm used is Cooley-Tukey.

            The resulting Faust is complex and has (n+1) sparse factors whose the n first
            has 2 non-zero elements per row and per column. The last factor is
            a permutation matrix.

            Args:
            n: the power of two exponent for a DFT of order 2^n and a
            factorization in n+1 factors.

            Returns:
            The Faust implementing the DFT of dimension 2^n.

            Examples:
                >>> from pyfaust import FaustFactory
                >>> FaustFactory.dft(10)
                Faust size 1024x1024, density 0.0205078, nnz_sum 21504, 11 factor(s):
                - FACTOR 0 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
                - FACTOR 1 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
                - FACTOR 2 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
                - FACTOR 3 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
                - FACTOR 4 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
                - FACTOR 5 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
                - FACTOR 6 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
                - FACTOR 7 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
                - FACTOR 8 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
                - FACTOR 9 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
                - FACTOR 10 (complex) SPARSE, size 1024x1024, density 0.000976562, nnz 1024
        """
        F = Faust(core_obj=FaustCorePy.FaustCore.fourierFaust(n))
        return F

    @staticmethod
    def eye(m,n=None,t='real'):
        """
            Identity matrix as a Faust object.

            Args:
              m: number of rows,
              n (optional): number of columns, set to m if not specified.
              t (optional): 'complex' to return a complex Faust otherwise (by default)
              it's a real Faust.
        """
        from scipy.sparse import eye
        if(not n):
            n = m
        e = eye(m,n).tocsr()
        if(t == 'complex'):
            e = e.astype(np.complex)
        elif(t != 'real'):
            raise ValueError("t must be 'real' or 'complex'")
        return Faust(e)



    @staticmethod
    def rand(num_factors, dim_sizes, density=None, fac_type="mixed",
                  field='real', per_row=True):
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
                default value is such that each factor will have 5 non-zero
                elements per row, if per_row is True, or per column otherwise.
                It should be a floating point number between 0 and 1.
                fac_type: (optional) the storage representation of factors. Must be
                "sparse", "dense" or "mixed" if you want a mix of dense and
                sparse matrices in the generated Faust (choice's done according
                to an uniform distribution).
                            The default value is "mixed".
                field: (optional) a str to set the Faust field: 'real' or
                'complex'. By default, the value is 'real'.
                per_row: if True the factor matrix is constructed per row
                applying the density to each row. If False the construction is
                made with the density applied on each column.

        Returns:
            the random Faust.

        Examples:
            >>> from pyfaust import FaustFactory
            >>> F = FaustFactory.rand(2, 10, .5, field='complex')
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
        if(field == 'real'):
            is_real = True
        elif(field == 'complex'):
            is_real = False
        else:
            raise ValueError('field argument must be either "real" or "complex".')
        DENSE=0
        SPARSE=1
        MIXED=2
        REAL=3
        COMPLEX=4
        # set repr. type of factors
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
        if(not isinstance(per_row, bool)):
           raise ValueError("FaustFact.rand(): per_row argument must be a "
                            "bool.")
        if(not density):
            density = -1
        elif(not isinstance(density, np.float)):
            raise ValueError("FaustFactory.rand(): density must be a float")
        rF = Faust(core_obj=FaustCorePy.FaustCore.randFaust(
            fac_type_map[fac_type], field, min_num_factors, max_num_factors,
            min_dim_size, max_dim_size, density, per_row))
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


pyfaust.Faust.__div__ = pyfaust.Faust.__truediv__

def version():
    """Returns the FAµST package version.
    """
    return "@CPACK_PACKAGE_VERSION@"

def faust_fact(*args, **kwargs):
    """
    This function is a shorthand for FaustFactory.fact_hierarchical().

    <b/> See also FaustFactory.fact_hierarchical
    """
    return FaustFactory.fact_hierarchical(*args, **kwargs)

