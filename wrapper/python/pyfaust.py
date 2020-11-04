# -*- coding: @PYTHON_ENCODING@ -*-
# @PYFAUST_LICENSE_HEADER@
# @TORCH_LIBS_LOADING@

## @package pyfaust @brief @b The @b FAuST @b Python @b Wrapper

import numpy as np, scipy
from scipy.io import loadmat
from scipy.sparse import csr_matrix, csc_matrix, dia_matrix
import _FaustCorePy
import pyfaust
import pyfaust.factparams
import warnings

class Faust:
    """<b>FAuST Python wrapper main class</b> for using multi-layer sparse transforms.

    This class provides a numpy-like interface for operations
    with FAuST data structures, which correspond to matrices that can be
    written exactly as the product of sparse matrices.

    A FAuST data structure is designed to allow fast matrix-vector multiplications
    together with reduced memory storage compared to what would be obtained by
    manipulating directly the corresponding dense matrix.

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
    almost handled as a native numpy array.

    The main exception is that contrary to a numpy array a Faust is immutable.
    It means that you cannot affect elements of a Faust using
    the assignment operator `=' like you do with a numpy array (e.g. `M[i,j] =
    2').
    That limitation is the reason why the Python built-in `__setitem__()' is not
    implemented in this class.

    Other notable limitations are that one cannot:
        - compute the real and imaginary parts of a Faust,
        - perform elementwise operations between two Fausts (e.g. elementwise
        multiplication), the addition and subtraction are available though,
        - reshape a Faust.

    A last but not least caveat is that Faust doesn't support numpy universal
    functions (ufuncs) except if the contrary is specified in the API doc. for
    a particular function.

    Primarily for convenience and test purposes, a Faust can be converted into
    the corresponding full matrix using the function Faust.toarray or
    Faust.toarray.

    Warning: using Faust.toarray is discouraged except for test purposes, as it
    loses the main potential interests of the FAuST structure: compressed
    memory storage and faster matrix-vector multiplication compared to its
    equivalent full matrix representation.

    In this documentation, the expression 'full matrix' designates the array
    Faust.toarray() obtained by the multiplication of the Faust factors.

    List of functions that are memory costly: Faust.toarray(), Faust.toarray(),
    Faust.pinv().

    For more information about FAuST take a look at http://faust.inria.fr.
    """

    def  __init__(F, factors=None, filepath=None, **kwargs):
        """ Creates a Faust from a list of factors or alternatively from a file.

            Other easy ways to create a Faust is to call one of the
            following functions: pyfaust.rand(),
            pyfaust.dft() or pyfaust.wht().

        Args:
            factors: list of numpy/scipy array/matrices or a single array/matrix.<br/>
                     The factors must respect the dimensions needed for
                     the product to be defined <code>(for i in range(0,len(factors)-1): factors[i].shape[1] == factors[i+1].shape[0])</code>.<br/>
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
        use filepath you must explicitly set argument with the keyword.

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

            >>> F.save("F.mat")
            >>> # define a Faust from file
            >>> H = Faust(filepath="F.mat")

            >>> Faust(np.random.rand(10,10)) # creating a Faust with only one
                                             # factor

        <b/> See also Faust.save, pyfaust.rand

        """
        is_on_gpu = False
        gpu_dev = 0
        if("scale" in kwargs.keys()):
            # scale hidden argument
            scale = kwargs['scale']
            if(not isinstance(scale, (float, int, np.complex))):
                raise Exception("Scale must be a number.")
        else:
            scale = 1.0
        if("dev" in kwargs.keys()):
            dev = kwargs['dev']
            if dev.startswith('gpu'):
                is_on_gpu = True
        if("core_obj" in kwargs.keys()):
            core_obj = kwargs['core_obj']
            if(core_obj):
                F.m_faust = core_obj
        else:
            if(isinstance(factors, str) and not filepath):
                filepath = factors
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
                if(is_on_gpu):
                    F.m_faust = _FaustCorePy.FaustCoreGPU(factors, scale);
                else:
                    F.m_faust = _FaustCorePy.FaustCore(factors, scale);
            else:
                raise Exception("Cannot create an empty Faust.")

    @property
    def nbytes(F):
        """
        Gives the memory size of the Faust in bytes.
        """
        return F.m_faust.nbytes()

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
            >>> from pyfaust import rand
            >>> F = rand(2, 50, field='complex')
            >>> nrows, ncols = F.shape
            >>> nrows = F.shape[0]
            >>> ncols = F.shape[1]

        <b/> See also Faust.nnz_sum, Faust.size
        """
        return F.m_faust.shape()

    @property
    def ndim(F):
        """
            Number of Faust dimensions (always 2).
        """
        return 2

    @property
    def size(F):
        """
        Gives the size of the Faust F (that is F.shape[0]*F.shape[1]) .

        It's equivalent to numpy.prod(F.shape)).

        This function is intended to be used as a property (see the examples).

        Args:
            F: the Faust object.

        Returns:
            The number of elements in the Faust F.

        Examples:
            >>> from pyfaust import rand
            >>> F = rand(2, 50, field='complex')
            >>> F.size
            2500

        <b/> See also Faust.shape

        """
        return np.prod(F.shape)

    @property
    def device(self):
        """
        Returns the device on which the Faust is located (cpu or gpu).
        """
        if(isinstance(self.m_faust, (_FaustCorePy.FaustCoreGPU,
                      _FaustCorePy.FaustCore))):
            return self.m_faust.device()
        else:
            raise TypeError('This Faust object is invalid.')

    def transpose(F):
        """
        Returns the transpose of F.

        Args:
            F: the Faust object.

        Returns:
            a Faust object implementing the transpose of F.toarray(), such
            that:
            <code>F.transpose().toarray() == F.toarray().transpose()</code>

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
            a Faust object implementing the transpose of F.toarray(), such
            that:
            <code>F.T.toarray() == F.toarray().T</code>

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
            a Faust object Fc implementing the conjugate of F.toarray() such
            that:
            <code>Fc.toarray() == F.toarray().conj()</code>


        Examples:
            >>> from pyfaust import rand
            >>> F = rand(5, 50, field='complex')
            >>> Fc = F.conj()

        <b/> See also Faust.transpose, Faust.getH, Faust.H
        """
        F_conj = Faust(core_obj=F.m_faust.conj())
        return F_conj

    def conjugate(F):
        """
        Returns the complex conjugate of F.

        Args:
            F: the Faust object.

        Returns:
            a Faust object Fc implementing the conjugate of F.toarray() such
            that:
            <code>Fc.toarray() == F.toarray().conjugate()</code>


        Examples:
            >>> from pyfaust import rand
            >>> F = rand(5, 50, field='complex')
            >>> Fc = F.conjugate()

        <b/> See also Faust.transpose, Faust.getH, Faust.H
        """
        return F.conj()

    def getH(F):
        """
        Returns the conjugate transpose of F.

        Args:
            F: the Faust object.

        Returns:
            a Faust object H implementing the conjugate transpose of
            F.toarray() such that:
            <code>H.toarray() == F.toarray().getH()</code>

        Examples:
            >>> from pyfaust import rand
            >>> F = rand(5, 50, field='complex')
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
            F.toarray() such that:
            <code>H.toarray() == F.toarray().H</code>

        Examples:
            >>> from pyfaust import rand
            >>> F = rand(5, 50, field='complex')
            >>> H1 = F.H
            >>> H2 = F.transpose()
            >>> H2 = H2.conj()
            >>> (H1.toarray() == H2.toarray()).all()
            True

        <b/> See also Faust.transpose, Faust.conj, Faust.getH
        """
        return F.getH()

    def pruneout(F, thres=0, **kwargs):
        """
        Returns a Faust optimized by removing useless zero rows and columns as many times as needed.

        Args:
            F: the Faust to optimize.
            thres: the threshold of number of nonzeros under what the
            rows/columns are removed.


        Returns:
            The optimized Faust.

        """
        # hidden parameters (useless unless for debug)
        #            only_forward: True for applying only the forward passes of removal.
        #            npasses: the number of passes to run, by default it goes until the
        #            optimal Faust is obtained.
        npasses = 'auto'
        only_forward = False
        if('npasses' in kwargs):
            npasses = kwargs['npasses']
        if('only_forward' in kwargs):
            only_forward = kwargs['only_forward']
        if(npasses == 'auto'):
            npasses = -1
        elif(not isinstance(npasses, int)):
            raise TypeError('npasses must be a int'
                            ' or \'auto\'')
        if(not isinstance(only_forward, bool)):
            raise TypeError('only_forward '
                            'must be a bool.')
        if(not isinstance(thres, int)):
            raise TypeError('thres '
                            'must be a int.')
        #print("only_forward=", only_forward, "npasses=", npasses)
        F_prunedout = Faust(core_obj=F.m_faust.zpruneout(thres, npasses,
                                                         only_forward))
        return F_prunedout

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
            >>> from pyfaust import rand
            >>> F = rand(2, 50)
            >>> F.__repr__()
            >>> # the same function is called when typing F in a terminal:
            >>> F

        <b/>See also Faust.nnz_sum, Faust.rcg, Faust.shape, Faust.factors,
        <b/>Faust.numfactors, Faust.display

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
            >>> from pyfaust import rand
            >>> F = rand([1, 2], [50, 100], .5)
            >>> F.display()
            Faust size 98x82, density 0.686909, nnz_sum 5520, 2 factors:<br/>
            FACTOR 0 (real) SPARSE, size 98x78, density 0.395081, nnz 3020<br/>
            FACTOR 1 (real) SPARSE, size 78x82, density 0.390869, nnz 2500<br/>

            >>> F
            Faust size 98x82, density 0.686909, nnz_sum 5520, 2 factors:<br/>
            FACTOR 0 (real) SPARSE, size 98x78, density 0.395081, nnz 3020<br/>
            FACTOR 1 (real) SPARSE, size 78x82, density 0.390869, nnz 2500<br/>
            <!-- >>> -->

        <b/>See also Faust.nnz_sum, Faust.density, Faust.shape, Faust.factors,
        <b/>Faust.numfactors, Faust.__repr__

        """
        print(F.__repr__())
        #F.m_faust.display()

    def __add__(F, *args):
        """
        Sums F to one or a sequence of variables. Faust objects, arrays or scalars.


        NOTE: This method overloads the Python function/operator +.

        Args:
            args: the list of variables to sum all together with F.
            These can be Faust objects, numpy arrays (and
            scipy.sparse.csr_matrix) or scalars.
            Faust and arrays/matrices must be of compatible size.

        Returns:
            the sum as a Faust object.

        Example:
            >>> from pyfaust import *
            >>> F = rand(5,10)
            >>> G = F+F
            >>> H = F+2

        See also Faust.__sub__
        """
        dev = 'cpu'
        if isinstance(F.m_faust, _FaustCorePy.FaustCoreGPU):
            dev = 'gpu'
        for i in range(0,len(args)):
            G = args[i]
            if(isinstance(G,Faust)):
                if(F.shape != G.shape):
                    raise Exception('Dimensions must agree')
                C = F.concatenate(G, axis=1)
                # Id = np.eye(int(C.shape[1]/2))
                Fid = eye(int(C.shape[1]/2), dev=dev)
                #F = C*Faust(np.concatenate((Id,Id)),axis=0)
                F = C*Fid.concatenate(Fid,axis=0)
            elif(isinstance(G,np.ndarray) or
                 isinstance(G,scipy.sparse.csr_matrix) or
                 isinstance(G,scipy.sparse.csc_matrix)):
                F = F+Faust(G)
            elif(isinstance(G,(int, float, np.complex))):
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

    def __radd__(F,lhs_op):
        """
        Returns lhs_op+F.

        <b/>See also Faust.__add__
        """
        return F.__add__(lhs_op)

    def __sub__(F,*args):
        """
        Subtracts from F one or a sequence of variables. Faust objects, arrays or scalars.

        NOTE: This method overloads the Python function/operator -.

        Args:
            args: the list of variables to compute the difference with F. These can
            be Faust objects, arrays (and scipy.sparse.csr_matrix) or scalars.
            Faust and arrays/matrices must be of compatible size.

        Returns:
            the difference as a Faust object.

        Example:
            >>> from pyfaust import *
            >>> F = rand(5,10)
            >>> G = F-F
            >>> H = F-2

        <b/>See also Faust.__add__
        """

        nargs = []
        for arg in args:
            nargs += [ arg*-1 ]
        return F.__add__(*nargs)

    def __rsub__(F,lhs_op):
        """
        Returns lhs_op-F.

        <b/>See also Faust.__sub__
        """
        return F.__mul__(-1).__radd__(lhs_op)

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
        if(isinstance(s, (float, np.complex, int))):
            return F*(1./s)
        else:
            raise Exception("unsupported operand type(s) for /: a Faust can only be "
                  "divided by a scalar.")

    def __matmul__(F, A):
        """
        Multiplies F by A which is a dense numpy.matrix/numpy.ndarray or a Faust object.

        @warning The operator @ is not supported in Python 2, you can use *.

        This method overloads a Python function/operator (@).

        <b>The primary goal</b> of this function is to implement “fast” multiplication by a
        Faust, which is the operation performed when A is a dense matrix.<br/>
        In the best case, F @ A is F.rcg() times faster than equivalent
        F.toarray() @ A.

        <b>Other use case</b> is available for this function:
        - If A is a Faust, no actual multiplication is performed, instead a
        new Faust is built to implement the multiplication.<br/>
        This Faust verifies that:<br/>
            <code>
            (F @ A).toarray() == F.toarray()@A.toarray()
            </code>
            <br/>N.B.: you could have an elementwise non-significant absolute
            difference between the two members.

        Args:
            F: the Faust object.
            A: a Faust object, a sparse matrix (scipy.sparse.csr_matrix or dia_matrix) or a 2D full matrix (numpy.ndarray, numpy.matrix).
            <br/> In the latter case, A must be Fortran contiguous (i.e. Column-major order;
                `order' argument passed to np.ndararray() must be equal to str
                'F').
            <br/> Note that performing a Faust-csr_matrix product is often
            slower than converting first the csr_matrix to a dense
            representation (toarray()) and then proceed to the
            Faust-dense_matrix multiplication. In some cases though, it stays
            quicker: moreover when the Faust is composed of a small number of
            factors.

        Returns: The result of the multiplication
            - as a numpy.ndarray if A is a ndarray,<br/>
            - as a Faust object if A is a Faust.
            - When either F or A is complex, G=F @ A is also complex.

        Raises: ValueError
            The multiplying operand A is a scalar:
            <code>
            >>> from pyfaust import rand
            >>> F = rand(5,50)
            >>> F@2
            ValueError Scalar operands are not allowed, use '*' instead
            </code>

        Examples:
            >>> from pyfaust import rand
            >>> import numpy as np
            >>> F = rand(5, [50, 100])
            >>> A = np.random.rand(F.shape[1], 50)
            >>> B = F@A # == F*A or F.dot(A)
            >>> # is equivalent to B = F.__matmul__(A)
            >>> G = rand(5, F.shape[1])
            >>> H = F@G
            >>> # H is a Faust because F and G are

        <b/>See also Faust.__init__, Faust.rcg, Faust.__mul__, Faust.__matmul__, Faust.dot
        """
        if(isinstance(A, Faust)):
            if(F.shape[1] != A.shape[0]): raise ValueError("The dimensions of "
                                                          "the two Fausts must "
                                                          "agree.")
            return Faust(core_obj=F.m_faust.multiply(A.m_faust))
        elif(isinstance(A, (float, int, np.complex))):
            raise ValueError("Scalar operands are not allowed, use '*'"
                             " instead")
        elif(isinstance(A, np.ndarray) and A.dtype == np.complex and F.dtype !=
            np.complex):
            j = np.complex(0,1)
            return F.m_faust.multiply(A.real).astype(np.complex) + j*F.m_faust.multiply(A.imag)
        elif(isinstance(A, scipy.sparse.csr_matrix)):
            if(A.dtype == np.complex and F.dtype != np.complex):
                j = np.complex(0,1)
                return (F.__matmul__(A.real))*np.complex(1,0) + \
                        (F.__matmul__(A.imag))*j
            else:
                return F.m_faust.multiply_csr_mat(A.astype(F.dtype))
        elif(isinstance(A, (dia_matrix, csc_matrix))):
            return F.__matmul__(A.tocsr())
        else:
            return F.m_faust.multiply(A)

    def dot(F, A):
        """
        Performs equivalent operation of numpy.dot() between the Faust F and A.

        More specifically:
            - Scalar multiplication if A is a scalar but F*A is preferred.
            - Matrix multiplication if A is a Faust or numpy.ndarray/numpy.matrix but F @ A is preferred.

        <b/>See also Faust.__init__, Faust.rcg, Faust.__mul__, Faust.__matmul__, Faust.dot
        """
        if(isinstance(A, (float, int, np.complex))):
            return F*A
        return F.__matmul__(A)

    def matvec(F, x):
        """
        This function implements the scipy.sparse.linalg.LinearOperator.matvec function such that scipy.sparse.linalg.aslinearoperator function works on a Faust object.
        """
        return F.dot(x)

    def __mul__(F, A):
        """
        Multiplies the Faust F by A.

        This method overloads a Python function/operator (*).

        More specifically:
        - It's a scalar multiplication if A is a scalar number.
        - If A is a Faust or a numpy.ndarray it performs a
        Faust.__matmul__() (@ operator).

        Args:
            F: the Faust object.
            A: is a scalar number, a Faust object or a numpy.ndarray
            or a sparse matrix (scipy.sparse.csr_matrix or dia_matrix).
            <br/> In the latter case, A must be Fortran contiguous (i.e. Column-major order;
                `order' argument passed to np.ndararray() must be equal to str
                'F').

        Returns: The result of the multiplication
            - as a numpy.ndarray if A is a ndarray,<br/>
            - as a Faust object if A is a Faust or a scalar number.
            - When either F or A is complex, G=F*A is also complex.

        Raises: TypeError
            If the type of A is not supported:
            <code>
            >>> import numpy
            >>> F*'a'
            TypeError ufunc 'multiply' did not contain a loop with signature matching types dtype('<U32') dtype('<U32') dtype('<U32')
            </code>

        Examples:
            >>> from pyfaust import rand
            >>> import numpy as np
            >>> F = rand(5, [50, 100])
            >>> A = np.random.rand(F.shape[1], 10)
            >>> B = F*A
            >>> # is equivalent to B = F.__mul__(A)
            >>> G = rand(5, F.shape[1])
            >>> H = F*G
            >>> # H is a Faust because F and G are
            >>> ((F*G).toarray() == (F@G).toarray()).all() #@ is not supported in python2
            True
            >>> F_times_two = F*2

        <b/>See also Faust.__init__, Faust.rcg, Faust.__mul__, Faust.__matmul__, Faust.dot
        """
        if(isinstance(A, (float, int, np.complex))):
            return Faust(core_obj=F.m_faust.multiply_scal(A))
        else: # A is a Faust, a numpy.ndarray (eg. numpy.matrix) or anything
            try:
                return F.__matmul__(A)
            except:
                raise TypeError("invalid type operand for Faust.__mul__.")

    def __rmul__(F, lhs_op):
        """ lhs_op*F

        <b/>See also Faust.__mul__
        """
        if(isinstance(lhs_op, (float, int, np.complex))):
            return Faust(core_obj=F.m_faust.multiply_scal(lhs_op))
        else:
            try:
                return F.__rmatmul__(lhs_op)
            except:
                raise TypeError("invalid type operand for Faust.__rmul__.")

    __array_ufunc__ = None # mandatory to override rmatmul
                           # it means Faust doesn't support ufuncs

    def __rmatmul__(F, lhs_op):
        """
        Returns lhs_op.__matmul__(F).

        <b>See also Faust.__matmul__</b>

        Examples:
            >>> from pyfaust import rand
            >>> import numpy as np
            >>> F = rand(5, [50, 100])
            >>> A = np.random.rand(50, F.shape[0])
            >>> B = A@F # == A*F or pyfaust.dot(A,F)


        """
        if(isinstance(lhs_op, (np.ndarray, csr_matrix, dia_matrix))):
            if(F.dtype == np.complex or lhs_op.dtype == np.complex):
                return (F.T.conj().__matmul__(lhs_op.T.conj())).T.conj()
            else: # real Faust and real lhs_op
                return (F.T.__matmul__(lhs_op.T)).T
        else:
            raise TypeError("invalid type operand for Faust.__matmul__.")

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
                >>> F = rand(2,51);
                >>> G = rand(2,25);
                >>> F.concatenate(G, 0)
                </code>
                ValueError
                Axis must be 0 or 1.
                <code>
                >>> F = rand(2,51);
                >>> G = rand(2,25);
                >>> F.concatenate(G, 5)
                </code>
                ValueError
                You can't concatenate a Faust with something that is not a
                Faust, a numpy array or a scipy sparse matrix.
                <code>
                >>> F = rand(2,51);
                >>> F.concatenate(['a', 'b', 'c'], 5)
                </code>

            Examples:
                >>> import pyfaust import rand
                >>> F = rand(5, 50)
                >>> G = rand(6, 50)
                >>> F.concatenate(G) # equivalent to F.concatenate(G, 0)
                Faust size 100x50, density 0.5634, nnz_sum 2817, 7 factor(s)<br/>
                FACTOR 0 (real) SPARSE, size 100x100, density 0.0473, nnz 47<br/>
                FACTOR 1 (real) SPARSE, size 100x100, density 0.0469, nnz 46<br/>
                FACTOR 2 (real) SPARSE, size 100x100, density 0.0504, nnz 50<br/>
                FACTOR 3 (real) SPARSE, size 100x100, density 0.0502, nnz 50<br/>
                FACTOR 4 (real) SPARSE, size 100x100, density 0.0482, nnz 48<br/>
                FACTOR 5 (real) SPARSE, size 100x100, density 0.0287, nnz 28<br/>
                FACTOR 6 (real) SPARSE, size 100x50, density 0.02, nnz 10<br/>
                >>> F.concatenate(G, axis=1)
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
                Faust size 384x50, density 0.677083, nnz_sum 13000, 12 factor(s):<br/>
                FACTOR 0 (real) SPARSE, size 384x400, density 0.0224609, nnz 3450<br/>
                FACTOR 1 (real) SPARSE, size 400x400, density 0.01125, nnz 1800<br/>
                FACTOR 2 (real) SPARSE, size 400x400, density 0.01125, nnz 1800<br/>
                FACTOR 3 (real) SPARSE, size 400x400, density 0.01125, nnz 1800<br/>
                FACTOR 4 (real) SPARSE, size 400x400, density 0.01125, nnz 1800<br/>
                FACTOR 5 (real) SPARSE, size 400x350, density 0.00714286, nnz 1000<br/>
                FACTOR 6 (real) SPARSE, size 350x300, density 0.00333333, nnz 350<br/>
                FACTOR 7 (real) SPARSE, size 300x250, density 0.004, nnz 300<br/>
                FACTOR 8 (real) SPARSE, size 250x200, density 0.005, nnz 250<br/>
                FACTOR 9 (real) SPARSE, size 200x150, density 0.00666667, nnz 200<br/>
                FACTOR 10 (real) SPARSE, size 150x100, density 0.01, nnz 150<br/>
                FACTOR 11 (real) SPARSE, size 100x50, density 0.02, nnz 100<br/>
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


        C=F
        for G in args:
           if(isinstance(G, (np.ndarray, csr_matrix, csc_matrix))):
               G = Faust([G])
           if(not isinstance(G, Faust)): raise ValueError("You can't concatenate a "
                                                           "Faust with something "
                                                           "that is not a Faust, "
                                                           "a numpy array or scipy "
                                                           "sparse matrix.")
           if(axis == 0 and F.shape[1] != G.shape[1] or axis == 1 and F.shape[0]
              != G.shape[0]): raise ValueError("The dimensions of "
                                               "the two Fausts must "
                                               "agree.")

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
            >>> from pyfaust import rand
            >>> F = rand(5, 10**6, .00001, 'sparse')
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
        return F.m_faust.get_product()


    def todense(F):
        """
        Returns a numpy matrix for the full matrix implemented by F.

        WARNING: Using this function is discouraged except for test purposes,
        as it loses the main potential interests of the Faust structure:
        compressed memory storage and faster matrix-vector multiplication
        compared to its equivalent full matrix representation.

        WARNING: this function is deprecated in favor to toarray function and
        will be removed in a future version. The cause behind is that this
        function returns a numpy.matrix which is deprecated too.

        Returns:
            A numpy matrix M such that M*x == F*x for any vector x.

        Raises:
            MemoryError

        WARNING: running the example below is likely to raise a memory
        error or freeze your computer for a certain amount of time.

        Examples:
            >>> from pyfaust import rand
            >>> F = rand(5, 10**6, .00001, 'sparse')
            >>> F
            Faust size 1000000x1000000, density 5e-05, nnz_sum 49999995, 5 factors:<br/>
            FACTOR 0 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999<br/>
            FACTOR 1 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999<br/>
            FACTOR 2 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999<br/>
            FACTOR 3 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999<br/>
            FACTOR 4 (real) SPARSE, size 1000000x1000000, density 1e-05, nnz 9999999<br/>
            >>> # an attempt to convert F to a dense matrix is most likely to raise a memory error
            >>> # the sparse format is the only way to handle such a large Faust
            >>> F.toarray()
            (...)
            MemoryError

        <b/>See also Faust.toarray
        """
        warnings.warn("Faust.todense() is deprecated and will be deleted in the "
             "future.")
        return np.matrix(F.toarray())


    def __getitem__(F, indices):
        """
        Indexes or slices a Faust.

        Returns a Faust representing a submatrix of F.toarray() or a scalar element if that Faust can be reduced to a single element.

        This function overloads a Python built-in.

        WARNING:
            - This function doesn't implement F[l1,l2] where l1 and l2 are
            integer lists, rather use F[l1][:,l2].
            - It is not advised to use this function as an element accessor
        (e.g. F[0,0]) because such a use induces to convert the Faust to its
        dense matrix representation and that is a very expensive computation if used
        repetitively.
            - Subindexing a Faust which would create an empty Faust will raise
            an error.
            - 'Fancy indexing' must be done with a list not a numpy array.

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
            >>> from pyfaust import rand
            >>> import numpy as np
            >>> from random import randint

            >>> F = rand(5, [50, 100])
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
            >>> F[[1,18,2]][:,[1,2]] # takes the rows 1, 18 and 2 but keeps only columns 1 and 2 in these rows
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
            if(len(indices) == 1):
                return F.__getitem__(indices[0])
            if(len(indices) == 2):
                out_indices = [0,0]
                if(isinstance(indices[0], int) and isinstance(indices[1],int)):
                    return F.toarray()[indices[0],indices[1]]
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
        compared to storing the corresponding dense matrix F.toarray(), as well
        as potentially faster matrix-vector multiplication when applying F * x
        instead of F.toarray()*x.

        NOTE: A density above one is possible but prevents any saving.

        Args:
            F: the Faust object.

        Returns:
            the density value (float).

        Examples:
        >>> from pyfaust import rand
        >>> F = rand(5, 50, .5)
        >>> dens = F.density()

        <b/> See also Faust.nnz_sum, Faust.rcg, Faust.size, Faust.toarray
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

    def norm(F, ord='fro', **kwargs): #**kwargs):
        """
        Computes the norm of F.

        Several types of norm are available: 1-norm, 2-norm, inf-norm and Frobenius norm.

        The norm of F is equal to the numpy.linalg.norm of F.toarray().

        WARNING: the computation time can be expected to be of order
        n*min(F.shape) with n the time for multipliying F by a vector.
        Nevertheless, the implementation ensures that memory usage remains
        controlled by avoiding to explicitly compute F.toarray() (at least for
        2-norm).

        Args:
            F: the Faust object.
            ord: (optional) the norm order (1, 2, numpy.inf) or "fro" for
            Frobenius norm (by default the Frobenius norm is computed).
            threshold: (optional) power iteration algorithm threshold (default
            to .001). Used only for norm(2).
            max_num_its: (optional) maximum number of iterations for
            power iteration algorithm. Used only for norm(2).

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
            >>> from pyfaust import rand
            >>> import numpy as np
            >>> F = rand([1, 2], [50, 100], .5)
            >>> F.norm()
            23.55588891399667
            >>> F.norm(2)
            5.929720822717308
            >>> F.norm('fro')
            23.55588891399667
            >>> F.norm(np.inf)
            18.509101197826254
        """
        if(ord not in [1, 2, "fro", np.inf]):
            raise ValueError("ord must have the value 1, 2, 'fro' or numpy.inf.")
        return F.m_faust.norm(ord, **kwargs)


    def normalize(F, ord='fro', axis=1):
        """
        Normalizes F along the axis dimension using the ord-norm.

        The function is able to normalize the columns (default axis=1):
        <br/><code>
            NF = F.normalize(ord) is such that for all i in
        range(0,F.shape[1]) NF.toarray()[:,i] ==
        F.toarray()[:,i]/norm(F.toarray(), ord)
        </code>

        Likewise for normalization of the rows (axis=0):
        <br/><code>
            NF = F.normalize(ord, axis=0) is such that for all i in
            range(0,F.shape[0]) NF.toarray()[i,:] ==
            F.toarray()[i,:]/norm(F.toarray(), ord)
        </code>

        The variable ord designates one of the Faust.norm() compatible norms.

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
            if(ord == -1):
                ord = 1
            elif(ord == 1):
                ord = -1
            NF = Faust(core_obj=tF.m_faust.normalize(ord))
            NF = NF.T
        else:
            NF = Faust(core_obj=F.m_faust.normalize(ord))
        return NF

    def numfactors(F):
        """
        Returns the number of factors of F.

        Returns:
            the number of factors.

        Examples:
        >>> from pyfaust import rand
        >>> F = rand(2, 100, .5)
        >>> nf = F.numfactors()
        >>> nf == len(F)
        True

        <b/> See also Faust.factors, Faust.__len__
        """
        return F.m_faust.get_nb_factors()

    def __len__(F):
        """
        Returns the number of factors of F.

        Returns:
            the number of factors.

        Examples:
        >>> from pyfaust import rand
        >>> F = rand(2, 100, .5)
        >>> nf = F.numfactors()
        >>> nf == len(F)
        True

        <b/> See also Faust.factors, Faust.numfactors
        """
        return F.numfactors()

    def factors(F, indices):
        """
        Returns the i-th factor of F.

        Args:
            F: the Faust object.
            indices: the factor contiguous indices.

        Returns:
            if indices is a single index: a copy of the i-th factor.
            Otherwise: a new Faust composed of copies of the contiguous factors of F
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
            >>> from pyfaust import rand
            >>> F = rand(5, [50, 100], .5)
            >>> f0 = F.factors(0)
            >>> G = F.factors(range(3,5)) # a new Faust composed of the two last factors of F

        <b/> See also Faust.numfactors, Faust.transpose
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

    def right(F, i):
        """
        Returns the right hand side factors of F from index i to F.numfactors()-1.

        Returns:
            a Faust if the size factor set to be returned is greater than 1, a
            numpy array otherwise.

        Examples:
            >>> from pyfaust import rand
            >>> F = rand(5, [5, 10])
            >>> RF = F.right(2)
            >>> print(F)
            Faust size 10x7, density 2.85714, nnz_sum 200, 5 factor(s):<br/>
                 FACTOR 0 (real) SPARSE, size 10x10, density 0.5, nnz 50<br/>
                 FACTOR 1 (real) DENSE,  size 10x10, density 0.5, nnz 50<br/>
                 FACTOR 2 (real) SPARSE, size 10x5, density 1, nnz 50<br/>
                 FACTOR 3 (real) DENSE,  size 5x5, density 1, nnz 25<br/>
                 FACTOR 4 (real) SPARSE, size 5x7, density 0.714286, nnz 25<br/>
            >>> print(RF)
            Faust size 10x7, density 1.42857, nnz_sum 100, 3 factor(s):<br/>
                 FACTOR 0 (real) SPARSE, size 10x5, density 1, nnz 50<br/>
                 FACTOR 1 (real) DENSE,  size 5x5, density 1, nnz 25<br/>
                 FACTOR 2 (real) SPARSE, size 5x7, density 0.714286, nnz 25<br/>


        See also:
            Faust.factors, Faust.left, Faust.numfactors

        """
        i = F._check_factor_idx(i)
        rF = Faust(core_obj=F.m_faust.right(i))
        if(len(rF) == 1):
            return rF.factors(0)
        return rF

    def left(F, i):
        """
        Returns the left hand side factors of F from index 0 to i included.

        Returns:
            a Faust if the size of factor set to be returned is greater than 1, a
            numpy array otherwise.

        Examples:
            >>> from pyfaust import rand
            >>> F = rand(5, [5, 10])
            >>> LF = F.left(3)
            >>> print(F)
            Faust size 6x10, density 3.25, nnz_sum 195, 5 factor(s):<br/>
                FACTOR 0 (real) DENSE,  size 6x8, density 0.625, nnz 30<br/>
                FACTOR 1 (real) DENSE,  size 8x9, density 0.555556, nnz 40<br/>
                FACTOR 2 (real) SPARSE, size 9x10, density 0.5, nnz 45<br/>
                FACTOR 3 (real) DENSE,  size 10x6, density 0.833333, nnz 50<br/>
                FACTOR 4 (real) DENSE,  size 6x10, density 0.5, nnz 30<br/>
            >>> print(LF)
            Faust size 6x6, density 4.58333, nnz_sum 165, 4 factor(s):<br/>
                FACTOR 0 (real) DENSE,  size 6x8, density 0.625, nnz 30<br/>
                FACTOR 1 (real) DENSE,  size 8x9, density 0.555556, nnz 40<br/>
                FACTOR 2 (real) SPARSE, size 9x10, density 0.5, nnz 45<br/>
                FACTOR 3 (real) DENSE,  size 10x6, density 0.833333, nnz 50<br/>


        See also:
            Faust.factors, Faust.right
        """
        i = F._check_factor_idx(i)
        lF = Faust(core_obj=F.m_faust.left(i))
        if(len(lF) == 1):
            return lF.factors(0)
        return lF

    def _check_factor_idx(F, i):
        if(not isinstance(i, (float, int))):
            raise TypeError('i must be an integer.')
        i = int(np.floor(i))
        if(i < 0 or i >= F.numfactors()):
            raise ValueError('i is out of range.')
        return i

    def get_factor_nonopt(F, i):
        """
        DEPRECATED: use Faust.factors
        Returns the i-th factor of F.

        Args:
            F: the Faust object.
            i: the factor index.

        Returns:
            a copy of the i-th factor as a dense matrix (of type numpy.ndarray).

        Raises:
            ValueError.


        Examples:
            >>> from pyfaust import rand
            >>> F = rand(5, [50, 100], .5)
            >>> f0 = F.factors(0)

        <b/> See also Faust.numfactors
        """
        fact = F.m_faust.get_fact(i)
        return fact

    def save(F, filepath, format="Matlab"):
        """
        Saves the Faust F into a file.

        The file is saved in Matlab format version 5 (.mat extension).

        NOTE: storing F should typically use rcg(F) times less disk space than
        storing F.toarray().

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
            >>> from pyfaust import rand, Faust
            >>> F = rand([2, 3], [10, 20],.5, field='complex')
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
            >>> from pyfaust import rand
            >>> F = rand([2, 3], [10, 20],.5, field='complex')
            >>> F.dtype
            dtype('complex128')
            >>> F = rand([2, 3], [10, 20],.5)
            >>> F.dtype
            dtype('float64')


        """
        if(F.m_faust.isReal()):
            return np.dtype(np.float64)
        else:
            return np.dtype(np.complex)

    def imshow(F, name='F'):
        """
        Displays image of F's full matrix and its factors.

        Args:
            F: the Faust object.
            name: (optional str) the displayed name on the plotted figure.


        Examples:
        >>> from pyfaust import rand
        >>> import matplotlib.pyplot as plt
        >>> F = rand([2, 3], [10, 20],.5, field='complex')
        >>> F.imshow()
        >>> plt.show()


        <b/> See also Faust.display.
        """
        import matplotlib.pyplot as plt
        if(not isinstance(name, str)): raise TypeError('name must be a str.')
        nf = F.numfactors()
        max_cols = 5
        ncols = min(nf,max_cols)
        nrows = int(nf/ncols)+1
        plt.subplot(nrows,ncols,nrows*ncols)
        plt.title(name+'.toarray()')
        if(F.dtype == np.complex):
            plt.imshow(abs(F.toarray()), aspect='equal')
        else:
            plt.imshow(F.toarray(), aspect='equal')
        #plt.xticks([]); plt.yticks([])
        for i in range(0,nf):
            plt.subplot(nrows,ncols,(i%ncols)+int(i/ncols)*ncols+1)
            plt.title(str(i))
            fac = F.factors(i)
            if(not isinstance(fac, np.ndarray)):
                fac = fac.toarray()
            plt.xticks([]); plt.yticks([])
            plt.suptitle('Factors of the Faust '+ name)
            if(fac.dtype == np.complex):
                plt.imshow(abs(fac),aspect='equal')
            else:
                plt.imshow(fac, aspect='equal')

    def pinv(F):
        """
        Computes the (Moore-Penrose) pseudo-inverse of a Faust full matrix.

        Warning: this function makes a call to Faust.toarray().

        Returns:
            The dense pseudo-inverse matrix.
        """
        from numpy.linalg.linalg import pinv
        return pinv(F.toarray())

    @staticmethod
    def isFaust(obj):
        """
        Returns True if obj is a Faust object, False otherwise.

        Examples:
            >>> from pyfaust import *
            >>> Faust.isFaust(2)
            False
            >>> Faust.isFaust(rand(5,10))
            True

        """
        return isinstance(obj, Faust)

    def optimize_memory(F):
        """
        Optimizes a Faust by changing the storage format of each factor in order to optimize the memory size.

        Returns:
            The optimized Faust.

        <b/> See also Faust.optimize
        """
        F_opt = Faust(core_obj=F.m_faust.optimize_storage(False))
        return F_opt

    def optimize(F, transp=False):
        """
        Returns a Faust optimized with Faust.pruneout, Faust.optimize_memory and Faust.optimize_time.

        Args:
            transp: True in order to optimize the Faust according to its transpose.

        Returns:
            The optimized Faust.

        Note: this function is still experimental, you might use manually
        Faust.optimize_time, Faust.optimize_memory or Faust.pruneout to be
        more specific about the optimization to proceed.

        <b/> See also Faust.optimize_time, Faust.optimize_memory, Faust.pruneout
        """
        F_opt = Faust(core_obj=F.m_faust.optimize(transp))
        return F_opt

    def optimize_time(F, transp=False, inplace=False, nsamples=1):
        """
        Returns a Faust configured with the quickest Faust-matrix multiplication mode (benchmark ran on the fly).

        NOTE: this function launches a small benchmark on the fly. Basically, the methods
        available differ by the order used to compute the matrix chain
        multiplication or by the use (or unuse) of libraries to performs the
        calculation.
        The evaluated methods in the benchmark are listed in
        pyfaust.FaustMulMode but note that FaustMulMode.CPP_PROD_PAR_REDUC and
        FaustMulMode.OMP_PROD_PAR_REDUC are excluded from the benchmark because
        it doesn't worth it in any case when Eigen multithread is enabled
        (which is the case in any package of pyfaust delivered).
        Although depending on the package you installed and the capability of your
        hardware the methods based on Torch or GPU Cuda library can be used.

        Args:
            inplace: to optimize the current Faust directly instead of returning a new
            Faust with the optimization enabled. If True, F is returned
            otherwise a new Faust object is returned.
            transp: True in order to optimize the Faust according to its transpose.
            nsamples: the number of Faust-Dense matrix products
            calculated in order to measure time taken by each method (it could matter
            to discriminate methods when the performances are similar). By default,
            only one product is computed to evaluate the method.

        Returns:
            The optimized Faust.

        <b/> See also Faust.optimize

        """
        if(inplace):
            F.m_faust.optimize_time(transp, inplace, nsamples)
            return F
        else:
            F_opt = Faust(core_obj=F.m_faust.optimize_time(transp, inplace,
                                                          nsamples))
            return F_opt

    def clone(F, dev='cpu'):
        """
        Clones the Faust (in a new memory space).

        Args:
            dev: 'cpu' to clone on CPU RAM, 'gpu:id' to clone on the GPU device
            #id (e.g. gpu:0).

        Returns:
            The Faust clone.
        """
        if isinstance(F.m_faust, _FaustCorePy.FaustCoreGPU):
             clone_F = Faust(core_obj=F.m_faust.clone(dev))
             return clone_F
        else:
            return Faust([F.factors(i) for i in range(F.numfactors())], dev=dev)

pyfaust.Faust.__div__ = pyfaust.Faust.__truediv__

def version():
    """Returns the FAuST package version.
    """
    return "@CPACK_PACKAGE_VERSION@"

def faust_fact(*args, **kwargs):
    """
    This function is a shorthand for pyfaust.fact.hierarchical.

    <b/> See also pyfaust.fact.hierarchical
    """
    import pyfaust.fact
    return pyfaust.fact.hierarchical(*args, **kwargs)

def license():
    """ Prints the FAuST license.
    """
    print("""@PYFAUST_LICENSE_HEADER@""")


def norm(F, ord='fro', **kwargs):
    """
        Returns Faust.norm(F, ord) or numpy.linalg.norm(F, ord) depending of F type.

    <b/>See also Faust.norm
    """
    if(Faust.isFaust(F)):
        return F.norm(ord, **kwargs)
    else: # if F is not a Faust, try to rely on numpy (not breaking possible
          # past import)
        axis = None
        keepdims = False
        if('axis' in kwargs.keys()): axis = kwargs['axis']
        if('keepdims' in kwargs.keys()): keepdims = kwargs['keepdims']
        return np.linalg.norm(F, ord, axis=axis,
                              keepdims=keepdims)

def dot(A, B):
    """
        Returns Faust.dot(A,B) if A or B is a Faust object, returns numpy.dot(A,B) ortherwise.

    <b/>See also Faust.norm
    """
    if(Faust.isFaust(A)):
        return A.dot(B)
    elif(Faust.isFaust(B)):
        return B.T.dot(A.T).T
    else: # if F is not a Faust, try to rely on numpy (not breaking possible
          # past import)
        return np.dot(A,B)

def pinv(F):
    """
    A package function alias for the member function Faust.pinv().

    Args:
        F: is a Faust object.
    """
    if(Faust.isFaust(F)):
        return F.pinv()
    else:
        return np.linalg.linalg.pinv(F)

def concatenate(_tuple, axis=0):
    """
    A package function alias for the member function Faust.concatenate.


    <b>See also</b> numpy.concatenate

    Example:
        >>> from pyfaust import *
        >>> F1 = rand(5,50)
        >>> F2 = rand(5,50)
        >>> concatenate((F1, F2), axis=0)

    """
    if(Faust.isFaust(_tuple[0])):
        return _tuple[0].concatenate(*_tuple[1:], axis=axis)
    else:
        return np.concatenate(_tuple, axis=axis)

def hstack(_tuple):
    """
    Concatenates horizontally Faust and/or numpy.ndarray objects using Faust.concatenate().

    <b>See also</b> numpy.hstack

    Example:
        >>> from pyfaust import *
        >>> F1 = rand(5,50)
        >>> F2 = rand(5,50)
        >>> hstack((F1, F2))


    """
    return pyfaust.concatenate(_tuple, axis=1)

def vstack(_tuple):
    """
    Concatenates vertically Faust and/or numpy.ndarray arrays using Faust.concatenate().

    <b>See also</b> numpy.vstack

    Example:
        >>> from pyfaust import *
        >>> F1 = rand(5,50)
        >>> F2 = rand(5,50)
        >>> vstack((F1, F2))

    """
    return pyfaust.concatenate(_tuple, axis=0)

def isFaust(obj):
    """
        Returns True if obj is a Faust object, False otherwise.

        Examples:
            >>> from pyfaust import isFaust
            >>> isFaust(2)
            False
            >>> isFaust(rand(5,10))
            True
    """
    return Faust.isFaust(obj)

def wht(n, normed=True, dev="cpu"):
    """
       Constructs a Faust implementing the Hadamard transform of dimension n.

       The resulting Faust has n sparse factors of order n, each one having
       2 nonzero elements per row and per column.

       Args:
           n: the power of two exponent for a Hadamard matrix of order n
           and a factorization in log2(n) factors.
           normed: default to True to normalize the Hadamard Faust as if you called
           Faust.normalize() and False otherwise.
       Returns:
           The Faust implementing the Hadamard transform of dimension n.

      Examples:
          >>> from pyfaust import wht
          >>> wht(1024)
          Faust size 1024x1024, density 0.0195312, nnz_sum 20480, 10 factor(s):
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
          >>> wht(1024, normed=True) # is equiv. to next call
          >>> wht(1024, normed=False).normalize() # which is less optimized though

    """
    log2n = np.floor(np.log2(n))
    if(n > 2**log2n): raise ValueError("n must be a power of 2.")
    if(not isinstance(normed, bool)):
        raise TypeError("normed must be True of False.")
    if dev == "cpu":
        H = Faust(core_obj=_FaustCorePy.FaustCore.hadamardFaust(log2n, normed))
    elif dev.startswith("gpu"):
        H = Faust(core_obj=_FaustCorePy.FaustCoreGPU.hadamardFaust(log2n, normed))
    return H

def dft(n, normed=True):
    """
        Constructs a Faust F such that F.toarray() is the Discrete Fourier Transform square matrix of order n.

        The factorization algorithm used is Cooley-Tukey.

        The resulting Faust is complex and has (log2(n)+1) sparse factors
        whose the log2(n) first has 2 non-zero elements per row and per column.
        The last factor is a permutation matrix.

        Args:
            n: the power of two exponent for a DFT of order n and a
            factorization in log2(n)+1 factors.
            normed: default to True to normalize the DFT Faust as if you called
            Faust.normalize() and False otherwise.

        Returns:
            The Faust implementing the DFT of dimension n.

        Examples:
            >>> from pyfaust import dft
            >>> dft(1024)
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
            >>> dft(1024, normed=True) # is equiv. to next call
            >>> dft(1024, normed=False).normalize() # which is less optimized though
    """
    log2n = np.floor(np.log2(n))
    if(n > 2**log2n): raise ValueError("n must be a power of 2.")
    if(not isinstance(normed, bool)):
        raise TypeError("normed must be True of False.")
    F = Faust(core_obj=_FaustCorePy.FaustCore.fourierFaust(log2n, normed))
    return F

def eye(m,n=None,t='real', dev="cpu"):
    """
        Identity matrix as a Faust object.

        Args:
          m: number of rows,
          n (optional): number of columns, set to m if not specified.
          t (optional): 'complex' to return a complex Faust otherwise (by default)
          it's a real Faust.

        Examples:
            >>> from pyfaust import eye
            >>> eye(5)
            Faust size 5x5, density 0.2, nnz_sum 5, 1 factor(s):<br/>
            FACTOR 0 (real) SPARSE, size 5x5, density 0.2, nnz 5<br/>
            >>> eye(5,4)
            Faust size 5x4, density 0.2, nnz_sum 4, 1 factor(s):<br/>
            FACTOR 0 (real) SPARSE, size 5x4, density 0.2, nnz 4<br/>
            >>> eye(5,t='complex')
            Faust size 5x4, density 0.2, nnz_sum 4, 1 factor(s):<br/>
            FACTOR 0 (complex) SPARSE, size 5x4, density 0.2, nnz 4<br/>
    """
    if(t not in ['complex', 'real']):
        raise ValueError("t must be 'real' or 'complex'")
    if(n == None): n = m
    if dev == "cpu":
        rF = Faust(core_obj=_FaustCorePy.FaustCore.eyeFaust(m, n, t))
    elif dev.startswith("gpu"):
        rF = Faust(core_obj=_FaustCorePy.FaustCoreGPU.eyeFaust(m, n, t))
    return rF
#    from scipy.sparse import eye
#    if(not n):
#        n = m
#    e = eye(m,n).tocsr()
#    if(t == 'complex'):
#        e = e.astype(np.complex)
#    elif(t != 'real'):
#        raise ValueError("t must be 'real' or 'complex'")
#    return Faust(e)

def rand(num_factors, dim_sizes, density=None, fac_type="mixed",
              field='real', per_row=True, dev="cpu"):
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
        >>> from pyfaust import rand
        >>> F = rand(2, 10, .5, field='complex')
        >>> G = rand([2, 5], [10, 20], .5, fac_type="dense")
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
        raise ValueError('rand(): argument fac_type must be a'
                         ' str equal to "sparse",'
                         ' "dense" or "mixed".')

    fac_type_map = {"sparse" : SPARSE, "dense" : DENSE, "mixed" : MIXED}
    # set field of matrices/factors
    if(not isinstance(is_real, bool)):
        raise ValueError('rand(): argument is_real must be a'
                         'boolean.')
    if(is_real):
        field = REAL
    else:
        field = COMPLEX
    if((isinstance(num_factors,(list, tuple))) and
    len(num_factors) == 2):
        min_num_factors = num_factors[0]
        max_num_factors = num_factors[1]
    elif(isinstance(num_factors, int)):
        min_num_factors = max_num_factors = num_factors
    else:
        raise ValueError("rand(): num_factors argument must be an "
                         "integer or a list/tuple of 2 integers.")
    if(isinstance(dim_sizes, (list, tuple)) and len(dim_sizes) == 2):
        min_dim_size = dim_sizes[0]
        max_dim_size = dim_sizes[1]
    elif(isinstance(dim_sizes, (int, np.long))):
        min_dim_size = max_dim_size = dim_sizes
    else:
        raise ValueError("rand(): dim_sizes argument must be an "
                         "integer or a list/tuple of 2 integers.")
    if(not isinstance(per_row, bool)):
       raise ValueError("FaustFact.rand(): per_row argument must be a "
                        "bool.")
    if(not density):
        density = -1
    elif(not isinstance(density, np.float)):
        raise ValueError("rand(): density must be a float")
    if dev == "cpu":
        rF = Faust(core_obj=_FaustCorePy.FaustCore.randFaust(
            fac_type_map[fac_type], field, min_num_factors, max_num_factors,
            min_dim_size, max_dim_size, density, per_row))
    elif dev.startswith("gpu"):
        rF = Faust(core_obj=_FaustCorePy.FaustCoreGPU.randFaust(
            fac_type_map[fac_type], field, min_num_factors, max_num_factors,
            min_dim_size, max_dim_size, density, per_row))

    return rF

def enable_gpu_mod(libpaths=None, backend='cuda', silent=False, fatal=False):
    """
    This function loads explicitely the gpu_mod library in memory.

    Normally it's not required to load the library manually, but it could be
    useful to set a non-default path or to diagnose a loading issue.

    Args:
        libpaths: the absolute or relative paths where to search the dynamic
        library (gm) to load. By default, it's none to auto-find the library
        (if possible).
        backend: the GPU backend to use, only cuda is available for now.
        silent: if True nothing or almost will be displayed on loading (e.g.
        silent errors), otherwise all messages are visible.
    """
    _FaustCorePy.FaustCore.enable_gpu_mod(libpaths, backend, silent, fatal)

# experimental block start
import torch

class FaustTorch:

    def __init__(self, F, device=None):
        """
        Creates a FaustTorch from a Faust or a list of numpy arrays/sparse matrices.

        Args:
            F: a Faust or a list of factors to create the FaustTorch from.
            device: "cpu" or "cuda", by default cpu. The device on which the
            Tensors are instantiated.
        """
        self.factors = []
        err = ValueError('F must be a Faust or a list of numpy.ndarrays or'
                         ' scipy.sparse.csr_matrix')
        self.device = device
        if(isinstance(F, Faust)):
            F = [F.factors(i) for i in range(0, F.numfactors())]
        #TODO: tests all factors are ndim == 1,2 long
        if(isinstance(F, (list, tuple))):
            for factor in F:
                if(isinstance(factor, np.ndarray)):
                    self.factors += [torch.from_numpy(factor)]
                elif(isinstance(factor, csr_matrix)):
                    coo_mat = factor.tocoo()
                    row = coo_mat.row
                    col = coo_mat.col
                    I = torch.from_numpy(np.stack((row, col), axis=0))
                    V = torch.from_numpy(coo_mat.data)
                    self.factors += [
                        torch.sparse_coo_tensor(I,V,dtype=torch.float64)] #.coalesce() ]
                else:
                    raise err
                if(self.device != None):
                    self.factors[-1] = self.factors[-1].to(device)
        else:
            raise err

    def __mul__(self, op=None, optimize_chain=False, allow_none_op=False):
        """
        Multiplies the chain of tensors into a single tensor (matrix product).

        """
        if(isinstance(op, np.ndarray)):
            res = torch.from_numpy(op)
            if(self.device != None):
                res = res.to(self.device)
            factors = self.factors
            if(optimize_chain):
                return self.mul_opt(res)
        elif(allow_none_op and op == None):
            factors = self.factors[:]
            if(factors[-1].is_sparse):
                factors.append(torch.from_numpy(np.eye(self.factors[-1].size()[1])))
                if(self.device != None):
                    factors[-1] = factors[-1].to(self.device)
            res = factors[-1].clone()
            factors = factors[:-1]
            if(optimize_chain):
                tmp = self.factors
                self.factors = factors
                res = self.mul_opt(res)
                self.factors = tmp
                return res
        else:
            raise TypeError('op must be a np.ndarray')

        # torch matmul
        #res = self.factors[0]

        #for f in self.factors[1:]:
        for f in reversed(factors[:]):
            if(f.is_sparse):
                res = torch.sparse.mm(f, res)
            else:
                res = torch.matmul(f, res)
        return res

    def totensor(self, optimize_chain=False):
        """
         See Faust.toarray()
        """
        return self.__mul__(allow_none_op=True, optimize_chain=optimize_chain)

    def mul_opt(self, op):
        """
        Computes the product self*op optimizing the order of matrix chain products.
        """
        def cost(a,b):
#            if(a.is_sparse):
#                a_cost = a._nnz() #len(a.values())
#            else:
            a_cost = a.size()[0]*a.size()[1]
            if(b.is_sparse):
                b_cost = b._nnz()
            else:
                b_cost = b.size()[1]
            return a_cost*b_cost


        factors = self.factors.copy() + [ op ]
        costs = [cost(factors[i], factors[i+1]) for i in range(len(factors)-1)
                 ]
        while len(factors) > 1:
            i = np.argmin(costs)
            f1 = factors[i]
            f2 = factors[i+1]
            if(f2.is_sparse):
                 if(f1.is_sparse):
                    res = torch.sparse.mm(f1, f2.to_dense())
                 else:
                    res = torch.t(torch.mm(torch.t(f2), torch.t(f1)))
            else:
                if(f1.is_sparse):
                    res = torch.sparse.mm(f1, f2)
                else:
                    res = torch.matmul(f1, f2)
            factors = factors[0:i]+[res]+factors[i+2:]
            if(len(factors) > i+1):
                costs = costs[0:i] + [ cost(res, factors[i+1]) ] + costs[i+2:]
            else:
                costs = costs[:-1]
        return factors[0]

    def __repr__(self):
        _str = ''
        for i, f in enumerate(self.factors):
            _size = f.size()
            _str += "- FACTOR "+str(i)
            if(f.is_sparse):
                _str += " SPARSE "
            else:
                _str += " DENSE "
            _str += " size "
            _str += str(_size[0])+'x'+str(_size[1]) + '\n'
        return _str
# experimental block end


class FaustMulMode:
    """
    <b/> Enumeration class of all matrix chain multiplication methods available to multiply a Faust to a matrix, to a vector or to compute Faust.toarray().

    These methods are used by Faust.optimize_time().

    NOTE: it's not advisable to use these methods directly. The user should use
    Faust.optimize_time() but the access is left open for experimental purpose.

    Examples:
        >>> from pyfaust import rand as frand
        >>> from numpy.random import rand
        >>> F = frand(5, [100, 1024])
        >>> F.m_faust.set_FM_mul_mode(FaustMulMode.ORDER_ALL_BEST_MIXED) # method used to compute Faust-matrix product or Faust.toarray()
        >>> F.m_faust.set_Fv_mul_mode(FaustMulMode.DEFAULT) # method used to compute Faust-vector product
        >>> F*rand(F.shape[1],1) # Faust-vector mul. using the DEFAULT method
        >>> F*rand(F.shape[1], 512) # Faust-matrix mul. using method ORDER_ALL_BEST_MIXED
        >>> F.toarray() # using the same method
    """
    ## \brief The default method. Multiplying from the right to the left.
    DEFAULT=0
    ## \brief This method computes the product by its ends.
    ##
    ## For each iteration/multiplication it chooses to multiply the most right
    ## or the most left pair of matrices (in order to decrease the computation cost).
    ## The computational cost depends on the matrix dimensions.
    ORDER_ALL_ENDS=1
    ## \brief This method computes the product starting by the pair of matrices whose the computation cost is the smallest.
    ##
    ## After this first multiplication the rest of the factors are multiplied
    ## on the left side (from the right to the left) and then the right factors are
    ## multiplied (from the left to the right).
    ## The computational cost depends on the matrix dimensions.
    ORDER_1ST_BEST=2
    ## \brief This method computes a chain of matrices by ordering the product according to the minimal computation cost order.
    ##
    ## The computational cost depends on the matrix dimensions.
    ORDER_ALL_BEST_CONVDENSE=3
    ## \brief This method follows the same principle as ORDER_ALL_BEST_CONVDENSE method but is capable to multiply dense matrices as well as sparse matrices.
    ##
    ## The computational cost depends on the matrix dimensions and the number
    ## of nonzeros (when a matrix is in sparse format).
    ORDER_ALL_BEST_MIXED=4
    ## \brief This method computes the product performing a parallel reduction of the product.
    ##
    ## It uses as many threads as C++ STL advises (std::thread::hardware_concurrency() -- https://en.cppreference.com/w/cpp/thread/thread/hardware_concurrency).
    ##
    ## Reference: https://en.wikipedia.org/wiki/Reduce_%28parallel_pattern%29
    CPP_PROD_PAR_REDUC=5
    ## \brief This method is equivalent to CPP_PROD_PAR_REDUC but is implemented using OpenMP.
    OMP_PROD_PAR_REDUC=6
    ## \brief This method computes the product of the matrix chain from the left to the right using Torch C++ library (CPU backend).
    ##
    ## This method is only available for the specific packages pyfaust_torch.
    TORCH_CPU=7
    ## \brief This method computes the product following the minimal cost order using Torch C++ library (CPU backend).
    ##
    ## The method is the same as the one used for ORDER_ALL_BEST_MIXED but is implemented with Torch library.
    ##
    ## This method is only available for the specific packages pyfaust_torch.
    TORCH_CPU_BEST_ORDER=8
    ## \brief The same as TORCH_CPU except that torch::chain_matmul is used to
    ## compute in one call the intermediary product of dense contiguous
    ## factors, then the result is multiplied by sparse factors if any remains.
    ##
    ## This method is only available for the specific packages pyfaust_torch.
    TORCH_CPU_DENSE_ROW_TORCH=9
    ## \brief Use the GPU module to compute the product on a NVIDIA GPU.
    ##
    ## Multiplying from the left to the right or in the way around in order to minimize the cost.
    GPU_MOD=10
