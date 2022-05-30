# -*- coding: @PYTHON_ENCODING@ -*-
# @PYFAUST_LICENSE_HEADER@
# @BLASOMP_LOADING@
# @OMP_LIB_LOADING@
# @TORCH_LIBS_LOADING@

## @package pyfaust @brief @b The @b FAuST @b Python @b Wrapper

import numpy as np, scipy
import scipy
from scipy.io import loadmat
from scipy.sparse import (csr_matrix, csc_matrix, dia_matrix, bsr_matrix,
                          coo_matrix, diags, eye as seye, kron, vstack as
                          svstack, hstack as shstack)
import _FaustCorePy
import pyfaust
import pyfaust.factparams
import warnings
import decimal
import numpy.lib.mixins
from os import environ

HANDLED_FUNCTIONS = {}

class Faust(numpy.lib.mixins.NDArrayOperatorsMixin):
    """<b>FAuST Python wrapper main class</b> for using multi-layer sparse transforms.

    This class provides a numpy-like interface for operations
    with FAuST data structures, which correspond to matrices that can be
    written exactly as the product of sparse matrices.

    The Faust class is designed to allow fast matrix-vector multiplications
    together with reduced memory storage compared to what would be obtained by
    manipulating directly the corresponding dense matrix.

    A particular example is the matrix associated to the discrete Fourier
    transform, which can be represented exactly as a Faust,
    leading to a fast and compact implementation (see pyfaust.dft).

    Although sparse matrices are more interesting for optimization it's not
    forbidden to define a Faust as a product of dense matrices or a mix of dense
    and sparse matrices.

    The matrices composing the Faust product, also called the factors, are
    defined on complex or real fields. Hence a Faust can be a complex Faust or a
    real Faust.

    Several Python built-ins have been overloaded to ensure that a Faust is
    almost handled as a native numpy array.

    The main exception is that contrary to a numpy array a Faust is immutable.
    It means that you cannot modify elements of a Faust using
    the assignment operator `=' like you do with a numpy array (e.g. `M[i,j] =
    2').
    That limitation is the reason why the Python built-in `__setitem__()' is not
    implemented in this class.
    Note however that you can optionally contravene the immuability in
    certain functions (e.g. with the `inplace` argument of the
    functions Faust.swap_rows, Faust.swap_cols, Faust.optimize_time).

    Other noticeable limitations are that one cannot:
        - perform elementwise operations between two Fausts (e.g. elementwise
        multiplication), the addition and subtraction are available though,
        - reshape a Faust.

    A last but not least caveat is that Faust doesn't support the numpy universal
    functions (ufuncs) except if the contrary is specified in the API doc. for
    a particular function.

    Mainly for convenience and test purposes, a Faust can be converted into
    the corresponding full matrix using the function Faust.toarray.

    Warning: using Faust.toarray is discouraged except for test purposes, as it
    loses the main potential interests of the FAuST structure: compressed
    memory storage and faster matrix-vector multiplication compared to its
    equivalent full matrix representation.

    In this documentation, the expression 'full matrix' designates the array
    Faust.toarray() obtained by the multiplication of the Faust factors.

    List of functions that are memory costly:
    - element indexing (F[i,j] / __getitem__, but note that slicing is memory efficient through memory views),
    - function Faust.toarray(),
    - function Faust.pinv().

    For more information about FAuST take a look at http://faust.inria.fr.
    """

    def  __init__(F, factors=None, filepath=None, **kwargs):
        """ Creates a Faust from a list of factors or alternatively from a file.
        Other easy ways to create a Faust is to call one of the following functions: pyfaust.rand, pyfaust.dft or pyfaust.wht.

        Args:
            factors: list of numpy/scipy array/matrices or a single array/matrix.<br/>
                     The factors must respect the dimensions needed for
                     the product to be defined <code>(for i in range(0,len(factors)-1): factors[i].shape[1] == factors[i+1].shape[0])</code>.<br/>
                     The factors can be sparse or dense matrices
                     (either scipy.sparse.csr_matrix/bsr_matrix or
                     numpy.ndarray with ndim == 2).
                     scipy.sparse.csc_matrix/coo_matrix are supported but
                     converted to csr_matrix on the fly.<br/>
                     The Faust will be in the same dtype as the factors
                     (only 'float32', 'float64'/'double' and 'complex128'/'complex' dtype are supported).
                     Passing only an array or a sparse matrix to the
                     constructor is equivalent to
                     passing a list of a single factor.
            filepath: the file from which a Faust is created.<br/>
                      The format is Matlab version 5 (.mat extension).<br/>
                      The file must have been saved before with Faust.save().
            dev: 'cpu' (by default) or 'gpu' to instantiate the Faust resp. on
            CPU or on NVIDIA GPU memory.
            kwargs: internal purpose arguments.

        NOTE: filepath and factors arguments are mutually exclusive. Either
        you specify one of them explicitly with the keyword or the first (positional) argument type
        will determine if it's a filepath (a str) or a factor list. If both of
        the keyword arguments are set then filepath will be prioritary.

        Examples:
            Example 1 -- Creating a Faust made of CSR matrices and numpy arrays:
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
                >>> H = Faust(filepath="F.mat") # F = Faust("F.mat") does the same

                >>> Faust(np.random.rand(10,10)) # creating a Faust with only one
                                                 # factor

            Example 2 -- Creating a Faust containing one BSR matrix:
                >>> from scipy.sparse import bsr_matrix
                >>> from pyfaust import Faust
                >>> from numpy import allclose
                >>> from numpy.random import rand
                >>> nzblocks = rand(3, 2, 3) # 3 blocks of size 2x3
                >>> # create a scipy BSR matrix
                >>> B = bsr_matrix((nzblocks, [0, 1, 2], [0, 1, 2, 3, 3, 3]), shape=(10,9))
                >>> # create the single factor Faust
                >>> F = Faust(B)
                >>> F
                Faust size 10x9, density 0.2, nnz_sum 18, 1 factor(s):<br/>
                 FACTOR 0 (double) BSR, size 10x9, density 0.2, nnz 18
                >>> allclose(F.toarray(), B.toarray())
                True
                >>> # of course it's ok to create a Faust with a BSR and another type of factors
                >>> Faust([B, rand(9, 18)])
                Faust size 10x18, density 1, nnz_sum 180, 2 factor(s): <br/>
                  FACTOR 0 (double) BSR, size 10x9, density 0.2, nnz 18<br/>
                  FACTOR 1 (double) DENSE, size 9x18, density 1, nnz 162
                >>>


        <b>See also</b> Faust.save, pyfaust.rand, pyfaust.dft, pyfaust.wht
        <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html">csr_matrix, </a>
        <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.bsr_matrix.html">bsr_matrix</a>
        <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html">csc_matrix, </a>
        <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html">coo_matrix </a>

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
            check_dev(dev)
        if("core_obj" in kwargs.keys()):
            core_obj = kwargs['core_obj']
            if(core_obj):
                F.m_faust = core_obj
        else:
            if(isinstance(factors, str) and not filepath):
                filepath = factors
            if(filepath and isinstance(filepath, str)):
                G = Faust.load(filepath)
                F.m_faust = G.m_faust
                F._is_real = G._is_real
                return
            if isinstance(factors,
                           (np.ndarray, csc_matrix, csr_matrix, coo_matrix, bsr_matrix)) and factors.ndim == 2:
                factors = [ factors ]
            if(not isinstance(factors, list)):
                raise Exception("factors must be a non-empty list of/or a numpy.ndarray, "
                                "scipy.sparse.csr_matrix/csc_matrix/bsr_matrix.")
            F._is_real = True
            # verify if all factors has the same dtype (mandatory)
            dtype = None
            if len(factors) == 0:
                raise ValueError("Empty list of matrices.")
            for f in factors:
                F._is_real = f.dtype != np.complex
                if dtype is None:
                    dtype = f.dtype
                elif dtype != f.dtype:
                    raise TypeError('All Faust factors must have the same dtype.')
            if dtype not in ['double', 'float32', 'complex128']:
                raise TypeError('Unmanaged factor dtype:'+str(dtype)+' (must be float32, double or complex128')
            if(factors is not None and len(factors) > 0):
                if(is_on_gpu):
                    if F._is_real:
                        if dtype == 'double':
                            F.m_faust = _FaustCorePy.FaustCoreGenDblGPU(factors, scale)
                        else: # if dtype == 'float32':
                            F.m_faust = _FaustCorePy.FaustCoreGenFltGPU(factors, scale)
                    else:
                        F.m_faust = _FaustCorePy.FaustCoreGenCplxDblGPU(factors, scale)
                elif F._is_real:
                    if dtype == 'double':
                        F.m_faust = _FaustCorePy.FaustCoreGenDblCPU(factors, scale)
                    else: # if dtype == 'float32':
                        F.m_faust = _FaustCorePy.FaustCoreGenFltCPU(factors, scale)
                else:
                    F.m_faust = _FaustCorePy.FaustCoreGenCplxDblCPU(factors, scale)
            else:
                raise Exception("Cannot create an empty Faust.")

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method == '__call__':
            if str(ufunc) == "<ufunc 'matmul'>" and len(inputs) >= 2 and \
               isFaust(inputs[1]):
                return inputs[1].__rmatmul__(inputs[0])
            elif str(ufunc) == "<ufunc 'multiply'>" and len(inputs) >= 2 and \
               isFaust(inputs[1]):
                return inputs[1].__rmul__(inputs[0])
            elif str(ufunc) == "<ufunc 'add'>" and len(inputs) >= 2 and \
                    isFaust(inputs[1]):
                return inputs[1].__radd__(inputs[0])
            elif str(ufunc) == "<ufunc 'subtract'>" and len(inputs) >= 2 and \
                    isFaust(inputs[1]):
                return inputs[1].__rsub__(inputs[0])
            N = None
            fausts = []
        elif method == 'reduce':
#            # not necessary numpy calls Faust.sum
#            if ufunc == "<ufunc 'add'>":
#                if len(inputs) == 1 and pyfaust.isFaust(inputs[0]):
#                    #return inputs[0].sum(*inputs[1:], **kwargs)
#                else:
            return NotImplemented

    def __array__(self, *args, **kwargs):
        return self

    def __array_function__(self, func, types, args, kwargs):
        print("__array_function__")
        if func not in HANDLED_FUNCTIONS:
            return NotImplemented
        # Note: this allows subclasses that don't override
        # __array_function__ to handle MyArray objects
        if not all(issubclass(t, Faust) for t in types):
            return NotImplemented
        return HANDLED_FUNCTIONS[func](*args, **kwargs)

    @property
    def nbytes(F):
        """
        Gives the memory size of the Faust in bytes.

        Example:
            >>> import pyfaust as pf
            >>> F = pf.rand(1024, 1024)
            >>> F
            Faust size 1024x1024, density 0.0244141, nnz_sum 25600, 5 factor(s):
            - FACTOR 0 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 1 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 2 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 3 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 4 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120

            >>> F.nbytes
            327700
            >>> F = pf.rand(1024, 1024, fac_type='dense')
            >>> F
            Faust size 1024x1024, density 0.0244141, nnz_sum 25600, 5 factor(s):
            - FACTOR 0 (double) DENSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 1 (double) DENSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 2 (double) DENSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 3 (double) DENSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 4 (double) DENSE, size 1024x1024, density 0.00488281, nnz 5120

            >>> F.nbytes
            41943040

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
            >>> F = rand(50, 50, field='complex')
            >>> nrows, ncols = F.shape
            >>> nrows = F.shape[0]
            >>> ncols = F.shape[1]

        <b>See also</b> Faust.nnz_sum, Faust.size
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
            >>> F = rand(50, 50, field='complex')
            >>> F.size
            2500

        <b>See also</b> Faust.shape

        """
        return np.prod(F.shape)

    @property
    def device(self):
        """
        Returns the device on which the Faust is located ('cpu' or 'gpu').

        Example:
            >>> import pyfaust as pf
            >>> cpuF = pf.rand(5, 5, dev='cpu')
            >>> cpuF.device
            'cpu'
            >>> gpuF = pf.rand(5, 5, dev='gpu')
            >>> gpuF.device
            'gpu'
        """
        return self.m_faust.device()

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

        <b>See also</b> Faust.conj, Faust.getH, Faust.H, Faust.T
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

        <b>See also</b> Faust.conj, Faust.getH, Faust.H, Faust.T
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
            >>> F = rand(50, 50, field='complex')
            >>> Fc = F.conj()

        <b>See also</b> Faust.transpose, Faust.getH, Faust.H
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
            >>> F = rand(50, 50, field='complex')
            >>> Fc = F.conjugate()

        <b>See also</b> Faust.transpose, Faust.getH, Faust.H
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
            >>> F = rand(50, 50, field='complex')
            >>> H1 = F.getH()
            >>> H2 = F.transpose()
            >>> H2 = H2.conj()
            >>> (H1.toarray() == H2.toarray()).all()
            True

        <b>See also</b> Faust.transpose, Faust.conj, Faust.H
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
            >>> F = rand(50, 50, field='complex')
            >>> H1 = F.H
            >>> H2 = F.transpose()
            >>> H2 = H2.conj()
            >>> (H1.toarray() == H2.toarray()).all()
            True

        <b>See also</b> Faust.transpose, Faust.conj, Faust.getH
        """
        return F.getH()

    def pruneout(F, thres=0, **kwargs):
        """Returns a Faust optimized by removing useless zero rows and columns as many times as needed.

        Args:
            F: the Faust to optimize.
            thres: the threshold of number of nonzeros under what the
            rows/columns are removed.

        Returns:
            The optimized Faust.

        Example:
			>>> from pyfaust import rand
			>>> F = rand(1024, 1024, dim_sizes=[1, 1024], num_factors=32, fac_type='mixed')
			>>> pF = F.pruneout()
			>>> F.nbytes
			42875160
			>>> pF.nbytes
			38594468

        <b>See also</b> Faust.optimize
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
        """Returns a str object representing the Faust object.
        This method overloads a Python function.

        NOTE: Ideally this function is intended to return a valid Python
        expression but here this is not the case. Only information is
        displayed.

        Args:
            F: the Faust object.

        Examples:
            >>> from pyfaust import rand
            >>> F = rand(50, 50)
            >>> F.__repr__()
            >>> F.__repr__()
            'Faust size 50x50, density 0.5, nnz_sum 1250, 5 factor(s): \n-
            FACTOR 0 (double) SPARSE, size 50x50, density 0.1, nnz 250\n-
            FACTOR 1 (double) SPARSE, size 50x50, density 0.1, nnz 250\n-
            FACTOR 2 (double) SPARSE, size 50x50, density 0.1, nnz 250\n-
            FACTOR 3 (double) SPARSE, size 50x50, density 0.1, nnz 250\n-
            FACTOR 4 (double) SPARSE, size 50x50, density 0.1, nnz 250\n'
            >>> # the same function is called when typing F in a terminal (but
            >>> # the output is properly formatted with line feeds):
            >>> F
            Faust size 50x50, density 0.5, nnz_sum 1250, 5 factor(s):
                - FACTOR 0 (double) SPARSE, size 50x50, density 0.1, nnz 250
                - FACTOR 1 (double) SPARSE, size 50x50, density 0.1, nnz 250
                - FACTOR 2 (double) SPARSE, size 50x50, density 0.1, nnz 250
                - FACTOR 3 (double) SPARSE, size 50x50, density 0.1, nnz 250
                - FACTOR 4 (double) SPARSE, size 50x50, density 0.1, nnz 250

        <b>See also</b> Faust.nnz_sum, Faust.rcg, Faust.shape, Faust.factors,
        <b/>Faust.numfactors, Faust.display

        """
        _str = str(F.m_faust.to_string())
        return _str

    def __str__(F):
        """
        Converts F to its str representation.

        Returns:
            The str conversion of F.

        """
        return F.__repr__()

    def display(F):
        """
        Displays information about F.

        Args:
            F: the Faust object.

        Examples:
            >>> from pyfaust import rand
            >>> F = rand(50, 100, [1, 2], [50, 100], .5)
            >>> F.display()
            Faust size 50x100, density 0.94, nnz_sum 4700, 2 factor(s): <br/>
            FACTOR 0 (real) SPARSE, size 50x63, density 0.492063, nnz 1550 <br/>
            FACTOR 1 (real) SPARSE, size 63x100, density 0.5, nnz 3150 <br/>
            >>> F
            Faust size 50x100, density 0.94, nnz_sum 4700, 2 factor(s): <br/>
            FACTOR 0 (real) SPARSE, size 50x63, density 0.492063, nnz 1550 <br/>
            FACTOR 1 (real) SPARSE, size 63x100, density 0.5, nnz 3150 <br/>
            >>> print(F)
            Faust size 50x100, density 0.94, nnz_sum 4700, 2 factor(s): <br/>
            FACTOR 0 (real) SPARSE, size 50x63, density 0.492063, nnz 1550 <br/>
            FACTOR 1 (real) SPARSE, size 63x100, density 0.5, nnz 3150 <br/>

        <b>See also</b> Faust.nnz_sum, Faust.density, Faust.shape, Faust.factors,
        <b/>Faust.numfactors, Faust.__repr__

        """
        print(F.__repr__())
        #F.m_faust.display()

    def __add__(F, *args, **kwargs):
        """
        Sums F to one or a sequence of variables (Faust objects, arrays or scalars).


        NOTE: This method overloads the Python function/operator +.

        Args:
            args: the list of variables to sum all together with F.
            These can be Faust objects, numpy arrays (and
            scipy.sparse.csr_matrix) or scalars.
            Faust and arrays/matrices must be of compatible sizes.

        Returns:
            the sum as a Faust object.

        Example:
			>>> import pyfaust
			>>> from numpy.linalg import norm
			>>> F = pyfaust.rand(10, 12)
			>>> F
			Faust size 10x12, density 2.08333, nnz_sum 250, 5 factor(s):
			- FACTOR 0 (double) SPARSE, size 10x11, density 0.454545, nnz 50
			- FACTOR 1 (double) SPARSE, size 11x10, density 0.5, nnz 55
			- FACTOR 2 (double) SPARSE, size 10x11, density 0.454545, nnz 50
			- FACTOR 3 (double) SPARSE, size 11x10, density 0.5, nnz 55
			- FACTOR 4 (double) SPARSE, size 10x12, density 0.333333, nnz 40

			>>> G = pyfaust.rand(10, 12)
			>>> G
			Faust size 10x12, density 2.025, nnz_sum 243, 5 factor(s):
			- FACTOR 0 (double) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 1 (double) SPARSE, size 10x11, density 0.454545, nnz 50
			- FACTOR 2 (double) SPARSE, size 11x10, density 0.5, nnz 55
			- FACTOR 3 (double) SPARSE, size 10x12, density 0.333333, nnz 40
			- FACTOR 4 (double) SPARSE, size 12x12, density 0.333333, nnz 48

			>>> F+G
			Faust size 10x12, density 4.475, nnz_sum 537, 7 factor(s):
			- FACTOR 0 (double) SPARSE, size 10x20, density 0.1, nnz 20
			- FACTOR 1 (double) SPARSE, size 20x21, density 0.238095, nnz 100
			- FACTOR 2 (double) SPARSE, size 21x21, density 0.238095, nnz 105
			- FACTOR 3 (double) SPARSE, size 21x21, density 0.238095, nnz 105
			- FACTOR 4 (double) SPARSE, size 21x22, density 0.205628, nnz 95
			- FACTOR 5 (double) SPARSE, size 22x24, density 0.166667, nnz 88
			- FACTOR 6 (double) SPARSE, size 24x12, density 0.0833333, nnz 24

			>>> norm((F+G).toarray() - F.toarray() - G.toarray())
			8.09693699147347e-15
			>>> F+2
			Faust size 10x12, density 2.88333, nnz_sum 346, 7 factor(s):
			- FACTOR 0 (double) SPARSE, size 10x20, density 0.1, nnz 20
			- FACTOR 1 (double) SPARSE, size 20x21, density 0.142857, nnz 60
			- FACTOR 2 (double) SPARSE, size 21x20, density 0.154762, nnz 65
			- FACTOR 3 (double) SPARSE, size 20x21, density 0.142857, nnz 60
			- FACTOR 4 (double) SPARSE, size 21x11, density 0.281385, nnz 65
			- FACTOR 5 (double) SPARSE, size 11x24, density 0.19697, nnz 52
			- FACTOR 6 (double) SPARSE, size 24x12, density 0.0833333, nnz 24

			>>> norm((F+2).toarray() - F.toarray() - 2)
			3.580361673049448e-15
			>>> F+G+2
			Faust size 10x12, density 5.44167, nnz_sum 653, 9 factor(s):
			- FACTOR 0 (double) SPARSE, size 10x20, density 0.1, nnz 20
			- FACTOR 1 (double) SPARSE, size 20x30, density 0.05, nnz 30
			- FACTOR 2 (double) SPARSE, size 30x31, density 0.11828, nnz 110
			- FACTOR 3 (double) SPARSE, size 31x31, density 0.119667, nnz 115
			- FACTOR 4 (double) SPARSE, size 31x31, density 0.119667, nnz 115
			- FACTOR 5 (double) SPARSE, size 31x32, density 0.105847, nnz 105
			- FACTOR 6 (double) SPARSE, size 32x25, density 0.1225, nnz 98
			- FACTOR 7 (double) SPARSE, size 25x24, density 0.06, nnz 36
			- FACTOR 8 (double) SPARSE, size 24x12, density 0.0833333, nnz 24

			>>> norm((F+G+2).toarray() - F.toarray() - G.toarray() - 2)
			1.0082879011490611e-14
			>>> F+G+F+2+F
			Faust size 10x12, density 11.175, nnz_sum 1341, 13 factor(s):
			- FACTOR 0 (double) SPARSE, size 10x20, density 0.1, nnz 20
			- FACTOR 1 (double) SPARSE, size 20x30, density 0.05, nnz 30
			- FACTOR 2 (double) SPARSE, size 30x40, density 0.0333333, nnz 40
			- FACTOR 3 (double) SPARSE, size 40x50, density 0.025, nnz 50
			- FACTOR 4 (double) SPARSE, size 50x51, density 0.0509804, nnz 130
			- FACTOR 5 (double) SPARSE, size 51x52, density 0.0659879, nnz 175
			- FACTOR 6 (double) SPARSE, size 52x51, density 0.0678733, nnz 180
			- FACTOR 7 (double) SPARSE, size 51x54, density 0.0744372, nnz 205
			- FACTOR 8 (double) SPARSE, size 54x54, density 0.0713306, nnz 208
			- FACTOR 9 (double) SPARSE, size 54x36, density 0.063786, nnz 124
			- FACTOR 10 (double) SPARSE, size 36x34, density 0.0743464, nnz 91
			- FACTOR 11 (double) SPARSE, size 34x24, density 0.0784314, nnz 64
			- FACTOR 12 (double) SPARSE, size 24x12, density 0.0833333, nnz 24

			>>> norm((F+G+F+2+F).toarray() - 3*F.toarray() - G.toarray() - 2)
			2.892210888531005e-14


        See also Faust.__sub__
        """
        array_types = (np.ndarray,
                       scipy.sparse.csr_matrix,
                       scipy.sparse.csc_matrix)
        dev = F.device
        brd_err = ("operands could not be broadcast"
                   " together with shapes ")
        # handle possible broadcasting of F
        # if dimensions are inconsistent it'll fail later
        if(F.shape[1] == 1):
            max_ncols = np.max([a.shape[1] for a in args if isinstance(a,
                                                                       (Faust,
                                                                        *array_types))])
            F = F@Faust(np.ones(1, max_ncols, dtype=F.dtype), dev=F.device)
            if(F.shape[0] == 1):
                max_nrows = np.max([a.shape[0] for a in args if isinstance(a,
                                                                           (Faust,
                                                                            *array_types))])
            F = Faust(np.ones(max_nrows, 1, dtype=F.dtype), dev=F.device)@F
        def scalar2Faust(G):
            if isinstance(G, int):
                G = float(G)
            elif not isinstance(G, (np.float, np.complex)):
                raise TypeError("scalar must be int, float or complex")
            G, Gdtype = (float(G), np.float) if (isinstance(G, np.float) and
                                                 F.dtype != 'complex') else (complex(G), np.complex)
            return Faust([np.ones((F.shape[0], 1), dtype=F.dtype)*G,
                          np.ones((1, F.shape[1]), dtype=F.dtype).astype(Gdtype)],
                         dev=F.device)
        def broadcast_to_F(G):
            if G.shape[0] == 1:
                if G.shape[1] != F.shape[1]:
                    raise ve
                G = Faust(np.ones((F.shape[0], 1), dtype=F.dtype), dev=F.device) @ G
            elif G.shape[1] == 1:
                if G.shape[0] != F.shape[0]:
                    raise ve
                G = G @ Faust(np.ones((1, F.shape[1]), dtype=F.dtype), dev=F.device)
            return G
        # prepare the list of Fausts
        largs = []
        for i in range(0,len(args)):
            G = args[i]
            if isinstance(G,Faust):
                ve = ValueError(brd_err, F.shape,
                                G.shape, " argument i=", i)
                G = broadcast_to_F(G)
                if F.shape != G.shape:
                    raise Exception('Dimensions must agree, argument i=', i)
            elif isinstance(G,
                            array_types):
                if(G.size == 1):
                    G = scalar2Faust(G.reshape(1,)[0])
                if G.ndim == 1:
                    G = Faust([np.ones((F.shape[0], 1), dtype=F.dtype), G.reshape(1, G.size)], dev=F.device)
                else:
                    G = broadcast_to_F(Faust(G, dev=F.device))
            elif isinstance(G,(int, float, np.complex)):
                G = scalar2Faust(G)
            largs.append(G)

        id = seye(F.shape[1], dtype=F.dtype)
        id_vstack = svstack([id for i in range(0,
                                               len(largs)+1)])
        C = F.concatenate(*largs, axis=1)
        F = C@Faust(id_vstack, dev=F.device)
        return F

    def __radd__(F,lhs_op):
        """Returns lhs_op+F.
        <b>See also</b> Faust.__add__
        """
        return F.__add__(lhs_op)

    def __sub__(F,*args):
        """
        Subtracts from F one or a sequence of variables. Faust objects, arrays or scalars.

        NOTE: This method overloads the Python function/operator -.

        Args:
            args: the list of variables to compute the difference with F. These can
            be Faust objects, arrays (and scipy.sparse.csr_matrix) or scalars.
            Faust and arrays/matrices must be of compatible sizes.

        Returns:
            the difference as a Faust object.

        Example:
            >>> import pyfaust
            >>> from numpy.linalg import norm
            >>> F = pyfaust.rand(10, 12)
            >>> F
            Faust size 10x12, density 2.075, nnz_sum 249, 5 factor(s):
                - FACTOR 0 (double) SPARSE, size 10x12, density 0.333333, nnz 40
                - FACTOR 1 (double) SPARSE, size 12x11, density 0.454545, nnz 60
                - FACTOR 2 (double) SPARSE, size 11x10, density 0.5, nnz 55
                - FACTOR 3 (double) SPARSE, size 10x11, density 0.454545, nnz 50
                - FACTOR 4 (double) SPARSE, size 11x12, density 0.333333, nnz 44

            >>> G = pyfaust.rand(10, 12)
            >>> G
            Faust size 10x12, density 1.96667, nnz_sum 236, 5 factor(s):
                - FACTOR 0 (double) SPARSE, size 10x12, density 0.333333, nnz 40
                - FACTOR 1 (double) SPARSE, size 12x10, density 0.5, nnz 60
                - FACTOR 2 (double) SPARSE, size 10x12, density 0.333333, nnz 40
                - FACTOR 3 (double) SPARSE, size 12x12, density 0.333333, nnz 48
                - FACTOR 4 (double) SPARSE, size 12x12, density 0.333333, nnz 48

            >>> F-G
            Faust size 10x12, density 4.40833, nnz_sum 529, 7 factor(s):
                - FACTOR 0 (double) SPARSE, size 10x20, density 0.1, nnz 20
                - FACTOR 1 (double) SPARSE, size 20x24, density 0.166667, nnz 80
                - FACTOR 2 (double) SPARSE, size 24x21, density 0.238095, nnz 120
                - FACTOR 3 (double) SPARSE, size 21x22, density 0.205628, nnz 95
                - FACTOR 4 (double) SPARSE, size 22x23, density 0.193676, nnz 98
                - FACTOR 5 (double) SPARSE, size 23x24, density 0.166667, nnz 92
                - FACTOR 6 (double) DENSE, size 24x12, density 0.0833333, nnz 24

            >>> norm((F-G).toarray() - F.toarray() + G.toarray())
            2.8527547341127496e-15
            >>> F-2
            Faust size 10x12, density 2.875, nnz_sum 345, 7 factor(s):
                - FACTOR 0 (double) SPARSE, size 10x20, density 0.1, nnz 20
                - FACTOR 1 (double) SPARSE, size 20x22, density 0.113636, nnz 50
                - FACTOR 2 (double) SPARSE, size 22x21, density 0.151515, nnz 70
                - FACTOR 3 (double) SPARSE, size 21x20, density 0.154762, nnz 65
                - FACTOR 4 (double) SPARSE, size 20x12, density 0.25, nnz 60
                - FACTOR 5 (double) SPARSE, size 12x24, density 0.194444, nnz 56
                - FACTOR 6 (double) DENSE, size 24x12, density 0.0833333, nnz 24

            >>> norm((F-2).toarray() - F.toarray() + 2)
            0.0
            >>> F-G-2
            Faust size 10x12, density 5.375, nnz_sum 645, 9 factor(s):
                - FACTOR 0 (double) SPARSE, size 10x20, density 0.1, nnz 20
                - FACTOR 1 (double) SPARSE, size 20x30, density 0.05, nnz 30
                - FACTOR 2 (double) SPARSE, size 30x34, density 0.0882353, nnz 90
                - FACTOR 3 (double) SPARSE, size 34x31, density 0.12334, nnz 130
                - FACTOR 4 (double) SPARSE, size 31x32, density 0.105847, nnz 105
                - FACTOR 5 (double) SPARSE, size 32x33, density 0.102273, nnz 108
                - FACTOR 6 (double) SPARSE, size 33x25, density 0.123636, nnz 102
                - FACTOR 7 (double) SPARSE, size 25x24, density 0.06, nnz 36
                - FACTOR 8 (double) DENSE, size 24x12, density 0.0833333, nnz 24

            >>> norm((F-G-2).toarray() - F.toarray() + G.toarray() + 2)
            3.748544367384395e-15
            >>> F-G-F-2-F
            Faust size 10x12, density 11.0917, nnz_sum 1331, 13 factor(s):
                - FACTOR 0 (double) SPARSE, size 10x20, density 0.1, nnz 20
                - FACTOR 1 (double) SPARSE, size 20x30, density 0.05, nnz 30
                - FACTOR 2 (double) SPARSE, size 30x40, density 0.0333333, nnz 40
                - FACTOR 3 (double) SPARSE, size 40x50, density 0.025, nnz 50
                - FACTOR 4 (double) SPARSE, size 50x54, density 0.0407407, nnz 110
                - FACTOR 5 (double) SPARSE, size 54x53, density 0.0628931, nnz 180
                - FACTOR 6 (double) SPARSE, size 53x53, density 0.0622998, nnz 175
                - FACTOR 7 (double) SPARSE, size 53x55, density 0.0696398, nnz 203
                - FACTOR 8 (double) SPARSE, size 55x56, density 0.0688312, nnz 212
                - FACTOR 9 (double) SPARSE, size 56x35, density 0.0678571, nnz 133
                - FACTOR 10 (double) SPARSE, size 35x35, density 0.0702041, nnz 86
                - FACTOR 11 (double) SPARSE, size 35x24, density 0.0809524, nnz 68
                - FACTOR 12 (double) DENSE, size 24x12, density 0.0833333, nnz 24

            >>> norm((F-G-F-2-F).toarray() - F.toarray() + 2*F.toarray() + G.toarray() + 2)
            1.2631196225257117e-14

        <b>See also</b> Faust.__add__
        """

        nargs = []
        for arg in args:
            nargs += [ arg*-1 ]
        return F.__add__(*nargs)

    def __rsub__(F,lhs_op):
        """Returns lhs_op-F.
        <b>See also</b> Faust.__sub__
        """
        return F.__mul__(-1).__radd__(lhs_op)

    def __truediv__(F,s):
        """Divides F by the scalar s.
        This method overloads the Python function/operator `/' (whether s is a float or an integer).

        Args:
        F: the Faust object.
        s: the scalar to divide the Faust object with.

        Returns:
            the division result as a Faust object.

        <b>See also</b> Faust.__mul__, Faust.__itruediv__
        """
        if isinstance(s, (float, np.complex, int)) or isinstance(s, np.ndarray):
            return F*(1./s)
        else:
            raise Exception("unsupported operand type(s) for /: a Faust can only be "
                  "divided by a scalar.")

    def __itruediv__(F, s):
        """Divides F by the scalar s inplace.
        This method overloads the Python function/operator `/=' (whether s is a
        float or an integer).

        Args:
        F: the Faust object.
        s: the scalar to divide the Faust object with.

        Returns:
            the division result as a Faust object.

        <b>See also</b> Faust.__mul__, Faust.__truediv__

        """
        F = F/s
        return F

    def __imatmul__(F, A):
        """
        Inplace operation: F = F @ A
        """
        F = F@A
        return F

    def __imul__(F, A):
        """
        Inplace operation: F = F * A
        """
        F = F*A
        return F

    def __iadd__(F, A):
        """
        Inplace operation: F = F + A
        """
        F = F+A
        return F

    def __isub__(F, A):
        """
        Inplace operation: F = F - A
        """
        F = F-A
        return F

    def __matmul__(F, A):
        """Multiplies F by A which is a dense numpy.matrix/numpy.ndarray or a Faust object.
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
            Faust-dense_matrix multiplication. In some cases though, it stays quicker: moreover when the Faust is composed of a small number of factors.

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
            >>> F = rand(50, 100)
            >>> A = np.random.rand(F.shape[1], 50)
            >>> B = F@A # == F*A or F.dot(A)
            >>> # is equivalent to B = F.__matmul__(A)
            >>> G = rand(F.shape[1], 5)
            >>> H = F@G
            >>> # H is a Faust because F and G are

        <b>See also</b> Faust.__init__, Faust.rcg, Faust.__mul__, Faust.__matmul__, Faust.dot
        """
        if(isinstance(A, Faust)):
            if(F.shape[1] != A.shape[0]): raise ValueError("The dimensions of "
                                                          "the two Fausts must "
                                                          "agree.")
            if F.dtype == np.complex and A.dtype != np.complex:
                A = A.astype(np.complex)
            elif F.dtype != np.complex and A.dtype == np.complex:
                F = F.astype(np.complex)
            return Faust(core_obj=F.m_faust.multiply_faust(A.m_faust))
        elif(isinstance(A, (float, int, np.complex))):
            raise ValueError("Scalar operands are not allowed, use '*'"
                             " instead")
        elif isinstance(A, np.ndarray):
            w = lambda: warnings.warn("The numpy array to multiply is not "
                                      "F_CONTIGUOUS, costly conversion on the "
                                      "fly. Please use np.asfortranarray by "
                                      "yourself.")
            if A.dtype == np.complex and F.dtype != np.complex:
                A_r = np.asfortranarray(A.real)
                A_i = np.asfortranarray(A.imag)
                if not A.imag.flags['F_CONTIGUOUS']:
                    w()
                    A_i = np.asfortranarray(A_i)
                    A_r = np.asfortranarray(A_r)
                G = F.m_faust.multiply(A_r) + 1j*F.m_faust.multiply(A_i)
                return G
            else:
                if not A.flags['F_CONTIGUOUS']:
                    w()
                    A = np.asfortranarray(A)
                if F.dtype == 'complex' and A.dtype != 'complex':
                    A = A.astype('complex')
                return F.m_faust.multiply(A)
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
            raise TypeError("can't multiply a Faust by something that is not a Faust, a np.ndarray, a csr_matrix or a dia_matrix.")

    def dot(F, A, *args, **kwargs):
        """Performs equivalent operation of numpy.dot() between the Faust F and A.
        More specifically:
            - Scalar multiplication if A is a scalar but F*A is preferred.
            - Matrix multiplication if A is a Faust or numpy.ndarray/numpy.matrix but F @ A is preferred.

        <b>See also</b> Faust.__init__, Faust.rcg, Faust.__mul__, Faust.__matmul__, Faust.dot
        """
        if(isinstance(A, (float, int, np.complex))):
            return F*A
        return F.__matmul__(A)

    def matvec(F, x):
        """
        This function implements the scipy.sparse.linalg.LinearOperator.matvec function such that scipy.sparse.linalg.aslinearoperator function works on a Faust object.

        Example:
            >>> import numpy as np
            >>> import pyfaust as pf
            >>> F = pf.rand(10, 10)
            >>> x = np.random.rand(F.shape[1])
            >>> F@x
            array([ 69.88213474,  80.99557416, 106.05274531,  55.55415099, 127.4286391 ,  43.72638629,  90.05175021, 105.05742936, 78.14695324,  98.79773085])
            >>> F.matvec(F, x)
            array([ 69.88213474,  80.99557416, 106.05274531,  55.55415099, 127.4286391 ,  43.72638629,  90.05175021, 105.05742936, 78.14695324,  98.79773085])


        <b>See also</b> Faust.dot, Faust.__matmul__
        """
        return F.dot(x)

    def __mul__(F, A):
        """Multiplies the Faust F by A.
        This method overloads a Python function/operator (*).

        More specifically:
        - It's a scalar multiplication if A is a scalar number.
        - if A is a vector of appropriate size, the function performs a Faust-vector broadcasting.

        Args:
            F: the Faust object.
            A: is a scalar number or a (vector) numpy.ndarray

        Returns: The result of the multiplication as a Faust (either A is a
        vector or a scalar).

        Raises: TypeError
            If the type of A is not supported:
            <code>
            >>> import numpy
            >>> F*'a'
            TypeError: * use is forbidden in this case. It is allowed only for Faust-scalar multiplication or Faust vector broadcasting.
            </code>

        Examples:
            >>> from pyfaust import rand
            >>> import numpy as np
            >>> F = rand(50, 100)
            >>> v = np.random.rand(F.shape[1])
            >>> B = F*v # Faust vector broadcasting
            >>> # is equivalent to B = F.__mul__(v)
            >>> F_times_two = F*2

        <b>See also</b> Faust.__rmul__
        """
        if(isinstance(A, (float, int, np.complex))):
            if isinstance(A, int):
                A = float(A)
            elif isinstance(A, np.complex):
                if F.dtype != np.complex:
                    F = F.astype(np.complex)
            else:
                # A is a float
                if F.dtype == np.complex:
                    A = np.complex(A)
            return Faust(core_obj=F.m_faust.multiply_scal(A))
        elif(isinstance(A, np.ndarray)):
            if(A.size == 1):
                if A.dtype == np.complex:
                    return F*(np.complex(A.squeeze()))
                else:
                    return F*(float(A.squeeze()))
            if A.ndim == 1 and A.size == F.shape[1] \
            or A.ndim == 2 and A.shape[0] == 1:
                return F@Faust(np.diag(A.squeeze()), dev=F.device)
        # A is a Faust, a numpy.ndarray (eg. numpy.matrix) or anything
        raise TypeError("* use is forbidden in this case. It is allowed only"
                        " for Faust-scalar multiplication or Faust vector"
                        " broadcasting.")

    def __rmul__(F, lhs_op):
        """ lhs_op*F
        <b>See also</b> Faust.__mul__
        """
        if isinstance(lhs_op, (float, int, np.complex)):
            return F.__mul__(lhs_op)
        elif(isinstance(lhs_op, np.ndarray)):
            #TODO: refactor with __mul__
            if(lhs_op.size == 1):
                if lhs_op.dtype == np.complex:
                    return F*(np.complex(lhs_op.squeeze()))
                else:
                    return F*(float(lhs_op.squeeze()))
            if lhs_op.ndim == 1 and lhs_op.size == F.shape[1] \
            or lhs_op.ndim == 2 and lhs_op.shape[0] == 1:
                return F@Faust(np.diag(lhs_op.squeeze()), dev=F.device)
        # A is a Faust, a numpy.ndarray (eg. numpy.matrix) or anything
        raise TypeError("* use is forbidden in this case. It is allowed only"
                        " for Faust-scalar multiplication or Faust vector"
                        " broadcasting.")

    def __rmatmul__(F, lhs_op):
        """Returns lhs_op.__matmul__(F).
        <b>See also Faust.__matmul__</b>

        Examples:
            >>> from pyfaust import rand
            >>> import numpy as np
            >>> F = rand(50, 100)
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
        """Concatenates F with len(args) Faust objects, numpy arrays or scipy sparse matrices.
        The resulting Faust:
                C = F.concatenate(G, H, ... ,axis)
        verifies that:
                C.toarray() == numpy.concatenate((F.toarray(), G.toarray(),
                H.toarray(),...), axis)

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
                >>> from pyfaust import rand
                >>> F = rand(2,51);
                >>> G = rand(2,25);
                >>> F.concatenate(G, 0)
                </code>
                ValueError
                Axis must be 0 or 1.
                <code>
                >>> from pyfaust import rand
                >>> F = rand(2,51);
                >>> G = rand(2,25);
                >>> F.concatenate(G, 5)
                </code>
                ValueError
                You can't concatenate a Faust with something that is not a
                Faust, a numpy array or a scipy sparse matrix.
                <code>
                >>> from pyfaust import rand
                >>> F = rand(2,51);
                >>> F.concatenate(['a', 'b', 'c'], 5)
                </code>

            Examples:
                >>> from pyfaust import rand
                >>> F = rand(50, 50)
                >>> G = rand(50, 50)
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
                Faust size 84x50, density 0.558333, nnz_sum 2345, 6 factor(s): <br/>
                FACTOR 0 (real) SPARSE, size 84x91, density 0.0549451, nnz 420<br/>
                FACTOR 1 (real) SPARSE, size 91x92, density 0.0543478, nnz 455<br/>
                FACTOR 2 (real) SPARSE, size 92x90, density 0.0555556, nnz 460<br/>
                FACTOR 3 (real) SPARSE, size 90x92, density 0.0543478, nnz 450<br/>
                FACTOR 4 (real) SPARSE, size 92x100, density 0.05, nnz 460<br/>
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
                Faust size 384x50, density 0.575521, nnz_sum 11050, 6 factor(s): <br/>
                FACTOR 0 (double) SPARSE, size 384x400, density 0.0224609, nnz 3450<br/>
                FACTOR 1 (double) SPARSE, size 400x400, density 0.01125, nnz 1800<br/>
                FACTOR 2 (double) SPARSE, size 400x400, density 0.01125, nnz 1800<br/>
                FACTOR 3 (double) SPARSE, size 400x400, density 0.01125, nnz 1800<br/>
                FACTOR 4 (double) SPARSE, size 400x400, density 0.01125, nnz 1800<br/>
                FACTOR 5 (double) SPARSE, size 400x50, density 0.02, nnz 400<br/>

        """
        if("axis" in kwargs.keys()):
            axis = kwargs['axis']
        else:
            axis = 0
        if(axis not in [0,1]):
            raise ValueError("Axis must be 0 or 1.")

        largs = []
        any_G_is_cplx = F.dtype == np.complex
        for i,G in enumerate(args):
            if(isinstance(G, (np.ndarray, csr_matrix, csc_matrix))):
                G = Faust(G, dev=F.device)
            elif(not isinstance(G, Faust)): raise ValueError("You can't concatenate a "
                                                           "Faust with something "
                                                           "that is not a Faust, "
                                                           "a numpy array or scipy "
                                                           "sparse matrix.")
            any_G_is_cplx |= G.dtype == np.complex
            largs.append(G)

            if(axis == 0 and F.shape[1] != G.shape[1] or axis == 1 and F.shape[0]
               != G.shape[0]): raise ValueError("The dimensions of "
                                                "the two Fausts must "
                                                "agree.")


        if any_G_is_cplx:
            # one Faust is complex convert all real Faust to complex
            for i in range(len(largs)):
                if largs[i].dtype != np.complex:
                    largs[i] = largs[i].astype(np.complex)
            if F.dtype != np.complex:
                F = F.astype(np.complex)


        if all([isFaust(G) for G in largs]) and not "iterative" in kwargs.keys() or kwargs['iterative']:
            # use iterative meth.
            if axis == 0:
                C = Faust(core_obj=F.m_faust.vertcatn([G.m_faust for G in largs]))
            else: # axis == 1
                C = Faust(core_obj=F.m_faust.horzcatn([G.m_faust for G in largs]))
            return C

        # use recursive meth.
        C=F
        for G in args:
           if axis == 0:
                C = Faust(core_obj=C.m_faust.vertcat(G.m_faust))
           elif axis == 1:
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
            >>> F = rand(10**5, 10**5, 2, 10**5, density=10**-4, fac_type='sparse')
            >>> F
            FACTOR 0 (real) SPARSE, size 100000x100000, density 9e-05, nnz 900000<br/>
            FACTOR 1 (real) SPARSE, size 100000x100000, density 9e-05, nnz 900000<br/>
            >>> # an attempt to convert F to a dense matrix is most likely to raise a memory error
            >>> # the sparse format is the only way to handle such a large Faust
            >>> F.toarray()
            ...
            MemoryError

        <b>See also</b> Faust.todense
        """
        return F.m_faust.get_product()


    def todense(F):
        """Returns a numpy matrix for the full matrix implemented by F.

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

        WARNING: this function is deprecated and might be deleted in future
        versions of pyfaust. Please use Faust.toarray instead.

        Examples:
            >>> from pyfaust import rand
            >>> F = rand(10**5, 10**5, 2, 10**5, density=10**-4, fac_type='sparse')
            >>> F
            FACTOR 0 (real) SPARSE, size 100000x100000, density 9e-05, nnz 900000<br/>
            FACTOR 1 (real) SPARSE, size 100000x100000, density 9e-05, nnz 900000<br/>
            >>> # an attempt to convert F to a dense matrix is most likely to raise a memory error
            >>> # the sparse format is the only way to handle such a large Faust
            >>> F.toarray()
            (...)
            MemoryError

        <b>See also</b> Faust.toarray
        """
        warnings.warn("Faust.todense() is deprecated and will be deleted in the "
             "future.")
        return np.matrix(F.toarray())


    def __getitem__(F, indices):
        """Indexes or slices a Faust.
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

            >>> F = rand(50, 100)
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
        if(isinstance(indices, np.ndarray)):
            indices = list(indices)
        if(indices == Ellipsis): # F[...]
            out_indices = [slice(0,F.shape[0]), slice(0, F.shape[1])]
        elif(isinstance(indices,int)): # F[i] # a line
            out_indices = [slice(indices, indices+1), slice(0, F.shape[1])]
        elif(isinstance(indices,slice)):
            #F[i:j] a group of contiguous lines
            out_indices = [indices, slice(0, F.shape[1])]
        elif(isinstance(indices, list)):
            out_indices = [indices, slice(0, F.shape[1])]
            #TODO: check indices are all integers lying into F shape
        elif(isinstance(indices, tuple)):
            if(len(indices) == 1):
                return F.__getitem__(indices[0])
            if(len(indices) == 2):
                out_indices = [0,0]
                if(isinstance(indices[0], int) and isinstance(indices[1],int)):
                    if 'OPT_GET_ITEM' in environ and environ['OPT_GET_ITEM'] == '0':
                        return F.toarray()[indices[0],indices[1]]
                    else:
                        return F.m_faust.get_item(indices[0], indices[1])
                if(isinstance(indices[0], np.ndarray)):
                    indices = (list(indices[0]), indices[1])
                if(isinstance(indices[1], np.ndarray)):
                    indices = (indices[0], list(indices[1]))
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
        """Gives the total number of non-zero elements in the factors of F.

        The function sums together the number of non-zero elements of
        each factor and returns the result. Note that for efficiency the sum is
        computed at Faust creation time and kept in cache.

        Returns:
            the number of non-zeros.

        Example:
			>>> import pyfaust as pf
			>>> F = pf.rand(1024, 1024)
			>>> F
			Faust size 1024x1024, density 0.0244141, nnz_sum 25600, 5 factor(s): 
			- FACTOR 0 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120
			- FACTOR 1 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120
			- FACTOR 2 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120
			- FACTOR 3 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120
			- FACTOR 4 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120

			>>> F.nnz_sum()
			25600

        <b>See also</b> Faust.rcg, Faust.density.
        """
        return F.m_faust.nnz()

    def density(F):
        """ Calculates the density of F such that F.nnz_sum() == F.density()/F.size.

        NOTE: A value of density below one indicates potential memory savings
        compared to storing the corresponding dense matrix F.toarray(), as well
        as potentially faster matrix-vector multiplication when applying F @ x
        instead of F.toarray() @ x.

        NOTE: A density above one is possible but prevents any saving.

        Args:
            F: the Faust object.

        Returns:
            the density value (float).

        Examples:
        >>> from pyfaust import rand
        >>> F = rand(5, 50, density=.5)
        >>> dens = F.density()

        <b>See also</b> Faust.nnz_sum, Faust.rcg, Faust.size, Faust.toarray
        """
        return float(F.nnz_sum())/F.size


    def rcg(F):
        """Computes the Relative Complexity Gain.

        The RCG is the theoretical gain brought by the Faust representation relatively to its dense
        matrix equivalent. <br/>The higher is the RCG, the more computational
        savings are made.
        This gain applies both for storage space and computation time.

        NOTE: F.rcg() == 1/F.density()

        Args:
            F: the Faust object.

        Returns:
            the RCG value (float).

        >>> from pyfaust import rand
        >>> F = rand(1024, 1024)
        >>> F.rcg()
        40.96
        >>> F
        Faust size 1024x1024, density 0.0244141, nnz_sum 25600, 5 factor(s):
            - FACTOR 0 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 1 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 2 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 3 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 4 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120

        >>> F.size/F.nnz_sum()
        40.96


        <b>See also</b>: Faust.density, Faust.nnz_sum, Faust.shape.
        """
        d = F.density()
        if(d > 0):
            return 1/d
        elif(d == 0):
            return float("inf")
        else:
            return float(-1)

    def norm(F, ord='fro', **kwargs): #**kwargs):
        """Computes the norm of F.
        Several types of norm are available: 1-norm, 2-norm, inf-norm and Frobenius norm.
        The norm of F is equal to the numpy.linalg.norm of F.toarray().

        WARNING:
            The norm computation time can be expected to be of order
        n*min(F.shape) with n the time for multipliying F by a vector.
        Nevertheless, the implementation allows that memory usage remains
        controlled by avoiding to explicitly compute F.toarray(). Please pay
        attention to the full_array (and batch_size) arguments for a better
        understanding.

        Args:
            F: the Faust object.
            ord: (optional) the norm order (1, 2, numpy.inf) or "fro" for
            Frobenius norm (by default the Frobenius norm is computed).
            threshold: (optional) power iteration algorithm threshold (default
            to .001). Used only for norm(2).
            max_num_its: (optional) maximum number of iterations for
            power iteration algorithm (default to 100). Used only for norm(2).
            full_array: (optional) this argument applies only for 1-norm,
            inf-norm and Frobenius norm. If True the Faust full array
            is computed before computing the norm otherwise it is not. By
            default it is set to False. Many configurations exist in which
            full_array == False can be more efficient but it needs to
            finetune the batch_size argument.
            batch_size: (optional) this argument applies only when the
            full_array argument is set to False (for the 1-norm, inf-norm and
            Frobenius norm). It determines the number of Faust columns (resp. rows)
            that are built in memory in order to compute the Frobenius norm and
            the 1-norm (resp. the inf-norm). This parameter is primary in the
            efficiency of the computation and memory consumption. By  default,
            it is set to 1 (which is certainly not the optimal configuration in
            many cases in matter of computation time but always the best in
            term of memory cost).

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
            >>> F = rand(50, 100, [1, 2], density=.5)
            >>> F.norm()
            23.55588891399667
            >>> F.norm(2)
            5.929720822717308
            >>> F.norm('fro')
            23.55588891399667
            >>> F.norm(np.inf)
            18.509101197826254

        <b>See also</b>: Faust.normalize
        """
        if ord not in [1, 2, "fro", np.inf]:
            raise ValueError("ord must have the value 1, 2, 'fro' or numpy.inf.")
        return F.m_faust.norm(ord, **kwargs)

    def power_iteration(self, threshold=1e-3, maxiter=100):
        """Performs the power iteration algorithm to compute the greatest eigenvalue of the Faust.
        For the algorithm to succeed the Faust should be diagonalizable
        (similar to a digonalizable Faust), ideally, a symmetric positive-definite Faust.

        Args:
            threshold: the precision required on the eigenvalue.
            maxiter: the number of iterations above what the algorithm will stop anyway.

        Returns:
            The greatest eigenvalue approximate.

        Examples:
            >>> from pyfaust import rand
            >>> F = rand(50, 50)
            >>> F = F@F.H
            >>> power_iteration(F)
            14630.932668438209

        """
        return self.m_faust.power_iteration(threshold=threshold,
                                            max_num_its=maxiter)

    def normalize(F, ord='fro', axis=1):
        """Normalizes F along the axis dimension using the ord-norm.
        The function is able to normalize the columns (default axis=1):

            NF = F.normalize(ord) is such that for all i in range(0,F.shape[1]) NF.toarray()[:,i] == F.toarray()[:,i]/norm(F.toarray(), ord)

        Likewise for normalization of the rows (axis=0):

            NF = F.normalize(ord, axis=0) is such that for all i in range(0,F.shape[0]) NF.toarray()[i,:] == F.toarray()[i,:]/norm(F.toarray(), ord)

        The variable ord designates one of the Faust.norm() compatible norms.

        Args:
            ord: the norm order to use (see Faust.norm).
            axis: if 1 the columns are normalized, if 0 the rows.

        Returns:
            the normalized Faust

        Example:
			>>> from numpy.linalg import norm
			>>> import numpy as np
			>>> import pyfaust as pf
			>>> F = pf.rand(10, 10)
			>>> F
			Faust size 10x10, density 2.5, nnz_sum 250, 5 factor(s):
			- FACTOR 0 (double) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 1 (double) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 2 (double) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 3 (double) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 4 (double) SPARSE, size 10x10, density 0.5, nnz 50

			>>> # normalize columns according to default fro-norm/2-norm
			>>> # then test the second column is properly normalized
			>>> nF2 = F.normalize()
			>>> nF2
			Faust size 10x10, density 2.5, nnz_sum 250, 5 factor(s):
			- FACTOR 0 (double) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 1 (double) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 2 (double) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 3 (double) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 4 (double) SPARSE, size 10x10, density 0.5, nnz 50

			>>> norm(nF2[:,1].toarray() - F[:, 1].toarray()/F[:,1].norm())
			1.2719202621569003e-16
			>>>
			>>> # this time normalize rows using 1-norm and test the third row
			>>> nF1 = F.normalize(ord=1, axis=0)
			>>> nF1
			Faust size 10x10, density 2.5, nnz_sum 250, 5 factor(s):
			- FACTOR 0 (double) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 1 (double) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 2 (double) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 3 (double) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 4 (double) SPARSE, size 10x10, density 0.5, nnz 50

			>>> norm(nF1[2, :].toarray() - F[2, :].toarray()/F[2, :].norm(ord=1))
			2.5438405243138006e-16
			>>>
			>>> # and the same with inf-norm
			>>> nFinf = F.normalize(ord=np.inf, axis=0)
			>>> norm(nFinf[2, :].toarray() - F[2, :].toarray()/F[2, :].norm(ord=np.inf))
			4.3885418357208765e-17


        <b>See also</b>: Faust.norm
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
        >>> F = rand(100,100, num_factors=2, density=.5)
        >>> nf = F.numfactors()
        >>> nf
        2
        >>> nf == len(F)
        True

        <b>See also</b> Faust.factors, Faust.__len__
        """
        return F.m_faust.get_nb_factors()

    def __len__(F):
        """
        Returns the number of factors of F.

        Returns:
            the number of factors.

        Examples:
        >>> from pyfaust import rand
        >>> F = rand(50, 100)
        >>> nf = F.numfactors()
        >>> nf
        2
        >>> nf == len(F)
        True

        <b>See also</b> Faust.factors, Faust.numfactors
        """
        return F.numfactors()

    def factors(F, indices):
        """
        Returns the i-th factor of F or a new Faust composed of F factors whose indices are listed in indices.

        Note: Factors are copied in memory.

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


        Example:
            >>> from pyfaust import rand
            >>> F = rand(5, 10)
            >>> f0 = F.factors(0)
            >>> G = F.factors(range(3,5)) # a new Faust composed of the two last factors of F

        <b>See also</b> Faust.numfactors, Faust.transpose
        """
        if(hasattr(indices, '__iter__')):
            indices = list(indices)
        else:
            indices = list([indices])
        factors = []
        oi = None
        for i in indices:
            if(not isinstance(i, int)):
                raise TypeError("Index must be an integer.")
            if(oi != None and i-oi != 1):
                raise Exception("Indices must be contiguous.")
            factors += [F.m_faust.get_fact_opt(i)]
            oi = i
        if(len(factors) == 1):
            return factors[0]
        else:
            return pyfaust.Faust(factors, dev=F.device)

    def right(F, i):
        """
        Returns the right hand side factors of F from index i to F.numfactors()-1.

        Returns:
            a Faust if the number of factors to be returned is greater than 1, a
            numpy array or a sparse matrix otherwise.

        Examples:
            >>> from pyfaust import rand
            >>> F = rand(5, 10, 5)
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
            a Faust if the number of factors to be returned is greater than 1, a
            numpy array or a sparse matrix otherwise.

        Examples:
            >>> from pyfaust import rand
            >>> F = rand(5, 10, 5)
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
            >>> F = rand(5, 10)
            >>> f0 = F.factors(0)

        <b>See also</b> Faust.numfactors
        """
        fact = F.m_faust.get_fact(i)
        return fact

    def save(F, filepath, format="Matlab"):
        """Saves the Faust F into a file.
        The file is saved in Matlab format version 5 (.mat extension).

        NOTE: storing F should typically use rcg(F) times less disk space than
        storing F.toarray(). See Faust.nbytes for a precise size.

        Args:
            F: the Faust object.
            filepath: the path for saving the Faust (should end with .mat
            if Matlab format is used).
            format: (optional) The format to use for
            writing. By default, it's "Matlab" to save the Faust in a .mat
            file (currently only that format is available).

        Raises:
            ValueError.


        Example:
            >>> from pyfaust import rand, Faust
            >>> F = rand(5, 10, field='complex')
            >>> F.save("F.mat")
            >>> G = Faust(filepath="F.mat")

        <b>See also</b> Faust.__init__ Faust.rcg Faust.load Faust.load_native
        """
        if(format not in ["Matlab"]):
            raise ValueError("Only Matlab or Matlab_core format is supported.")
        if(format == "Matlab"):
            F.m_faust.save_mat_file(filepath)

    @staticmethod
    def load(filepath):
        """Loads a Faust from a MAT file.

        The format is Matlab format version 5 and the filepath should end with
        a .mat extension.

        The Faust must have been saved before with Faust.save.

        Args:
            filepath: the filepath of the .mat file.

        Example:
            >>> from pyfaust import rand, Faust
            >>> F = rand(5, 10, field='complex')
            >>> F.save("F.mat")
            >>> G = Faust.load(filepath="F.mat") # equiv. to Faust("F.mat")

        <b>See also</b> Faust.__init__, Faust.save
        """
        contents = loadmat(filepath)
        factors = contents['faust_factors'][0].tolist()
        # check if any BSR matrix exists here

        for i in range(len(factors)):
            if factors[i][0][0].dtype == '<U3' and factors[i][0][0][0] == 'bsr':
                nrows, ncols, bnnz = factors[i][0][1][0]
                bcolinds = factors[i][0][2][0]
                browptr = factors[i][0][3][0]
                bdata = factors[i][0][4][0]
                bnrows = int(nrows/(browptr.shape[0]-1))
                bncols = int(bdata.shape[0]/bnrows/bnnz)
                bdata_ = np.empty((bnnz, bnrows, bncols))
                for bi in range(bnnz): # .mat is in col-major order for blocks
                    bdata_[bi] = bdata.reshape(bnnz, bncols, bnrows)[bi].T
                # override the factor with the corresponding scipy bsr matrix
                factors[i] = bsr_matrix((bdata_, bcolinds, browptr), shape=(nrows,
                                                                            ncols))
        return Faust(factors)

    def load_native(filepath):
        """
        The format is Matlab format version 5 and the filepath should end with
        a .mat extension (native C++ version).

        The Faust must have been saved before with Faust.save.

        Args:
            filepath: the filepath of the .mat file.

        Example:
            >>> from pyfaust import rand, Faust
            >>> F = rand(5, 10, field='complex')
            >>> F.save("F.mat")
            >>> G = Faust.load_native(filepath="F.mat") # equiv. to Faust("F.mat")

        <b>See also</b> Faust.__init__, Faust.save
        """
        _type = _FaustCorePy.FaustCoreGenDblCPU.get_mat_file_type(filepath)
        if _type == -1:
            raise Exception("Invalid .mat file")
        elif _type == 0:
            F = Faust( core_obj=_FaustCorePy.FaustCoreGenFltCPU.read_from_mat_file(filepath))
        elif _type == 1:
            F = Faust( core_obj=_FaustCorePy.FaustCoreGenDblCPU.read_from_mat_file(filepath))
        elif _type == 2:
            F = Faust( core_obj=_FaustCorePy.FaustCoreGenCplxDblCPU.read_from_mat_file(filepath))
        return F

    def astype(F, dtype):
        """
        Converts F to the dtype passed as argument in a new Faust.

        Args:
            dtype: np.float or np.complex.

        Returns:
            A Faust copy of F converted to dtype.

        Example:
            >>> from pyfaust import rand
            >>> F = rand(10, 10, dtype='double')
            >>> F.astype('float32')
            Faust size 10x10, density 2.5, nnz_sum 250, 5 factor(s):
                - FACTOR 0 (float) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 1 (float) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 2 (float) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 3 (float) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 4 (float) SPARSE, size 10x10, density 0.5, nnz 50

            >>> F.astype("complex")
            Faust size 10x10, density 2.5, nnz_sum 250, 5 factor(s):
                - FACTOR 0 (complex) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 1 (complex) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 2 (complex) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 3 (complex) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 4 (complex) SPARSE, size 10x10, density 0.5, nnz 50

        """
        #TODO: full list of numpy args or **kw_unknown
        if np.dtype(dtype) not in [np.complex, np.float32, np.double]:
            raise TypeError("Faust supports only complex, double or float32 dtype-s")
        if dtype == F.dtype:
            return F.clone(dev=F.device)
        if F.dtype == np.complex:
            if dtype == np.float32:
                return F.real.astype('float32')
            else:
                # dtype is double
                return F.real
        else:
            return Faust([F.factors(i).astype(dtype) for i in
                          range(F.numfactors())], dev=F.device)

    @property
    def real(F):
        """
        Returns the real part of F as a Faust.
        """
        if F.dtype != np.complex:
            return F
        else:
            #  return 1/2 * (F + F.conj())
            return _cplx2real_op(F)[:F.shape[0],
                                     :F.shape[1]].pruneout()
        #return 1/2 * (F + F.conj())

    @property
    def imag(F):
        """
        Returns the imaginary part of F as a Faust.
        """
        if F.dtype != np.complex:
            # return Faust(csr_matrix(F.shape)) # TODO: debug pyx code
            return Faust(csr_matrix((np.array([0.]).astype(F.dtype),
                                     ([0],[0])), (F.shape)))
        else:
            # return 1/2j * (F + F.conj())
            return _cplx2real_op(F)[F.shape[0]:2*F.shape[0],
                                     :2*F.shape[1]].pruneout()


    def asarray(F, *args, **kwargs):
        return F

    @property
    def dtype(F):
        """
        Returns the dtype of the Faust.

        This function is intended to be used as a property (see the examples).

        Args:
            F: the Faust object.

        Returns:
            the dtype of F, which can be float32, float64 or complex128.

        Examples:
            >>> from pyfaust import rand
            >>> F = rand(5, 10, field='complex')
            >>> F.dtype
            dtype('complex128')
            >>> F = rand(5, 10)
            >>> F.dtype
            dtype('float64')
            >>> G = rand(5, 5, dtype='float32')
            >>> G.dtype
            dtype('float32')


        """
        return F.m_faust.dtype()

    def imshow(F, name='F'):
        """
        Displays image of F's full matrix and its factors.

        Args:
            F: the Faust object.
            name: (optional str) the displayed name on the plotted figure.


        Examples:
        >>> from pyfaust import rand
        >>> import matplotlib.pyplot as plt
        >>> F = rand(10, 20, density=.5, field='complex')
        >>> F.imshow()
        >>> plt.show()


        <b>See also</b> Faust.display
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
        """Computes the (Moore-Penrose) pseudo-inverse of F.toarray().

        Warning: this function makes a call to Faust.toarray().

        Returns:
            The dense pseudo-inverse matrix.

        Example:
            >>> from pyfaust import rand
            >>> import numpy as np
            >>> from numpy.linalg import pinv
            >>> F = rand(128, 32)
            >>> M = F.toarray()
            >>> np.allclose(F.pinv(), pinv(M))
            True

        See also <a href="https://numpy.org/doc/stable/reference/generated/numpy.linalg.pinv.html">numpy.linalg.pinv</a>
        """
        from numpy.linalg.linalg import pinv
        return pinv(F.toarray())

    @staticmethod
    def isFaust(obj):
        """
        Returns True if obj is a Faust object, False otherwise.

        Examples:
            >>> from pyfaust import *
            >>> Faust.isFaust(2) # isFaust(2) works as well
            False
            >>> Faust.isFaust(rand(5,10))
            True

        """
        return isinstance(obj, Faust)

    def issparse(F, csr=True, bsr=False):
        """
        Returns True if all F factors are sparse False otherwise.

        What a sparse factor is, depends on csr and bsr arguments.

        Args:
            csr: True to consider CSR matrices in F as sparse matrices, False otherwise.
            bsr: True to consider BSR matrices in F as sparse matrices, False otherwise.

        Example:
			>>> import pyfaust as pf
			>>> F = pf.rand(10, 10, fac_type='sparse')
			>>> F
			Faust size 10x10, density 2.5, nnz_sum 250, 5 factor(s):
			- FACTOR 0 (double) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 1 (double) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 2 (double) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 3 (double) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 4 (double) SPARSE, size 10x10, density 0.5, nnz 50

			>>> F.issparse()
			True
			>>> F.issparse(csr=True)
			True
			>>> F.issparse(csr=False)
			ValueError: It doesn't make sense to set csr=False and bsr=False as the function will always return False

			>>> F = pf.rand_bsr(10, 10, 2, 2, 2, .1)
			>>> F
			Faust size 10x10, density 0.24, nnz_sum 24, 2 factor(s):
			- FACTOR 0 (double) BSR, size 10x10 (blocksize = 2x2), density 0.12, nnz 12 (nnz blocks: 3)
			- FACTOR 1 (double) BSR, size 10x10 (blocksize = 2x2), density 0.12, nnz 12 (nnz blocks: 3)

			>>> F.issparse() # default config. recognizes only csr 
            False
            >>> F.issparse(bsr=True)
            True
            >>> F = pf.rand(10, 10, fac_type='dense')
            >>> F.issparse()
            False

        <b>See also</b> Faust.isdense, pyfaust.rand, pyfaust.rand_bsr
        """
        if not csr and not bsr:
            raise ValueError('It doesn\'t make sense to set csr=False and'
                             ' bsr=False as the function will always return'
                             ' False')
        return F.m_faust.is_all_sparse(csr, bsr)

    def isdense(F):
        """
        Returns True if all factors are dense arrays (as np.ndarray-s) False otherwise.

        Example:
            >>> import pyfaust as pf
            >>> F = pf.rand(10, 10, fac_type='sparse')
            >>> F
            Faust size 10x10, density 2.5, nnz_sum 250, 5 factor(s):
            - FACTOR 0 (double) SPARSE, size 10x10, density 0.5, nnz 50
            - FACTOR 1 (double) SPARSE, size 10x10, density 0.5, nnz 50
            - FACTOR 2 (double) SPARSE, size 10x10, density 0.5, nnz 50
            - FACTOR 3 (double) SPARSE, size 10x10, density 0.5, nnz 50
            - FACTOR 4 (double) SPARSE, size 10x10, density 0.5, nnz 50

            >>> F.isdense()
            False
            >>> F = pf.rand(10, 10, fac_type='dense')
            >>> F
            Faust size 10x10, density 2.5, nnz_sum 250, 5 factor(s):
            - FACTOR 0 (double) DENSE, size 10x10, density 0.5, nnz 50
            - FACTOR 1 (double) DENSE, size 10x10, density 0.5, nnz 50
            - FACTOR 2 (double) DENSE, size 10x10, density 0.5, nnz 50
            - FACTOR 3 (double) DENSE, size 10x10, density 0.5, nnz 50
            - FACTOR 4 (double) DENSE, size 10x10, density 0.5, nnz 50

            >>> F.isdense()
            True
            >>> F = pf.rand(10, 10, fac_type='mixed')
            >>> F
            Faust size 10x10, density 2.5, nnz_sum 250, 5 factor(s):
            - FACTOR 0 (double) DENSE, size 10x10, density 0.5, nnz 50
            - FACTOR 1 (double) DENSE, size 10x10, density 0.5, nnz 50
            - FACTOR 2 (double) SPARSE, size 10x10, density 0.5, nnz 50
            - FACTOR 3 (double) DENSE, size 10x10, density 0.5, nnz 50
            - FACTOR 4 (double) SPARSE, size 10x10, density 0.5, nnz 50

            >>> F.isdense()
            False

        <b>See also</b> Faust.issparse, pyfaust.rand, pyfaust.rand_bsr
        """
        return F.m_faust.is_all_dense()

    def swap_cols(F, id1, id2, permutation=False, inplace=False):
        """
        Swaps F columns of indices id1 and id2.

        Args:
            id1: index of the first column of the swap.
            id2: index of the second column of the swap.
            permutation: if True then the swap is performed by inserting a permutation
            matrix to the output Faust. If False, the last matrix
            in the Faust F sequence is edited to swap the columns.
            inplace: if True then F is modified instead of returning a new Faust.
            Otherwise, by default, a new Faust is returned.

        Returns:
            The column swapped Faust.

		Example:
			>>> from pyfaust import rand as frand
			>>> F = frand(10, 10)
			>>> G = F.swap_cols(2,4)
			>>> G
			Faust size 10x10, density 2.5, nnz_sum 250, 5 factor(s): 
			- FACTOR 0 (real) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 1 (real) DENSE,  size 10x10, density 0.5, nnz 50
			- FACTOR 2 (real) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 3 (real) DENSE,  size 10x10, density 0.5, nnz 50
			- FACTOR 4 (real) DENSE,  size 10x10, density 0.5, nnz 50

			>>> G[:, 2].toarray() == F[:, 4].toarray()
			array([[ True],
				   [ True],
				   [ True],
				   [ True],
				   [ True],
				   [ True],
				   [ True],
				   [ True],
				   [ True],
				   [ True]])
			>>> G[:, 4].toarray() == F[:, 2].toarray()
			array([[ True],
				   [ True],
				   [ True],
				   [ True],
				   [ True],
				   [ True],
				   [ True],
				   [ True],
				   [ True],
				   [ True]])
			>>> H  = F.swap_cols(2,4, permutation=True)
			>>> H
			Faust size 10x10, density 2.6, nnz_sum 260, 6 factor(s): 
			- FACTOR 0 (real) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 1 (real) DENSE,  size 10x10, density 0.5, nnz 50
			- FACTOR 2 (real) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 3 (real) DENSE,  size 10x10, density 0.5, nnz 50
			- FACTOR 4 (real) DENSE,  size 10x10, density 0.5, nnz 50
			- FACTOR 5 (real) SPARSE, size 10x10, density 0.1, nnz 10

        <b>See also</b> Faust.swap_rows
        """
        if(inplace):
            F.m_faust.swap_cols(id1, id2, permutation,
                                inplace)
            return F
        F_swapped = Faust(core_obj=F.m_faust.swap_cols(id1, id2, permutation,
                                                       inplace))
        return F_swapped

    def swap_rows(F, id1, id2, permutation=True, inplace=False):
        """
        Swaps F rows of indices id1 and id2.

        Args:
            id1: index of the first row of the swap.
            id2: index of the second row of the swap.
            permutation: if True then the swap is performed by inserting a permutation
            matrix to the output Faust. If False, the last matrix
            in the Faust F sequence is edited to swap the rows.
            inplace: if True then F is modified instead of returning a new Faust.
            Otherwise, by default, a new Faust is returned.

        Returns:
            The rows swapped Faust.

		Example:
			>>> from pyfaust import rand as frand
			>>> F = frand(10, 10)
			>>> G = F.swap_rows(2,4)
			>>> G
			Faust size 10x10, density 2.5, nnz_sum 250, 5 factor(s): 
			- FACTOR 0 (real) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 1 (real) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 2 (real) DENSE,  size 10x10, density 0.5, nnz 50
			- FACTOR 3 (real) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 4 (real) DENSE,  size 10x10, density 0.5, nnz 50

			>>> G[2,:].toarray() == F[4,:].toarray()
			array([[ True,  True,  True,  True,  True,  True,  True,  True,  True,
					 True]])
			>>> G[4,:].toarray() == F[2,:].toarray()
			array([[ True,  True,  True,  True,  True,  True,  True,  True,  True,
					 True]])
			>>> H  = F.swap_rows(2,4, permutation=True)
			>>> H
			Faust size 10x10, density 2.6, nnz_sum 260, 6 factor(s): 
			- FACTOR 0 (real) SPARSE, size 10x10, density 0.1, nnz 10
			- FACTOR 1 (real) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 2 (real) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 3 (real) DENSE,  size 10x10, density 0.5, nnz 50
			- FACTOR 4 (real) SPARSE, size 10x10, density 0.5, nnz 50
			- FACTOR 5 (real) DENSE,  size 10x10, density 0.5, nnz 50

        <b>See also</b> Faust.swap_cols
        """
        if(inplace):
            F.m_faust.swap_rows(id1, id2, permutation,
                                inplace)
            return F
        F_swapped = Faust(core_obj=F.m_faust.swap_rows(id1, id2, permutation,
                                                       inplace))
        return F_swapped

    def optimize_memory(F):
        """Optimizes a Faust by changing the storage format of each factor in order to optimize the memory size.

        Returns:
            The optimized Faust.

        Example:
            >>> from pyfaust import rand
            >>> F = rand(1024, 1024, fac_type='mixed')
            >>> F
            Faust size 1024x1024, density 0.0244141, nnz_sum 25600, 5 factor(s):
            - FACTOR 0 (double) DENSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 1 (double) DENSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 2 (double) DENSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 3 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 4 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120

            >>> F.nbytes
            25296904
            >>> pF = F.optimize_memory()
            >>> pF
            Faust size 1024x1024, density 0.0244141, nnz_sum 25600, 5 factor(s):
            - FACTOR 0 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 1 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 2 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 3 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120
            - FACTOR 4 (double) SPARSE, size 1024x1024, density 0.00488281, nnz 5120

            >>> pF.nbytes
            327700


        <b>See also</b> Faust.optimize
        """
        F_opt = Faust(core_obj=F.m_faust.optimize_storage(False))
        return F_opt

    def optimize(F, transp=False):
        """Returns a Faust optimized with Faust.pruneout, Faust.optimize_memory and Faust.optimize_time.

        Args:
            transp: True in order to optimize the Faust according to its transpose.

        Returns:
            The optimized Faust.

        Note: this function is still experimental, you might use manually
        Faust.optimize_time, Faust.optimize_memory or Faust.pruneout to be
        more specific about the optimization to proceed.

        Example:
			This example shows how Faust.optimize can diminish the memory size and
			speed up the Faust.toarray() and Faust-vector multiplication.
			>>> from pyfaust import rand
			>>> F = rand(1024, 1024, dim_sizes=[1, 1024], num_factors=32, fac_type='mixed')
			>>> F.nbytes
			24014520
			>>> pF = F.optimize()
			>>> pF.nbytes
			911816
			>>> from time import time
			>>> t1 = time(); F.toarray(); print("F.toarray() time:", time()-t1)
			array([[0.00000000e+00, 9.47768987e+07, 3.54841410e+07, ...,
					5.02648254e+08, 1.69126326e+08, 6.60395426e+08],
				   [0.00000000e+00, 4.53823851e+07, 1.69910070e+07, ...,
					2.40685014e+08, 8.09834069e+07, 3.16219690e+08],
				   [0.00000000e+00, 7.05413763e+07, 2.64104431e+07, ...,
					3.74115410e+08, 1.25878813e+08, 4.91524836e+08],
				   ...,
				   [0.00000000e+00, 6.19890617e+07, 2.32084866e+07, ...,
					3.28758359e+08, 1.10617488e+08, 4.31933316e+08],
				   [0.00000000e+00, 1.05375223e+08, 3.94521173e+07, ...,
					5.58856367e+08, 1.88038697e+08, 7.34243484e+08],
				   [0.00000000e+00, 1.04672252e+08, 3.91889277e+07, ...,
					5.55128177e+08, 1.86784271e+08, 7.29345296e+08]])
			F.toarray() time: 0.3663349151611328
			>>> t1 = time(); pF.toarray(); print("pF.toarray() time:", time()-t1)
			array([[0.00000000e+00, 9.47768987e+07, 3.54841410e+07, ...,
					5.02648254e+08, 1.69126326e+08, 6.60395426e+08],
				   [0.00000000e+00, 4.53823851e+07, 1.69910070e+07, ...,
					2.40685014e+08, 8.09834069e+07, 3.16219690e+08],
				   [0.00000000e+00, 7.05413763e+07, 2.64104431e+07, ...,
					3.74115410e+08, 1.25878813e+08, 4.91524836e+08],
				   ...,
				   [0.00000000e+00, 6.19890617e+07, 2.32084866e+07, ...,
					3.28758359e+08, 1.10617488e+08, 4.31933316e+08],
				   [0.00000000e+00, 1.05375223e+08, 3.94521173e+07, ...,
					5.58856367e+08, 1.88038697e+08, 7.34243484e+08],
				   [0.00000000e+00, 1.04672252e+08, 3.91889277e+07, ...,
					5.55128177e+08, 1.86784271e+08, 7.29345296e+08]])
			pF.toarray() time: 0.10646939277648926
			>>> import numpy as np
			>>> v = np.random.rand(F.shape[1])
			>>> t1 = time(); F@v;print("F@v time:", time()-t1)
			array([2.79905047e+11, 1.34028007e+11, 2.08330156e+11, ...,
				   1.83072593e+11, 3.11205142e+11, 3.09129061e+11])
			F@v time: 0.003871440887451172
			>>> t1 = time(); pF@v; print("pF@v time:", time()-t1)
			array([2.79905047e+11, 1.34028007e+11, 2.08330156e+11, ...,
				   1.83072593e+11, 3.11205142e+11, 3.09129061e+11])
			pF@v time: 0.0013530254364013672

        <b>See also</b> Faust.optimize_time, Faust.optimize_memory, Faust.pruneout, Faust.nbytes, <a href="https://faust.inria.fr/tutorials/pyfaust-jupyter-notebooks/faust-optimizations-with-pyfaust/">Jupyter notebook about Faust optimizations</a>
        """
        F_opt = Faust(core_obj=F.m_faust.optimize(transp))
        return F_opt

    def optimize_time(F, transp=False, inplace=False, nsamples=1, mat=None):
        """Returns a Faust configured with the quickest Faust-matrix multiplication mode (benchmark ran on the fly).

        NOTE: this function launches a small benchmark on the fly. Basically, the methods
        available differ by the order used to compute the matrix chain.

        The evaluated methods in the benchmark are listed in pyfaust.FaustMulMode.
        Although depending on the package you installed and the capability of your
        hardware the methods based on Torch library can be used.

        Args:
            inplace: to optimize the current Faust directly instead of returning a new
            Faust with the optimization enabled. If True, F is returned
            otherwise a new Faust object is returned.
            transp: True in order to optimize the Faust according to its transpose.
            nsamples: the number of Faust-Dense matrix products
            calculated in order to measure time taken by each method (it could matter
            to discriminate methods when the performance is similar). By default,
            only one product is computed to evaluate the method.
            mat: if not None must be a numpy.ndarray or a
            scipy.sparse.csr_matrix. Use this argument to run the benchmark on
            the Faust multiplication by the matrix mat instead of Faust.toarray() (if mat
            is None). Note that mat must be of the same dtype as F.

        Returns:
            The optimized Faust.

        Example:
			>>> from pyfaust import rand
			>>> from time import time
			>>> F = rand(1024, 1024, dim_sizes=[1, 1024], num_factors=32, fac_type='dense')
			>>> oF = F.optimize_time()
			best method measured in time on operation Faust-toarray is: DYNPROG
			>>> t1 = time(); F.toarray(); print("F.toarray() time:", time()-t1)
			array([[1.08456293e+08, 9.06809253e+08, 2.14092678e+08, ...,
					2.39420597e+08, 1.73687901e+08, 1.48511336e+08],
				   [1.45686856e+08, 1.21809603e+09, 2.87585796e+08, ...,
					3.21608188e+08, 2.33310996e+08, 1.99491872e+08],
				   [1.37868515e+08, 1.15272644e+09, 2.72152397e+08, ...,
					3.04348964e+08, 2.20790272e+08, 1.88786066e+08],
				   ...,
				   [1.60235010e+08, 1.33973397e+09, 3.16303847e+08, ...,
					3.53723686e+08, 2.56609207e+08, 2.19412950e+08],
				   [1.24235524e+08, 1.03874022e+09, 2.45240878e+08, ...,
					2.74253715e+08, 1.98957633e+08, 1.70118148e+08],
				   [1.19470236e+08, 9.98897401e+08, 2.35834204e+08, ...,
					2.63734204e+08, 1.91326235e+08, 1.63592947e+08]])
			F.toarray() time: 0.2891380786895752
			>>> t1 = time(); oF.toarray(); print("oF.toarray() time:", time()-t1)
			array([[1.08456293e+08, 9.06809253e+08, 2.14092678e+08, ...,
					2.39420597e+08, 1.73687901e+08, 1.48511336e+08],
				   [1.45686856e+08, 1.21809603e+09, 2.87585796e+08, ...,
					3.21608188e+08, 2.33310996e+08, 1.99491872e+08],
				   [1.37868515e+08, 1.15272644e+09, 2.72152397e+08, ...,
					3.04348964e+08, 2.20790272e+08, 1.88786066e+08],
				   ...,
				   [1.60235010e+08, 1.33973397e+09, 3.16303847e+08, ...,
					3.53723686e+08, 2.56609207e+08, 2.19412950e+08],
				   [1.24235524e+08, 1.03874022e+09, 2.45240878e+08, ...,
					2.74253715e+08, 1.98957633e+08, 1.70118148e+08],
				   [1.19470236e+08, 9.98897401e+08, 2.35834204e+08, ...,
					2.63734204e+08, 1.91326235e+08, 1.63592947e+08]])
			oF.toarray() time: 0.0172119140625


        <b>See also</b> Faust.optimize, pyfaust.FaustMulMode

        """
        if(inplace):
            F.m_faust.optimize_time(transp, inplace, nsamples, M=mat)
            return F
        else:
            F_opt = Faust(core_obj=F.m_faust.optimize_time(transp, inplace,
                                                          nsamples, M=mat))
            return F_opt

    def copy(F, dev='cpu'):
        """Clone alias function (here just to mimic numpy API).
        <b>See also</b> Faust.clone
        """
        check_dev(dev)
        return F.clone(dev)

    def clone(F, dev=None):
        """Clones the Faust (in a new memory space).

        Args:
            dev (optional): 'cpu' to clone on CPU RAM, 'gpu' to clone on
            the GPU device. By default (None), the device is
            the F.device.

        Returns:
            The Faust clone.

        Example:
            >>> from pyfaust import rand
            >>> F = rand(10, 10)
            >>> F.clone()
            Faust size 10x10, density 2.5, nnz_sum 250, 5 factor(s):
                - FACTOR 0 (double) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 1 (double) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 2 (double) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 3 (double) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 4 (double) SPARSE, size 10x10, density 0.5, nnz 50

            >>> F.clone(dev='gpu') # only if a NVIDIA compatible GPU is available
			- GPU FACTOR 0 (double) SPARSE size 10 x 10, addr: 0x2c85390, density 0.500000, nnz 50
			- GPU FACTOR 1 (double) SPARSE size 10 x 10, addr: 0x7ff00c0, density 0.500000, nnz 50
			- GPU FACTOR 2 (double) SPARSE size 10 x 10, addr: 0x977f280, density 0.500000, nnz 50
			- GPU FACTOR 3 (double) SPARSE size 10 x 10, addr: 0x9780120, density 0.500000, nnz 50
			- GPU FACTOR 4 (double) SPARSE size 10 x 10, addr: 0x9780fc0, density 0.500000, nnz 50

        """
        if dev == None:
            dev = F.device
        check_dev(dev)
        # dev is 'gpu[:id]' or 'cpu'
        if F.device.startswith('gpu'):
            if F.dtype == 'double':
                clone_F = \
                        Faust(core_obj=_FaustCorePy.FaustCoreGenNonMemberFuncsDblGPU.clone(F.m_faust,
                                                                                           dev))
            else: # F.dtype == np.complex
                clone_F = \
                        Faust(core_obj=_FaustCorePy.FaustCoreGenNonMemberFuncsCplxDblGPU.clone(F.m_faust,
                                                                                            dev))
#            clone_F = Faust(core_obj=F.m_faust.clone(dev))
        elif F.device == 'cpu':
            if dev == 'cpu':
                clone_F = Faust(core_obj=F.m_faust.clone(-1))
            else:
                if F.dtype == 'double':
                    clone_F = \
                            Faust(core_obj=_FaustCorePy.FaustCoreGenNonMemberFuncsDblCPU.clone(F.m_faust, dev))
                else:
                    clone_F = \
                            Faust(core_obj=_FaustCorePy.FaustCoreGenNonMemberFuncsCplxDblCPU.clone(F.m_faust,
                                                                                                   dev))
        else:
            raise ValueError("F.device is not valid")
        return clone_F

    def sum(F, axis=None, **kwargs):
        """
        Sums Faust elements over a given axis.

        Args:
            axis (optional): None or int  or tuple of ints.
            Axis or axes along which the the sum is performed

        Returns:
            The Faust sum.

        Example:
            >>> from pyfaust import rand as frand
            >>> F = frand(10, 10)
            >>> F.sum()
            Faust size 1x10, density 26, nnz_sum 260, 6 factor(s): 
                - FACTOR 0 (real) DENSE,  size 1x10, density 1, nnz 10
                - FACTOR 1 (real) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 2 (real) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 3 (real) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 4 (real) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 5 (real) SPARSE, size 10x10, density 0.5, nnz 50

            >>> F.sum(axis=0).toarray()
            array([[135.62885806,  86.91727358, 212.40112068, 186.06476227,
                             68.74097449,  88.63216727,  64.38497784,  58.51572934,
                             91.05523782, 196.79601819]])
            >>> F.sum(axis=1)
            Faust size 10x1, density 26, nnz_sum 260, 6 factor(s): 
                - FACTOR 0 (real) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 1 (real) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 2 (real) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 3 (real) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 4 (real) SPARSE, size 10x10, density 0.5, nnz 50
                - FACTOR 5 (real) DENSE,  size 10x1, density 1, nnz 10
            >>> F.sum(axis=1).toarray()
            array([[ 78.7358253 ],
                  [122.06237092],
                  [171.60995629],
                  [110.86003948],
                  [147.82414116],
                  [100.35211187],
                  [123.56100581],
                  [104.49754233],
                  [ 99.99809178],
                  [129.63603461]])

        """
        if axis == None:
            axis = (0,1)
        is_tuple = isinstance(axis, tuple)
        is_int = isinstance(axis, int)
        is_tuple_or_int = is_tuple or is_int
        if not is_tuple_or_int or is_tuple and \
           (not isinstance(axis[0], int) or not isinstance(axis[1], int)):
            raise TypeError("axis must be int or tuple of ints")
        if axis == None or axis == 0 or is_tuple and 0 in axis:
            F = Faust(np.ones((1, F.shape[0])), dev=F.device)@F
        if axis == 1 or is_tuple and 1 in axis:
            F = F@Faust(np.ones((F.shape[1], 1)), dev=F.device)
        if is_tuple and len([i for i in axis if i < 0
                             or i > 1]) or is_int and (axis < 0 or axis > 1):
            raise ValueError("axis "+str(axis)+" is out of bounds for a Faust "
                             " (only two dimensions)")
        return F

    def sliceMultiply(F, start_col_id, end_col_id, vec):
        return F.m_faust.colSliceMultiply(start_col_id, end_col_id, vec)

    def indexMultiply(F, x, d0_ids=None, d1_ids=None):
        if d0_ids is not None:
            d0_ids = d0_ids if isinstance(d0_ids, np.ndarray) else np.array(d0_ids)
        if d1_ids is not None:
            d1_ids = d1_ids if isinstance(d1_ids, np.ndarray) else np.array(d1_ids)
        return F.m_faust.indexMultiply(d0_ids, d1_ids, x)

    def average(F, axis=None, weights=None, returned=False):
        """
        Computes the weighted average of F along the specified axis.

        Args:
            axis (optional): None or int  or tuple of ints.
            Axis or axes along which to average the Faust F.
			The default, axis=None, will average over all of the elements of the input array.
			If axis is a tuple of ints, averaging is performed on all of the axes specified in the tuple
			weights: an array of weights associated with the values in F. Each value in F contributes to the average according to its associated weight. The weights array can either be 1-D (in which case its length must be the size of a along the given axis) or of the same shape as a. If weights=None, then all data in F are assumed to have a weight equal to one. The 1-D calculation is:

			avg = sum(F @ weights) / sum(weights)

			The only constraint on weights is that sum(weights) must not be 0.

        Returns:
            The Faust average.

        Example:
			>>> from pyfaust import Faust
			>>> import numpy as np
			>>> data = np.arange(1, 5)
			>>> data
			array([1, 2, 3, 4])
			>>> F = Faust(data.reshape(1,data.shape[0]))
			>>> FA = F.average()
			>>> FA
			Faust size 1x1, density 9, nnz_sum 9, 3 factor(s): 
			- FACTOR 0 (real) DENSE,  size 1x1, density 1, nnz 1
			- FACTOR 1 (real) DENSE,  size 1x4, density 1, nnz 4
			- FACTOR 2 (real) DENSE,  size 4x1, density 1, nnz 4

			>>> FA.toarray()
			array([[2.5]])
			>>> 
			>>> data2 = np.arange(6).reshape((3,2))
			>>> F2 = Faust(data2)
			>>> F2.average(axis=1, weights=[1./4, 3./4]).toarray()
			array([[0.75],
				   [2.75],
				   [4.75]])

        """
        if axis == None:
            axis = (0,1)
        is_tuple = isinstance(axis, tuple)
        is_int = isinstance(axis, int)
        is_tuple_or_int = is_tuple or is_int
        if not is_tuple_or_int or is_tuple and \
           (not isinstance(axis[0], int) or not isinstance(axis[1], int)):
            raise TypeError("axis must be int or tuple of ints")
        if not isinstance(weights, np.ndarray) and weights == None:
            def_rweights = np.ones(F.shape[0])
            def_cweights = np.ones(F.shape[1])
            if isinstance(axis, int):
                if axis == 0:
                    weights = def_rweights
                elif axis == 1:
                    weights = def_cweights
            elif isinstance(axis, tuple):
                if 0 in axis:
                    F = F.average(axis=0,
                                 weights=def_rweights,
                                 returned=returned)
                if 1 in axis:
                    F = F.average(axis=1,
                                  weights=def_cweights,
                                  returned=returned)
                return F

        if not isinstance(weights, np.ndarray):
            weights = np.array(weights)
        if weights.shape != F.shape and weights.ndim != 1:
            raise TypeError("1D weights expected when shapes of a and weights"
                            " differ.")

        if weights.shape == F.shape:
            aF = Faust([(weights[:,0].T).reshape(1,F.shape[0])], dev=F.device)@F[:,0]
            for i in range(1,F.shape[1]):
                aF = pyfaust.hstack((aF, Faust([(weights[:,i].T).reshape(1,
                                                                         F.shape[0])], dev=F.device)@F[:,i]))
            sum_weights = 1/np.sum(weights, axis=axis)
            aFw =  aF[:,0]@Faust([sum_weights[0].reshape(1,1)], dev=F.device)
            for i in range(1, sum_weights.shape[0]):
                aFw = pyfaust.hstack((aFw,
                                      aF[:,i]@Faust([sum_weights[i].reshape(1,1)], dev=F.device)))
            if(returned):
                return (aFw, sum_weights)
            return aFw

        if axis == 1 or isinstance(axis, tuple) and 1 in axis:
            if weights.shape[0] == F.shape[1]:
                aF = F@Faust(weights.reshape(weights.size, 1), dev=F.device)
            else:
                raise ValueError("ValueError: Length of weights not compatible"
                                 " with specified axis 1.")
            sum_weights = np.sum(weights.reshape(weights.size, 1), axis=0)[0]

        if axis == 0 or isinstance(axis, tuple) and 0 in axis:
            weightsM = weights
            if weights.ndim == 1:
                weightsM = weights.reshape(1, weights.size)
            if weightsM.shape[1] == F.shape[0]:
                aF = Faust(weightsM, dev=F.device)@F
            else:
                raise ValueError("ValueError: Length of weights not compatible"
                                 " with axis 0.")
            sum_weights = np.sum(weightsM.reshape(1, weights.size),
                                 axis=1)[0]

        if sum_weights != 0:
            aF = aF * (1/sum_weights)
        else:
            raise decimal.DivisionByZero("Weights sum to zero, can't be normalized")
        if(returned):
            return (aF, sum_weights)
        return aF

pyfaust.Faust.__div__ = pyfaust.Faust.__truediv__

def implements(numpy_function):
    """Register an __array_function__ implementation for MyArray
    objects."""
    def decorator(func):
        HANDLED_FUNCTIONS[numpy_function] = func
        return func
    return decorator

def version():
    """Returns the FAuST package version.
    """
    return __version__

__version__ =  "@CPACK_PACKAGE_VERSION@"

def faust_fact(*args, **kwargs):
    """
    This function is a shorthand for pyfaust.fact.hierarchical.

    <b>See also</b> pyfaust.fact.hierarchical
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

    <b>See also</b> Faust.norm
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

@implements(np.dot)
def dot(A, B, **kwargs):
    """Returns Faust.dot(A,B) if A or B is a Faust object, returns numpy.dot(A,B) ortherwise.
    <b>See also</b> Faust.norm
    """
    if(Faust.isFaust(A)):
        return A.dot(B)
    elif(Faust.isFaust(B)):
        return B.T.dot(A.T).T
    else: # if F is not a Faust, try to rely on numpy (not breaking possible
          # past import)
        return np.dot(A,B)

def pinv(F):
    """A package function alias for the member function Faust.pinv().
    """
    if(Faust.isFaust(F)):
        return F.pinv()
    else:
        return np.linalg.pinv(F)


@implements(np.concatenate)
def concatenate(F, *args, axis=0, **kwargs):
    """
    A package function alias for the member function Faust.concatenate.

    Example:
        >>> from pyfaust import *
        >>> F1 = rand(5,50)
        >>> F2 = rand(5,50)
        >>> concatenate((F1, F2), axis=0)

    <b>See also</b> numpy.concatenate
    """
    if not isinstance(F, tuple):
        raise TypeError("first arg must be a tuple")
    _tuple = F
    if(Faust.isFaust(_tuple[0])):
        return _tuple[0].concatenate(*_tuple[1:], axis=axis)
    else:
        return np.concatenate(_tuple, axis=axis)

def hstack(_tuple):
    """
    Concatenates horizontally Faust-s and/or numpy.ndarray objects using Faust.concatenate().

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
    Concatenates vertically Faust-s and/or numpy.ndarray arrays using Faust.concatenate().

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
    Package alias function of Faust.isFaust.

    <b>See also:</b> Faust.__init__, Faust.isFaust.
    """
    return Faust.isFaust(obj)

def wht(n, normed=True, dev="cpu", dtype='double'):
    """
    Constructs a Faust implementing the Walsh-Hadamard transform of order n.

    The resulting Faust has log2(n) sparse factors of order n, each one having
    2 nonzero per row and per column.

    Args:
       n: the power of two exponent for a Hadamard matrix of order n
       and a factorization in log2(n) factors.
       normed: default to True to normalize the Hadamard Faust as if you called
       Faust.normalize() and False otherwise.
       dev: device on which to create the Faust ('cpu' or 'gpu').
       dtype: the Faust dtype, it must be 'double' or 'complex'.

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
    if dtype not in ['double', 'complex', 'float']:
        raise ValueError('dtype argument must be double or complex')
    check_dev(dev)
    log2n = np.floor(np.log2(n))
    if(n > 2**log2n): raise ValueError("n must be a power of 2.")
    if(not isinstance(normed, bool)):
        raise TypeError("normed must be True of False.")
    if dev == "cpu":
        if dtype == 'double':
            H = Faust(core_obj=_FaustCorePy.FaustAlgoGenCPUDbl.hadamardFaust(log2n, normed))
        elif dtype == 'float':
            H = Faust(core_obj=_FaustCorePy.FaustAlgoGenCPUFlt.hadamardFaust(log2n, normed))
        else: # dtype == 'complex'
            H = Faust(core_obj=_FaustCorePy.FaustAlgoGenCPUCplxDbl.hadamardFaust(log2n, normed))
    elif dev.startswith("gpu"):
        if dtype == 'double':
            H = Faust(core_obj=_FaustCorePy.FaustAlgoGenGPUDbl.hadamardFaust(log2n, normed))
        elif dtype == 'float':
            raise TypeError("float is not yet supported by GPU wht")
        else: # dtype == 'complex'
            H = Faust(core_obj=_FaustCorePy.FaustAlgoGenGPUCplxDbl.hadamardFaust(log2n, normed))
    return H

def dft(n, normed=True, dev='cpu'):
    """
    Constructs a Faust F such that F.toarray() is the Discrete Fourier Transform square matrix of order n.

    The factorization algorithm used is Cooley-Tukey.

    The resulting Faust is complex and has (log2(n)+1) sparse factors
    whose the log2(n) first has 2 nonzeros per row and per column.
    The last factor is a permutation matrix.

    Args:
        n: the power of two exponent for a DFT of order n and a
        factorization in log2(n)+1 factors.
        normed: default to True to normalize the DFT Faust as if you called
        Faust.normalize() and False otherwise.
        dev: device to create the Faust on ('cpu' or 'gpu').

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
    if dev == "cpu":
        F = Faust(core_obj=_FaustCorePy.FaustAlgoCplxDblGenCPU.fourierFaust(log2n, normed))
    elif dev.startswith("gpu"):
        F = Faust(core_obj=_FaustCorePy.FaustAlgoCplxDblGenGPU.fourierFaust(log2n, normed))
    return F

def dct(n, dev='cpu'):
    """Returns the Direct Cosine Transform (Type II) Faust of order n.

    The analytical formula of DCT II used here is:
        \f$2 \sum_{n=0}^{N-1} x_n cos \left( {\pi k (2n + 1)} \over {2N} \right)\f$


    Args:
        n: the order of the DCT (must be a power of two).
        dev: the device on which the Faust is created.

    Example:
        >>> from pyfaust import dct
        >>> from scipy.fft import dct as scipy_dct
        >>> import numpy as np
        >>> DCT8 = dct(8)
        >>> x = np.arange(8)
        >>> np.real(DCT8@x)
        array([ 56.        , -25.76929209,   0.        ,  -2.6938192 ,
               0.        ,  -0.80361161,   0.        ,  -0.20280929])
        >>> scipy_dct(x)
        array([ 56.        , -25.76929209,   0.        ,  -2.6938192 ,
               0.        ,  -0.80361161,   0.        ,  -0.20280929])
        >>> np.allclose(DCT8@x, scipy_dct(x))
        True

    See also pyfaust.dft, pyfaust.dst, <a
    href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.dct.html">
    scipy.fft.dct</a>
    """
    DFT = pyfaust.dft(n, dev='cpu', normed=False)
#    P_ = np.zeros((n, n))
#    P_[np.arange(0, n//2), np.arange(0, n, 2)] = 1
#    P_[np.arange(n//2, n), np.arange(1, n, 2)[::-1]] = 1
    # P_ as sparse matrix
    P_row_inds = np.arange(0, n)
    P_col_inds = np.hstack((np.arange(0, n, 2), np.arange(1, n, 2)[::-1]))
    P_ = csr_matrix((np.ones(n), (P_row_inds, P_col_inds)), shape=(n, n))
    E = diags([2*np.exp(-1j*np.pi*k/2/n) for k in range(0, n)])
    f0 = csr_matrix(E @ DFT.factors(0))
    f_end = csr_matrix(DFT.factors(len(DFT)-1) @ P_)
    mid_factors = DFT.factors(range(1, len(DFT)-1))
    if pyfaust.isFaust(mid_factors):
        mid_F = mid_factors
    else:
        mid_F = Faust(mid_factors)
    DCT = (Faust(f0) @ mid_F @ Faust(f_end)).real
    return DCT

# experimental block start
def dst2(n, dev='cpu'):
    """
    See issue #265 paper
    Learning_Fast_Algorithms_for_Linear_Transforms_Using_Butterfly.pdf
    This implementation doesn't work! wasn't able to fix it.
    """
    DFT = pyfaust.dft(n, dev='cpu', normed=False)
    P_ = np.zeros((n, n))
    P_[np.arange(0, n//2), np.arange(0, n, 2)] = 1
    P_[np.arange(n//2, n), np.arange(1, n, 2)[::-1]] = 1
    I_N2 = seye(n//2, format='csr')
    I_2 = np.array([ [1, 0], [0, -1]])
    D = kron(I_2, I_N2, format='csr')
    D = D @ diags([np.exp(-1j*2*np.pi*k/n) for k in range(n)])
    E = diags([2*1j*np.exp(-1j*np.pi*k/2/n) for k in range(n)])
#    E = diags([2*1j*np.exp(-1j*5*np.pi*k/2/n) for k in range(n)])
    f0 = csr_matrix(E @ DFT.factors(0))
    f_end = csr_matrix(DFT.factors(len(DFT)-1) @ D @ P_)
    mid_factors = DFT.factors(range(1, len(DFT)-1))
    if pyfaust.isFaust(mid_factors):
        mid_F = mid_factors
    else:
        mid_F = Faust(mid_factors)
    DST = Faust(f0) @ mid_F @ Faust(f_end)
    return DST

def dst3(n, dev='cpu'):
    """
    This one works (cf. DST.pdf issue #265) except for the last
    frequency of the DST. dst function below fixes it rewriting the FFT
    with shifting.
    """
    DFT = pyfaust.dft(n, normed=False)
    D1 = csr_matrix(-2*diags([- 1j * np.exp(-1j * np.pi / 2 / n * k) for k in range(0,
                                                                            n)]))
    D2 = csr_matrix(diags([np.exp(-1j * np.pi / n * k) for k in range(0, n)]))
    S = csr_matrix(diags(np.ones(n-1), 1)).astype('complex')
#    S[-1, 0] = 1
    P_ = np.zeros((n*2, n))
    P_[np.arange(0, n//2), np.arange(0, n, 2)] = 1
    P_[np.arange(n, n + n // 2), np.arange(1, n, 2)] = 1
    SD1 = Faust([S, D1])
    F_even = SD1 @ DFT
    F_odd = SD1 @ Faust(D2) @ DFT
    F = pyfaust.hstack((F_even, F_odd))
    F = F @ Faust(P_)
    return F.real
# experimental block end

def dst(n, dev='cpu'):
    """
    Returns the Direct Sine Transform (Type II) Faust of order n.

    The analytical formula of DST II used here is:
        \f$2 \sum_{n=0}^{N-1} x_n sin \left( {\pi (k+1) (2n + 1)} \over {2N} \right)\f$

    Args:
        n: the order of the DST (must be a power of two).
        dev: the device on which the Faust is created.

    Example:
        >>> from pyfaust import dst
        >>> from scipy.fft import dst as scipy_dst
        >>> import numpy as np
        >>> DST8 = dst(8)
        >>> x = np.ones(8)
        >>> np.real(DST8@x)
        array([1.02516618e+01, 4.93038066e-32, 3.59990489e+00, 1.97215226e-31,
               2.40537955e+00, 9.86076132e-32, 2.03918232e+00,
               0.00000000e+00])
        >>> scipy_dst(x)
        array([10.25166179,  0.        ,  3.59990489,  0.        ,  2.40537955,
                       0.        ,  2.03918232,  0.        ])
        >>> np.allclose(DST8@x, scipy_dst(x))
        True

    See also pyfaust.dft, pyfaust.dct, <a
    href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.dst.html">
    scipy.fft.dst</a>
    """
    def bitrev(inds):
        """
        Bitreversal indices.
        """
        n = len(inds)
        if n == 1:
            return inds
        else:
            even = bitrev(inds[np.arange(0, n, 2, dtype='int')])
            odd = bitrev(inds[np.arange(1, n, 2, dtype='int')])
            return np.hstack((even, odd))

    def bitrev_perm(N):
        """
        Bitreversal permutation.
        """
        row_inds = np.arange(0, N, dtype='int')
        col_inds = bitrev(row_inds)
        ones = np.ones((N), dtype='float')
        return csr_matrix((ones, (row_inds, col_inds)), shape=(N, N))

    def omega(N):
        """
        Returns the list of n-th root of unity raised to the power of -(k+1) (instead of
        k as in the FFT, the purpose here is to write the DST).
        """
        omega = np.exp(np.pi*1j/N) # exp(-i*2*pi/2N)
        return diags([omega**-(k+1) for k in range(0, N)]).tocsr()

    def B(N):
        """
        Butterfly factor of order N.
        """
        I_N2 = seye(N//2, format='csr')
        O_N2 = omega(N//2)
        return svstack((shstack((I_N2, O_N2)), shstack((I_N2, - O_N2))))

    def mod_fft(N, dev='cpu'):
        """
        Modified FFT for DST computation (it would have been the standard FFT
        if omega was raised to the power of -k instead of -(k+1)).
        """
        N_ = N
        Bs = []
        while N_ != 1:
            B_ = B(N_)
            diag_B = kron(seye(N//N_, format='csr'), B_).tocsr() # to avoid BSR
            # matrices because FAµST doesn't support BSR matrices concatenation
            # TODO: let BSR format after concatenation is supported
            Bs += [diag_B]
            N_ //= 2
        return Faust(Bs+[bitrev_perm(N).astype(Bs[-1].dtype)], dev=dev)

	# compute the DST (look at issue #265 for doc, S permutation was replaced by mod_fft to fix missing last frequency)
    MDFT = mod_fft(n, dev=dev)
    D1 = csr_matrix(-2*diags([- 1j * np.exp(-1j * np.pi / 2 / n * (k+1)) for k in range(0,
                                                                            n)]))
    D2 = csr_matrix(diags([np.exp(-1j * np.pi / n * (k+1)) for k in range(0, n)]))
    #    P_ = np.zeros((n*2, n))
    #    P_[np.arange(0, n//2), np.arange(0, n, 2)] = 1
    #    P_[np.arange(n, n + n // 2), np.arange(1, n, 2)] = 1
    # P_ as as sparse matrix
    P_row_inds = np.hstack((np.arange(0, n//2), np.arange(n, n + n // 2)))
    P_col_inds = np.hstack((np.arange(0, n, 2), np.arange(1, n, 2)))
    P_ = csr_matrix((np.ones(n), (P_row_inds, P_col_inds)), shape=(n*2, n))
    F_even = Faust(D1, dev=dev) @ MDFT
    F_odd = Faust(D1, dev=dev) @ Faust(D2) @ MDFT
    F = pyfaust.hstack((F_even, F_odd))
    F = F @ Faust(P_, dev=dev)
    return F.real

def circ(c, **kwargs):
    """Returns a circulant Faust C defined by the vector c (which is the first column of C.toarray()).

    Example:
        >>> from pyfaust import circ
        >>> import numpy as np
        >>> c = np.random.rand(8)
        >>> C = circ(c)
        >>> C
            Faust size 8x8, density 1.75, nnz_sum 112, 8 factor(s):
            - FACTOR 0 (complex) SPARSE, size 8x8, density 0.25, nnz 16
            - FACTOR 1 (complex) SPARSE, size 8x8, density 0.25, nnz 16
            - FACTOR 2 (complex) SPARSE, size 8x8, density 0.25, nnz 16
            - FACTOR 3 (complex) SPARSE, size 8x8, density 0.125, nnz 8
            - FACTOR 4 (complex) SPARSE, size 8x8, density 0.125, nnz 8
            - FACTOR 5 (complex) SPARSE, size 8x8, density 0.25, nnz 16
            - FACTOR 6 (complex) SPARSE, size 8x8, density 0.25, nnz 16
            - FACTOR 7 (complex) SPARSE, size 8x8, density 0.25, nnz 16

        >>> np.allclose(C.toarray()[:, 0], c)
        True
        >>> c
        array([0.98193782, 0.08954973, 0.48411884, 0.71426712, 0.25681982,
                      0.32582343, 0.91189492, 0.77211408])
        >>> np.real(C.toarray())
		array([[0.98193782, 0.77211408, 0.91189492, 0.32582343, 0.25681982, 0.71426712, 0.48411884, 0.08954973],
       [0.08954973, 0.98193782, 0.77211408, 0.91189492, 0.32582343, 0.25681982, 0.71426712, 0.48411884],
       [0.48411884, 0.08954973, 0.98193782, 0.77211408, 0.91189492, 0.32582343, 0.25681982, 0.71426712],
       [0.71426712, 0.48411884, 0.08954973, 0.98193782, 0.77211408, 0.91189492, 0.32582343, 0.25681982],
       [0.25681982, 0.71426712, 0.48411884, 0.08954973, 0.98193782, 0.77211408, 0.91189492, 0.32582343],
       [0.32582343, 0.25681982, 0.71426712, 0.48411884, 0.08954973, 0.98193782, 0.77211408, 0.91189492],
       [0.91189492, 0.32582343, 0.25681982, 0.71426712, 0.48411884, 0.08954973, 0.98193782, 0.77211408],
       [0.77211408, 0.91189492, 0.32582343, 0.25681982, 0.71426712, 0.48411884, 0.08954973, 0.98193782]])

    See also <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.circulant.html">scipy.linalg.circulant</a>, pyfaust.anticirc, pyfaust.toeplitz.
    """
    if isinstance(c, list):
        c = np.array(c)
    # hidden option
    if 'diag_factor' in kwargs.keys():
        diag_factor = kwargs['diag_factor']
    else:
        diag_factor = 'multiplied'
    if diag_factor not in ['multiplied', 'csr']:
        raise ValueError('option diag_factor must be \'csr\' or'
                         ' \'multiplied\'')
    log2c = np.log2(len(c))
    if  log2c != np.floor(log2c):
        raise ValueError('Only power of two length vectors are supported')
    if not isinstance(c, np.ndarray) or c.ndim != 1:
        raise TypeError('c must be a vector of numpy array type')
    n = len(c)
    F = dft(n, normed=False)
    FH = F.H
    S = csr_matrix(diags(FH@(c/n)))
#    S = csr_matrix(diags(np.sqrt(n)*FH@c)) # if it was normed==True
    if diag_factor == 'csr':
        C = F @ Faust(S) @ FH
    elif diag_factor == 'multiplied':
        C = F @ Faust(S@FH.factors(0)) @ FH.right(1)
    return C

def anticirc(c):
    """Returns an anti-circulant Faust A defined by the vector c (which is the last column of A.toarray()).

    Example:
        >>> from pyfaust import anticirc
        >>> import numpy as np
        >>> c = np.random.rand(8)
        >>> A = anticirc(c)
        >>> A
        Faust size 8x8, density 1.75, nnz_sum 112, 8 factor(s):
        - FACTOR 0 (complex) SPARSE, size 8x8, density 0.25, nnz 16
        - FACTOR 1 (complex) SPARSE, size 8x8, density 0.25, nnz 16
        - FACTOR 2 (complex) SPARSE, size 8x8, density 0.25, nnz 16
        - FACTOR 3 (complex) SPARSE, size 8x8, density 0.125, nnz 8
        - FACTOR 4 (complex) SPARSE, size 8x8, density 0.125, nnz 8
        - FACTOR 5 (complex) SPARSE, size 8x8, density 0.25, nnz 16
        - FACTOR 6 (complex) SPARSE, size 8x8, density 0.25, nnz 16
        - FACTOR 7 (complex) SPARSE, size 8x8, density 0.25, nnz 16

        >>> np.allclose(A.toarray()[:, -1], c)
        True
        >>> c
        array([0.62656024, 0.19827102, 0.54318481, 0.59119001, 0.28540344,
               0.00098907, 0.14192036, 0.35979329])
        >>> np.real(A.toarray())
        array([[0.19827102, 0.54318481, 0.59119001, 0.28540344, 0.00098907, 0.14192036, 0.35979329, 0.62656024],
               [0.54318481, 0.59119001, 0.28540344, 0.00098907, 0.14192036, 0.35979329, 0.62656024, 0.19827102],
               [0.59119001, 0.28540344, 0.00098907, 0.14192036, 0.35979329, 0.62656024, 0.19827102, 0.54318481],
               [0.28540344, 0.00098907, 0.14192036, 0.35979329, 0.62656024, 0.19827102, 0.54318481, 0.59119001],
               [0.00098907, 0.14192036, 0.35979329, 0.62656024, 0.19827102, 0.54318481, 0.59119001, 0.28540344],
               [0.14192036, 0.35979329, 0.62656024, 0.19827102, 0.54318481, 0.59119001, 0.28540344, 0.00098907],
               [0.35979329, 0.62656024, 0.19827102, 0.54318481, 0.59119001, 0.28540344, 0.00098907, 0.14192036],
               [0.62656024, 0.19827102, 0.54318481, 0.59119001, 0.28540344, 0.00098907, 0.14192036, 0.35979329]])

    See also pyfaust.circ, pyfaust.toeplitz.
    """
    G = circ(c)
    I = np.arange(len(c)-1, -1, -1)
    J = np.arange(0, len(c))
    P = csr_matrix((np.ones(J.size), (I, J)))
    return G.left(len(G)-2)@Faust(G.factors(len(G)-1)@P)

def toeplitz(c, r=None):
    """Constructs a toeplitz Faust whose the first column is c and the first row r.

    Args:
        c: the first column of the toeplitz Faust.
        r: the first row of the toeplitz Faust. If none then r =
        np.conjugate(c). r[0] is ignored, the first row is always [r[0],
        r[1:]].

    Example:
        >>> from pyfaust import toeplitz
        >>> import numpy as np
        >>> c = np.random.rand(10)
        >>> T = toeplitz(c)
        >>> T
            Faust size 10x10, density 6.16, nnz_sum 616, 12 factor(s):
            - FACTOR 0 (complex) SPARSE, size 10x32, density 0.0625, nnz 20
            - FACTOR 1 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 2 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 3 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 4 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 5 (complex) SPARSE, size 32x32, density 0.03125, nnz 32
            - FACTOR 6 (complex) SPARSE, size 32x32, density 0.03125, nnz 32
            - FACTOR 7 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 8 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 9 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 10 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 11 (complex) SPARSE, size 32x10, density 0.0625, nnz 20

        >>> np.allclose(T.toarray()[:, 0], c)
        True
        >>> np.allclose(T.toarray()[0, :], c)
        True
        >>> r = np.random.rand(10)
        >>> T2 = toeplitz(c, r)
        >>> T2
        Faust size 10x10, density 6.16, nnz_sum 616, 12 factor(s):
            - FACTOR 0 (complex) SPARSE, size 10x32, density 0.0625, nnz 20
            - FACTOR 1 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 2 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 3 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 4 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 5 (complex) SPARSE, size 32x32, density 0.03125, nnz 32
            - FACTOR 6 (complex) SPARSE, size 32x32, density 0.03125, nnz 32
            - FACTOR 7 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 8 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 9 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 10 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 11 (complex) SPARSE, size 32x10, density 0.0625, nnz 20

        >>> np.allclose(T2.toarray()[0, :], np.hstack((c[0:1], r[1:])))
        True
        >>> np.allclose(T2.toarray()[:, 0], c)
        True

    See also <a
    href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.toeplitz.htm">scipy.linalg.toeplitz</a>, pyfaust.circ, pyfaust.anticirc
    """
    if r is None:
        r = np.conjugate(c)
    m = len(c)
    n = len(r)
    N = int(2**np.ceil(np.log2(max(m, n))))
    c_ = np.hstack((c, np.zeros(N-m+1+N-n), r[:0:-1]))
    C = circ(c_)
    return C[:m, :n]

def eye(m, n=None, t='real', dev="cpu"):
    """
        Faust identity.

        Args:
          m: number of rows,
          n: number of columns, set to m by default.
          t: 'complex' to return a complex Faust, otherwise it's a real Faust.

        Examples:
            >>> from pyfaust import eye
            >>> eye(5)
            Faust size 5x5, density 0.2, nnz_sum 5, 1 factor(s):<br/>
            FACTOR 0 (real) SPARSE, size 5x5, density 0.2, nnz 5<br/>
            >>> eye(5, 4)
            Faust size 5x4, density 0.2, nnz_sum 4, 1 factor(s):<br/>
            FACTOR 0 (real) SPARSE, size 5x4, density 0.2, nnz 4<br/>
            >>> eye(5, t='complex')
            Faust size 5x4, density 0.2, nnz_sum 4, 1 factor(s):<br/>
            FACTOR 0 (complex) SPARSE, size 5x4, density 0.2, nnz 4<br/>
    """
    check_dev(dev)
    if(t not in ['complex', 'real', 'double', 'float']):
        raise ValueError("t must be 'real' (or 'double', 'float') or 'complex'")
    if(n == None): n = m
    if dev == "cpu":
        if t in ['real', 'double']:
            rF = Faust(core_obj=_FaustCorePy.FaustAlgoGenCPUDbl.eyeFaust(m, n))
        elif t == 'float':
            rF = Faust(core_obj=_FaustCorePy.FaustAlgoGenCPUFlt.eyeFaust(m, n))
        else:
            rF = Faust(core_obj=_FaustCorePy.FaustAlgoGenCPUCplxDbl.eyeFaust(m, n))
    elif dev.startswith("gpu"):
        if t in ['real', 'double']:
            rF = Faust(core_obj=_FaustCorePy.FaustAlgoGenGPUDbl.eyeFaust(m, n))
        elif t == 'float':
            raise TypeError("float is not yet supported by GPU wht")
        else:
            rF = Faust(core_obj=_FaustCorePy.FaustAlgoGenGPUCplxDbl.eyeFaust(m, n))
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

def rand_bsr(num_rows, num_cols, bnrows, bncols, num_factors=None, density=.1,
            dev='cpu', dtype='double'):
    """
    Generates a random Faust composed only of BSR matrices.

    Args:
        num_rows: the Faust number of rows.
        num_cols: the Faust number of columns.
        bnrows: the nonzero block number of rows (must divide num_rows).
        bncols: the nonzero block number of columns (must divide num_cols).
        num_factors: If it's an integer it will be the number of random factors to set in the Faust.
                    If num_factors is a tuple of 2 integers then the
                    number of factors will be set randomly between
                    num_factors[0] and num_factors[1] (inclusively).
                    If num_factors is None then 5 factors are generated.
        density: the Faust factor density (it determines the number of nonzero blocks). It must be between 0 and 1.
        dev: the device on which the Faust is created.
        dtype: the numpy dtype of the Faust.

    Example:
        >>> from pyfaust import rand_bsr
        >>> rand_bsr(100,100, 20, 10, num_factors=6) 
        Faust size 100x100, density 0.6, nnz_sum 6000, 6 factor(s): 
            - FACTOR 0 (double) BSR, size 100x100, density 0.1, nnz 1000
            - FACTOR 1 (double) BSR, size 100x100, density 0.1, nnz 1000
            - FACTOR 2 (double) BSR, size 100x100, density 0.1, nnz 1000
            - FACTOR 3 (double) BSR, size 100x100, density 0.1, nnz 1000
            - FACTOR 4 (double) BSR, size 100x100, density 0.1, nnz 1000
            - FACTOR 5 (double) BSR, size 100x100, density 0.1, nnz 1000

    <b>See also:</b> <a
    href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.bsr_matrix.html">bsr_matrix</a>, Faust.__init__, pyfaust.rand
    """
    if num_factors is None:
        min_num_factors = max_num_factors = 5
    elif isinstance(num_factors, int):
        min_num_factors = max_num_factors = num_factors
    elif isinstance(num_factors, tuple) and len(num_factors) == 2 and isinstance(num_factors[0], int) and isinstance(num_factors[1], int):
        min_num_factors = num_factors[0]
        max_num_factors = num_factors[1]
    else:
        raise ValueError('num_factors must be None, a int or a tuple of int')
    # sanity checks
    if num_rows != num_cols:
        raise ValueError('currently only random square BSR Fausts can be'
                         ' generated.')
    if num_rows % bnrows != 0 or num_cols % bncols != 0:
        raise ValueError('the size of blocks must evenly divide the size of Faust matrices')
    if dev.startswith('gpu') and bnrows != bncols:
        raise ValueError('currently only square blocks are supported on GPU.')
    if dev == "cpu":
        if dtype in ['float64', 'double']:
            rF = Faust(core_obj=_FaustCorePy.FaustAlgoGenCPUDbl.randBSRFaust(num_rows,
                                                                             num_cols,
                                                                             min_num_factors,
                                                                             max_num_factors,
                                                                             bnrows,
                                                                             bncols,
                                                                             density))
        elif dtype in ['float32', 'float']: # type == 'float'
            rF = Faust(core_obj=_FaustCorePy.FaustAlgoGenCPUFlt.randBSRFaust(num_rows,
                                                                             num_cols,
                                                                             min_num_factors,
                                                                             max_num_factors,
                                                                             bnrows,
                                                                             bncols,
                                                                             density))
        elif dtype in ['complex', 'complex128']:
            rF = Faust(core_obj=_FaustCorePy.FaustAlgoGenCPUCplxDbl.randBSRFaust(num_rows,
                                                                                 num_cols,
                                                                                 min_num_factors,
                                                                                 max_num_factors,
                                                                                 bnrows,
                                                                                 bncols,
                                                                                 density))
    elif dev.startswith("gpu"):
        if dtype in ['float64', 'double']:
            rF = Faust(core_obj=_FaustCorePy.FaustAlgoGenGPUDbl.randBSRFaust(num_rows,
                                                                             num_cols,
                                                                             min_num_factors,
                                                                             max_num_factors,
                                                                             bnrows,
                                                                             bncols,
                                                                             density))
        elif dtype in ['float32', 'float']: # type == 'float'
            rF = Faust(core_obj=_FaustCorePy.FaustAlgoGenGPUFlt.randBSRFaust(num_rows,
                                                                             num_cols,
                                                                             min_num_factors,
                                                                             max_num_factors,
                                                                             bnrows,
                                                                             bncols,
                                                                             density))
        elif dtype in ['complex', 'complex128']:
            rF = Faust(core_obj=_FaustCorePy.FaustAlgoGenGPUCplxDbl.randBSRFaust(num_rows,
                                                                                 num_cols,
                                                                                 min_num_factors,
                                                                                 max_num_factors,
                                                                                 bnrows,
                                                                                 bncols,
                                                                                 density))
    else:
        raise ValueError('Invalid device')
    return rF


def rand(num_rows, num_cols, num_factors=None, dim_sizes=None,
         density=None, fac_type='sparse',
         field='real', per_row=True, dev='cpu', dtype='double'):
    """
    Generates a random Faust.

        Args:
            num_rows: the Faust number of rows.
            num_cols: the Faust number of columns.
            num_factors: if it's an integer it is the number of random factors to set in the Faust.
                        If num_factors is a list or tuple of 2 integers then the
                        number of factors is set randomly between
                        num_factors[0] and num_factors[1] (inclusively).
                        Defaultly, num_factors is None, it means a 5 factors
                        long Faust is generated.
                        dim_sizes: if it's an integer all Faust factors
                        are square matrices (except maybe the first and last
                        ones, depending on num_rows and num_cols). The size of
                        the intermediary square factors is size_dims**2.
                        If it's a list or tuple of 2 integers then the
                        number of rows and columns are both
                        a random number between size_dims[0] and
                        size_dims[1] (inclusively).
                        Note that the first factor number of rows and the last
                        factor number of columns are always fixed (they are the
                        dimension sizes of the Faust: num_rows, num_cols arguments).
                        if dim_sizes is None then dim_sizes is defaultly [num_rows,
                        num_cols].
            density: the approximate density of factors. The
                    default value is such that each factor gets 5 non-zero
                    elements per row, if per_row is True, or per column otherwise.
                    It should be a floating point number greater than 0 and
                    lower or equal to 1.
                    A density of zero is equivalent to the default case.
            fac_type: the storage representation of factors. It must be
                    'sparse', 'dense' or 'mixed'. The latter designates a mix of dense and
                    sparse matrices in the generated Faust (the choice is made according
                    to a uniform distribution).
            field:  a str to set the Faust field: 'real' or 'complex'.
            per_row: if True the factor matrix is constructed per row
                    applying the density to each row. If False the construction is
                    made with the density applied on each column.
            dev: the device on which to create the Faust ('cpu' or 'gpu').
            dtype: the dtype of the Faust ('float32' and 'float128' == 'double'
                   are supported for fiel == 'real', only 'double' for field ==
                   'complex'.

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

    <b>See also</b> Faust.__init__, pyfaust.rand_bsr
    """
    check_dev(dev)
    field = field.lower()
    if field == 'real':
        is_real = True
        if dtype in [np.float64, 'float', 'double']:
            type = 'double'
        elif dtype in [np.float32, 'float32']:
            type = 'float'
    elif field == 'complex':
        if not dtype in [np.float64, 'float', 'double']:
            raise TypeError('Invalid dtype: only double are handled for for'
                            ' complex field')
        is_real = False
    else:
        raise ValueError('field argument must be either \'real\' or \'complex\'.')
    DENSE=0
    SPARSE=1
    MIXED=2
    REAL=3
    COMPLEX=4
    # set repr. type of factors
    if(not isinstance(fac_type, str) or fac_type not in ['sparse',
                                                         'dense',
                                                         'mixed']):
        raise ValueError('rand(): argument fac_type must be a'
                         ' str equal to \'sparse\','
                         ' \'dense\' or \'mixed\'.')

    fac_type_map = {"sparse" : SPARSE, "dense" : DENSE, "mixed" : MIXED}
    # set field of matrices/factors
    if(not isinstance(is_real, bool)):
        raise ValueError('rand(): argument is_real must be a'
                         'boolean.')
    if is_real:
        field = REAL
    else:
        field = COMPLEX
    if num_factors == None:
        num_factors = 5
    if((isinstance(num_factors,(list, tuple))) and
    len(num_factors) == 2):
        min_num_factors = num_factors[0]
        max_num_factors = num_factors[1]
    elif(isinstance(num_factors, int)):
        min_num_factors = max_num_factors = num_factors
    else:
        raise ValueError("rand(): num_factors argument must be an "
                         "integer or a list/tuple of 2 integers.")
    if dim_sizes == None:
        dim_sizes = [num_rows, num_cols]
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
    if not density:
        density = -1
    elif not isinstance(density, (float, int)):
        raise ValueError("rand(): density must be a float")
    density = float(density)
    if dev == "cpu":
        if field == REAL:
            if type == 'double':
                rF = Faust(core_obj=_FaustCorePy.FaustAlgoGenCPUDbl.randFaust(num_rows,
                                                                              num_cols,
                                                                              fac_type_map[fac_type], min_num_factors, max_num_factors,
                                                                              min_dim_size, max_dim_size, density, per_row))
            else: # type == 'float'
                rF = Faust(core_obj=_FaustCorePy.FaustAlgoGenCPUFlt.randFaust(num_rows,
                                                                              num_cols,
                                                                              fac_type_map[fac_type], min_num_factors, max_num_factors,
                                                                              min_dim_size, max_dim_size, density, per_row))
        elif field == COMPLEX:
            rF = Faust(core_obj=_FaustCorePy.FaustAlgoGenCPUCplxDbl.randFaust(num_rows,
                                                                           num_cols,
                                                                           fac_type_map[fac_type], min_num_factors, max_num_factors,
                                                                           min_dim_size, max_dim_size, density, per_row))
        # no else possible (see above)
    elif dev.startswith("gpu"):
        if field == REAL:
            if type == 'double':
                rF = Faust(core_obj=_FaustCorePy.FaustAlgoGenGPUDbl.randFaust(num_rows,
                                                                              num_cols,
                                                                              fac_type_map[fac_type], min_num_factors, max_num_factors,
                                                                              min_dim_size, max_dim_size, density, per_row))
            else: # type == float:
                rF = Faust(core_obj=_FaustCorePy.FaustAlgoGenGPUFlt.randFaust(num_rows,
                                                                           num_cols,
                                                                           fac_type_map[fac_type], min_num_factors, max_num_factors,
                                                                           min_dim_size, max_dim_size, density, per_row))
        elif field == COMPLEX:
            rF = Faust(core_obj=_FaustCorePy.FaustAlgoGenGPUCplxDbl.randFaust(num_rows,
                                                                           num_cols,
                                                                           fac_type_map[fac_type], min_num_factors, max_num_factors,
                                                                           min_dim_size, max_dim_size, density, per_row))
    return rF

def enable_gpu_mod(libpaths=None, backend='cuda', silent=False, fatal=False):
    """
    This function loads explicitly the gpu_mod library in memory.

    Normally it's not required to load the library manually, but it could be
    useful to set a non-default path or to diagnose a loading issue.

    Args:
        libpaths: the absolute or relative paths where to search the dynamic
                library (gm) to load. By default, it's none to auto-find the library
                (if possible).
        backend: the GPU backend to use, only 'cuda' is available for now.
        silent: if True nothing or almost will be displayed on loading,
                otherwise all messages are visible.
    """
    return _FaustCorePy.enable_gpu_mod(libpaths, backend, silent, fatal)

def is_gpu_mod_enabled():
    """
    Returns True if the gpu_mod plug-in has been loaded correctly, False otherwise.
    """
    return _FaustCorePy._is_gpu_mod_enabled()

def check_dev(dev):
    if dev.startswith('gpu'):
        if not is_gpu_mod_enabled():
            raise Exception('GPU device is not available on your environment.')
    elif dev != 'cpu':
        raise ValueError("dev must be 'cpu' or 'gpu[:id]'")

def _cplx2real_op(op):
    if pyfaust.isFaust(op):
        return Faust([_cplx2real_op(op.factors(i)) for i in range(op.numfactors())])
    else:
        rop = np.real(op)
        iop = np.imag(op)
        if isinstance(op, (csr_matrix, csc_matrix, coo_matrix, bsr_matrix)):
            vertcat = svstack
            horzcat = shstack
        elif isinstance(op, np.ndarray):
            vertcat = vstack
            horzcat = hstack
        else:
            raise TypeError('op must be a scipy sparse matrix or a np.ndarray')
        real_part = horzcat((rop, - iop))
        imag_part = horzcat((iop, rop))
        return vertcat((real_part, imag_part))


# experimental block start
# @PYTORCH_EXP_CODE@
# experimental block end

class FaustMulMode:
    """
    <b/> Enumeration class of all matrix chain multiplication methods available to multiply a Faust to a matrix or to compute Faust.toarray().

    These methods are used by Faust.optimize_time().

    NOTE: it's not advisable to use these methods directly. The user should use
    Faust.optimize_time() but the access is left open for experimental purpose.

    Examples:
        >>> from pyfaust import rand as frand
        >>> from numpy.random import rand
        >>> F = frand(100, 100, 5, [100, 1024])
        >>> F.m_faust.set_FM_mul_mode(FaustMulMode.DYNPROG) # method used to compute Faust-matrix product or Faust.toarray()
        >>> F*rand(F.shape[1], 512) # Faust-matrix mul. using method DYNPROG
        >>> F.toarray() # using the same method
    """
    ## \brief The default method, it computes the product from the right to the left.
    DEFAULT_L2R=0
    ## \brief This method implements the classic dynamic programming solution to the chain matrix problem.
    ##
    ## See https://en.wikipedia.org/wiki/Matrix_chain_multiplication#A_dynamic_programming_algorithm.
    ## Note that the standard method is extended in order to take into account the complexity of multiplications including a sparse matrix (because that's not the same cost than multiplying dense matrices).
    DYNPROG=5
    ## \brief This method computes the product of the matrix chain from the left to the right using the Torch C++ library (CPU backend).
    ##
    ## This method is only available for the specific packages pyfaust_torch.
    TORCH_CPU_L2R=8
    ## \brief This method is implemented using the Torch library and follows a greedy principle: it chooses to multiply the less costly product of two matrices at each step of the whole product computation.
    ##
    ## The computational cost depends on the matrix dimensions and the number of nonzeros (when a matrix is in sparse format).
	## This method is only available for the specific packages pyfaust_torch.
    TORCH_CPU_GREEDY=9
    ## \brief The same as TORCH_CPU_L2R except that torch::chain_matmul is used to
    ## compute in one call the intermediary product of dense contiguous
    ## factors, then the result is multiplied by sparse factors if any remains.
    ##
    ## torch::chain_matmul follows the dynamic programming principle as DYNPROG method does (but the former handles only dense matrices).
    ##
    ## References:
    ## https://pytorch.org/cppdocs/api/function_namespaceat_1aee491a9ff453b6033b4106516bc61a9d.html?highlight=chain_matmul
    ## https://pytorch.org/docs/stable/generated/torch.chain_matmul.html?highlight=chain_matmul#torch.chain_matmul
    ## This method is only available for the specific packages pyfaust_torch.
    TORCH_CPU_DENSE_DYNPROG_SPARSE_L2R=10
