## @package pyfaust.lazylinop @brief The pyfaust module for lazy linear operators.

import numpy as np
from scipy.sparse.linalg import LinearOperator

HANDLED_FUNCTIONS = {}

def isLazyLinearOp(obj):
    """
    Returns True if obj is a LazyLinearOp, False otherwise.
    """
    return LazyLinearOp.isLazyLinearOp(obj)

def aslazylinearoperator(obj):
    """
    Creates a LazyLinearOp based on the object obj which must be of a linear operator compatible type.

    NOTE: obj must support operations and attributes defined in the
    LazyLinearOp class.
    Any operation not supported would raise an exception at the evaluation
    time.

    Args:
        obj: the root object on which the LazyLinearOp is based (it could
        be a numpy array, a scipy matrix, a Faust object or almost any
        object that supports the same kind of functions).


    Returns:
        a LazyLinearOp instance based on obj.

    Example:
        >>> from pyfaust.lazylinop import asLazyLinearOp
        >>> import numpy as np
        >>> M = np.random.rand(10, 12)
        >>> lM = azlazylinearoperator(M)
        >>> twolM = lM + lM
        >>> twolM
        <pyfaust.lazylinop.LazyLinearOp at 0x7fcd7d7750f0>
        >>> import pyfaust as pf
        >>> F = pf.rand(10, 12)
        >>> lF = azlazylinearoperator(F)
        >>> twolF = lF + lF
        >>> twolF
        <pyfaust.lazylinop.LazyLinearOp at 0x7fcd7d774730>


    <b>See also:</b> pyfaust.rand.
    """
    return LazyLinearOp.create_from_op(obj)

def hstack(tup):
    """
    Concatenates tup objects horizontally.

    Args:
        tup: a tuple whose first argument is a LazyLinearOp and other must be
        compatible objects (numpy array, matrix, LazyLinearOp).

    Return:
        A LazyLinearOp resulting of the concatenation.
    """
    lop = tup[0]
    if isLazyLinearOp(lop):
        return lop.concatenate(*tup[1:], axis=1)
    else:
        raise TypeError('lop must be a LazyLinearOp')

def vstack(tup):
    """
    Concatenates tup objects vertically.

    Args:
        tup: a tuple whose first argument is a LazyLinearOp and other must be
        compatible objects (numpy array, matrix, LazyLinearOp).

    Return:
        A LazyLinearOp resulting of the concatenation.
    """
    lop = tup[0]
    if isLazyLinearOp(lop):
        return lop.concatenate(*tup[1:], axis=0)
    else:
        raise TypeError('lop must be a LazyLinearOp')

class LazyLinearOp(LinearOperator):
    """
    This class implements a lazy linear operator. A LazyLinearOp is a
    specialization of a <a
    href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.LinearOperator.html">scipy.linalg.LinearOperator</a>.

    The evaluation of any defined operation is delayed until a multiplication
    by a dense matrix/vector, a call of LazyLinearOp.toarray.

    To instantiate a LazyLinearOp look at
    pyfaust.lazylinop.aslazylinearoperator or
    pyfaust.lazylinop.LazyLinearOperator to instantiate from matmat/matvec
    functions.

    Warning: This code is in a beta status.
    """
    def __init__(self, lambdas, shape, root_obj=None, dtype=None):
        """
        Constructor. Not meant to be used directly.

        Args:
            lambdas: starting operations.
            shape: the initial shape of the operator.
            root_obj: the initial object the operator is based on.

        <b>See also:</b> pyfaust.lazylinop.aslazylinearoperator.
        """
        if root_obj is not None:
            if not hasattr(root_obj, 'ndim'):
                raise TypeError('The starting object to initialize a'
                                ' LazyLinearOperator must possess a ndim'
                                ' attribute.')
            if root_obj.ndim != 2:
                raise ValueError('The starting object to initialize a LazyLinearOp '
                                 'must have two dimensions, not: '+str(root_obj.ndim))
        self.lambdas = lambdas
        self._check_lambdas()
        self.shape = shape
        self._root_obj = root_obj
        self.dtype = dtype
        super(LazyLinearOp, self).__init__(self.dtype, self.shape)

    def _check_lambdas(self):
        if not isinstance(self.lambdas, dict):
            raise TypeError('lambdas must be a dict')
        keys = self.lambdas.keys()
        for k in ['@', 'H', 'T', 'slice']:
            if k not in keys:
                raise ValueError(k+' is a mandatory lambda, it must be set in'
                                 ' self.lambdas')

    @staticmethod
    def create_from_op(obj):
        """
        Alias of pyfaust.lazylinop.aslazylinearoperator.
        """
        lambdas = {'@': lambda op: obj @ op}
        lambdasT = {'@': lambda op: obj.T @ op}
        lambdasH = {'@': lambda op: obj.T.conj() @ op}
        lambdasC = {'@': lambda op: obj.conj() @ op}
        # set lambdas temporarily to None (to satisfy the ctor)
        # they'll be initialized later
        for l in [lambdas, lambdasT, lambdasH, lambdasC]:
            l['T'] = None
            l['H'] = None
            l['slice'] = None #TODO: rename slice to index
        lop = LazyLinearOp(lambdas, obj.shape, obj, dtype=obj.dtype)
        lopT = LazyLinearOp(lambdasT, (obj.shape[1], obj.shape[0]), obj, dtype=obj.dtype)
        lopH = LazyLinearOp(lambdasH, (obj.shape[1], obj.shape[0]), obj, dtype=obj.dtype)
        lopC = LazyLinearOp(lambdasC, obj.shape, obj, dtype=obj.dtype)

        # TODO: refactor with create_from_funcs (in ctor)
        lambdas['T'] = lambda: lopT
        lambdas['H'] = lambda: lopH
        lambdas['slice'] = lambda indices: LazyLinearOp._index_lambda(lop,
                                                                       indices)()
        lambdasT['T'] = lambda: lop
        lambdasT['H'] = lambda: lopC
        lambdasT['slice'] = lambda indices: LazyLinearOp._index_lambda(lopT,
                                                                        indices)()
        lambdasH['T'] = lambda: lopC
        lambdasH['H'] = lambda: lop
        lambdasH['slice'] = lambda indices: LazyLinearOp._index_lambda(lopH,
                                                                        indices)()
        lambdasC['T'] = lambda: lopH
        lambdasC['H'] = lambda: lopT
        lambdasC['slice'] = lambda indices: LazyLinearOp._index_lambda(lopC,
                                                                        indices)()

        return lop

    #TODO: function not used anywhere, delete it?
    @staticmethod
    def create_from_scalar(s, shape=None):
        """
        Returns a LazyLinearOp L created from the scalar s, such as each L[i, i] == s.
        """
        if shape is None:
            shape = (1, 1)
        a = np.ones(shape) * s
        return create_from_op(a)

    @staticmethod
    def create_from_funcs(matmat, rmatmat, shape, dtype=None):
        """
        Alias of ctor pyfaust.lazylinop.LazyLinearOperator.
        """
        from scipy.sparse import eye as seye

        #MX = lambda X: matmat(np.eye(shape[1])) @ X
        MX = lambda X: matmat(X)
        #MTX = lambda X: rmatmat(X.T).T
        MHX = lambda X: rmatmat(X)

        lambdas = {'@': MX}
        lambdasT = {'@': lambda op: rmatmat(op.real).conj() -
                    rmatmat(1j * op.imag).conj()}
        lambdasH = {'@': MHX}
        lambdasC = {'@': lambda op: matmat(op.real).conj() -
                    matmat(1j * op.imag)}
        # set lambdas temporarily to None (to satisfy the ctor)
        # they'll be initialized later
        for l in [lambdas, lambdasT, lambdasH, lambdasC]:
            l['T'] = None
            l['H'] = None
            l['slice'] = None

        lop = LazyLinearOp(lambdas, shape, dtype)
        lopT = LazyLinearOp(lambdasT, (shape[1], shape[0]), dtype)
        lopH = LazyLinearOp(lambdasH, (shape[1], shape[0]), dtype)
        lopC = LazyLinearOp(lambdasC, shape, dtype)

        lambdas['T'] = lambda: lopT
        lambdas['H'] = lambda: lopH
        lambdas['slice'] = lambda indices: LazyLinearOp._index_lambda(lop,
                                                                       indices)()
        lambdasT['T'] = lambda: lop
        lambdasT['H'] = lambda: lopC
        lambdasT['slice'] = lambda indices: LazyLinearOp._index_lambda(lopT,
                                                                        indices)()
        lambdasH['T'] = lambda: lopC
        lambdasH['H'] = lambda: lop
        lambdasH['slice'] = lambda indices: LazyLinearOp._index_lambda(lopH,
                                                                        indices)()
        lambdasC['T'] = lambda: lopH
        lambdasC['H'] = lambda: lopT
        lambdasC['slice'] = lambda indices: LazyLinearOp._index_lambda(lopC,
                                                                        indices)()

        return lop

    def _checkattr(self, attr):
        if self._root_obj is not None and not hasattr(self._root_obj, attr):
            raise TypeError(attr+' is not supported by the root object of this'
                            ' LazyLinearOp')

    def _index_lambda(lop, indices):
        from scipy.sparse import eye as seye
        s = lambda: \
                LazyLinearOp.create_from_op(seye(lop.shape[0],
                                                  format='csr')[indices[0]]) \
                @ lop @ LazyLinearOp.create_from_op(seye(lop.shape[1], format='csr')[:, indices[1]])
        return s

    @property
    def ndim(self):
        """
        Returns the number of dimensions of the LazyLinearOp (it is always 2).
        """
        return 2

    def transpose(self):
        """
        Returns the LazyLinearOp transpose.
        """
        self._checkattr('transpose')
        return self.lambdas['T']()

    @property
    def T(self):
        """
        Returns the LazyLinearOp transpose.
        """
        return self.transpose()

    def conj(self):
        """
        Returns the LazyLinearOp conjugate.
        """
        self._checkattr('conj')
        return self.H.T

    def conjugate(self):
        """
        Returns the LazyLinearOp conjugate.
        """
        return self.conj()

    def getH(self):
        """
        Returns the LazyLinearOp adjoint/transconjugate.
        """
        #self._checkattr('getH')
        return self.lambdas['H']()

    @property
    def H(self):
        """
        The LazyLinearOp adjoint/transconjugate.
        """
        return self.getH()

    def _adjoint(self):
        """
        Returns the LazyLinearOp adjoint/transconjugate.
        """
        return self.H

    def _slice(self, indices):
        return self.lambdas['slice'](indices)

    def __add__(self, op):
        """
        Returns the LazyLinearOp for self + op.

        Args:
            op: an object compatible with self for this binary operation.

        """
        self._checkattr('__add__')
        if not LazyLinearOp.isLazyLinearOp(op):
            op = LazyLinearOp.create_from_op(op)
        lambdas = {'@': lambda o: self @ o + op @ o,
                   'H': lambda: self.H + op.H,
                   'T': lambda: self.T + op.T,
                   'slice': lambda indices: self._slice(indices) + op._slice(indices)
                  }
        new_op = LazyLinearOp(lambdas=lambdas,
                              shape=tuple(self.shape),
                              root_obj=None)
        return new_op

    def __radd__(self, op):
        """
        Returns the LazyLinearOp for op + self.

        Args:
            op: an object compatible with self for this binary operation.

        """
        return self.__add__(op)

    def __iadd__(self, op):
        """
        Not Implemented self += op.
        """
        raise NotImplementedError(LazyLinearOp.__name__+".__iadd__")
# can't do as follows, it recurses indefinitely because of self.eval
#        self._checkattr('__iadd__')
#        self = LazyLinearOp(init_lambda=lambda:
#                              (self.eval()).\
#                              __iadd__(LazyLinearOp._eval_if_lazy(op)),
#                              shape=(tuple(self.shape)),
#                              root_obj=self._root_obj)
#        return self


    def __sub__(self, op):
        """
        Returns the LazyLinearOp for self - op.

        Args:
            op: an object compatible with self for this binary operation.

        """
        self._checkattr('__sub__')
        if not LazyLinearOp.isLazyLinearOp(op):
            op = LazyLinearOp.create_from_op(op)
        lambdas = {'@': lambda o: self @ o - op @ o,
                   'H': lambda: self.H - op.H,
                   'T': lambda: self.T - op.T,
                   'slice': lambda indices: self._slice(indices) - op._slice(indices)
                  }
        new_op = LazyLinearOp(lambdas=lambdas,
                              shape=tuple(self.shape),
                              root_obj=None)
        return new_op


    def __rsub__(self, op):
        """
        Returns the LazyLinearOp for op - self.

        Args:
            op: an object compatible with self for this binary operation.

        """
        self._checkattr('__rsub__')
        if not LazyLinearOp.isLazyLinearOp(op):
            op = LazyLinearOp.create_from_op(op)
        lambdas = {'@': lambda o: op @ o - self @ o,
                   'H': lambda: op.H - self.H,
                   'T': lambda: op.T - self.T,
                   'slice': lambda indices: op._slice(indices) - self._slice(indices)
                  }
        new_op = LazyLinearOp(lambdas=lambdas,
                              shape=self.shape,
                              root_obj=None)
        return new_op

    def __isub__(self, op):
        """
        Not implemented self -= op.
        """
        raise NotImplementedError(LazyLinearOp.__name__+".__isub__")
# can't do as follows, it recurses indefinitely because of self.eval
#        self._checkattr('__isub__')
#        self = LazyLinearOp(init_lambda=lambda:
#                              (self.eval()).\
#                              __isub__(LazyLinearOp._eval_if_lazy(op)),
#                              shape=(tuple(self.shape)),
#                              root_obj=self._root_obj)
#        return self


    def __truediv__(self, s):
        """
        Returns the LazyLinearOp for self / s.

        Args:
            s: a scalar.

        """
        new_op = self * 1/s
        return new_op

    def __itruediv__(self, op):
        """
        Not implemented self /= op.
        """
        raise NotImplementedError(LazyLinearOp.__name__+".__itruediv__")
# can't do as follows, it recurses indefinitely because of self.eval
#
#        self._checkattr('__itruediv__')
#        self = LazyLinearOp(init_lambda=lambda:
#                              (self.eval()).\
#                              __itruediv__(LazyLinearOp._eval_if_lazy(op)),
#                              shape=(tuple(self.shape)),
#                              root_obj=self._root_obj)
#        return self

    def _sanitize_matmul(self, op, swap=False):
        self._checkattr('__matmul__')
        if not hasattr(op, 'shape'):
            raise TypeError('op must have a shape attribute')
        if not hasattr(op, 'ndim'):
            raise TypeError('op must have a ndim attribute')
        if op.ndim == 1 and (self.shape[0] if swap else self.shape[1]) != op.size or op.ndim == 2 and (swap and
                                                                       op.shape[1]
                                                                       !=
                                                                       self.shape[0]
                                                                       or not
                                                                       swap and
                                                                       self.shape[1]
                                                                       !=
                                                                       op.shape[0]):
            raise ValueError('dimensions must agree')

    def __matmul__(self, op):
        """
        Computes self @ op as a sparse matrix / dense array or as a LazyLinearOp depending on the type of op.

        Args:
            op: an object compatible with self for this binary operation.

        Returns:
            If op is an numpy array or a scipy matrix the function returns (self @
            op) as a numpy array or a scipy matrix. Otherwise it returns the
            LazyLinearOp for the multiplication self @ op.

        """
        from scipy.sparse import issparse
        self._sanitize_matmul(op)
        if isinstance(op, np.ndarray) or issparse(op):
            if op.ndim == 1 and self._root_obj is not None:
                res = self.lambdas['@'](op.reshape(op.size, 1)).ravel()
            else:
                res = self.lambdas['@'](op)
        else:
            if not LazyLinearOp.isLazyLinearOp(op):
                op = LazyLinearOp.create_from_op(op)
            lambdas = {'@': lambda o: self @ (op @ o),
                       'H': lambda: op.H @ self.H,
                       'T': lambda: op.T @ self.T,
                       'slice': lambda indices: self._slice((indices[0],
                                                            slice(0,
                                                                  self.shape[1])))\
                       @ op._slice((slice(0, op.shape[0]), indices[1]))
                      }
            res = LazyLinearOp(lambdas=lambdas,
                               shape=(self.shape[0], op.shape[1]),
                               root_obj=None)
#            res = LazyLinearOp.create_from_op(super(LazyLinearOp,
#                                                     self).__matmul__(op))
        return res

    def dot(self, op):
        """
        Alias of LazyLinearOp.__matmul__.
        """
        return self.__matmul__(op)

    def matvec(self, op):
        """
        This function is an alias of self @ op, where the multiplication might
        be specialized for op a vector (depending on how self has been defined
        ; upon on a operator object or through a matvec/matmat function).


        <b>See also:</b> pyfaust.lazylinop.LazyLinearOperator.
        """
        if not hasattr(op, 'shape') or not hasattr(op, 'ndim'):
            raise TypeError('op must have shape and ndim attributes')
        if op.ndim > 2 or op.ndim == 0:
            raise ValueError('op.ndim must be 1 or 2')
        if op.ndim != 1 and op.shape[0] != 1 and op.shape[1] != 1:
            raise ValueError('op must be a vector -- attribute ndim to 1 or'
                             ' shape[0] or shape[1] to 1')
        return self.__matmul__(op)

    def _rmatvec(self, op):
        """
        Returns self^H @ op, where self^H is the conjugate transpose of A.

        Returns:
            It might be a LazyLinearOp or an array depending on the op type
            (cf. pyfaust.lazylinop.LazyLinearOp.__matmul__).
        """
        # LinearOperator need.
        return self.T.conj() @ op

    def _matmat(self, op):
        """
        Alias of LazyLinearOp.__matmul__.
        """
        # LinearOperator need.
        if not hasattr(op, 'shape') or not hasattr(op, 'ndim'):
            raise TypeError('op must have shape and ndim attributes')
        if op.ndim > 2 or op.ndim == 0:
            raise ValueError('op.ndim must be 1 or 2')
        return self.__matmul__(op)

    def _rmatmat(self, op):
        """
        Returns self^H @ op, where self^H is the conjugate transpose of A.

        Returns:
            It might be a LazyLinearOp or an array depending on the op type
            (cf. pyfaust.lazylinop.LazyLinearOp.__matmul__).
        """
        # LinearOperator need.
        return self.T.conj() @ op

    def __imatmul__(self, op):
        """
        Not implemented self @= op.
        """
        raise NotImplementedError(LazyLinearOp.__name__+".__imatmul__")

    def __rmatmul__(self, op):
        """
        Returns op @ self which can be a LazyLinearOp or an array depending on op type.

        Args:
            op: an object compatible with self for this binary operation.

        <b>See also:</b> pyfaust.lazylinop.LazyLinearOp.__matmul__)
        """
        self._checkattr('__rmatmul__')
        from scipy.sparse import issparse
        self._sanitize_matmul(op, swap=True)
        if isinstance(op, np.ndarray) or issparse(op):
            res = op @ self.toarray()
        else:
            if not LazyLinearOp.isLazyLinearOp(op):
                op = LazyLinearOp.create_from_op(op)
            lambdas = {'@': lambda o: (op @ self) @ o,
                       'H': lambda: self.H @ op.H,
                       'T': lambda: self.T @ op.T,
                       'slice': lambda indices: (op @ self)._slice(indices)
                      }
            res = LazyLinearOp(lambdas=lambdas,
                               shape=(op.shape[0], self.shape[1]),
                               root_obj=None)
        return res

    def __mul__(self, other):
        """
        Returns the LazyLinearOp for self * other if other is a scalar
        otherwise returns self @ other.

        Args:
            other: a scalar or a vector/array.

        <b>See also:</b> pyfaust.lazylinop.LazyLinearOp.__matmul__)
        """
        self._checkattr('__mul__')
        from scipy.sparse import eye, diags
        if np.isscalar(other):
            S = eye(self.shape[1], format='csr') * other
            lop = LazyLinearOp.create_from_op(S)
            new_op = self @ lop
        else:
            new_op = self @ other
        return new_op

    def __rmul__(self, other):
        """
        Returns other * self.

        Args:
            other: a scalar or a vector/array.

        """
        if np.isscalar(other):
            return self * other
        else:
            return other @ self


    def __imul__(self, op):
        """
        Not implemented self *= op.
        """
        raise NotImplementedError(LazyLinearOp.__name__+".__imul__")

    def toarray(self):
        """
        Returns self as a numpy array.
        """
        #from scipy.sparse import eye
        #return self @ eye(self.shape[1], self.shape[1], format='csr')
        # don't use csr because of function based LazyLinearOp
        # (e.g. scipy fft receives only numpy array)
        return self @ np.eye(self.shape[1])

    def __getitem__(self, indices):
        """
        Returns the LazyLinearOp for slicing/indexing.

        Args:
            indices: array of length 1 or 2 which elements must be slice, integer or
            Ellipsis (...). Note that using Ellipsis for more than two indices
            is normally forbidden.

        """
        self._checkattr('__getitem__')
        if isinstance(indices, tuple) and len(indices) == 2 and isinstance(indices[0], int) and isinstance(indices[1], int):
            return self.toarray().__getitem__(indices)
        elif isinstance(indices, slice) or isinstance(indices[0], slice) and \
                isinstance(indices[0], slice):
            return self._slice(indices)
        else:
            return self._slice(indices)

    def _newshape_getitem(self, indices):
        empty_lop_except = Exception("Cannot create an empty LazyLinearOp.")
        if isinstance(indices, (np.ndarray, list)):
            return (len(indices), self.shape[1])
        elif indices == Ellipsis:
            return self.shape
        elif isinstance(indices,int):
            # self[i] is a row
            return (1, self.shape[1])
        elif isinstance(indices, slice):
            #self[i:j] a group of contiguous rows
            start, stop, step = indices.start, indices.stop, indices.step
            if stop is None:
                stop = self.shape[0]
            if start is None:
                start = 0
            if step is None:
                step = 1
            return ((stop - start) // step, self.shape[1])
        elif isinstance(indices, tuple):
            if len(indices) == 1:
                return self._newshape_getitem(indices[0])
            elif len(indices) == 2:
                if(isinstance(indices[0], int) and isinstance(indices[1],int)):
                    # item
                    return (1, 1)
            else:
                raise IndexError('Too many indices.')

            if indices[0] == Ellipsis:
                if indices[1] == Ellipsis:
                    raise IndexError('an index can only have a single ellipsis '
                                     '(\'...\')')
                else:
                    # all rows
                    new_shape = self.shape
            elif isinstance(indices[0], int):
                # line F[i]
                new_shape = (1, self.shape[1])
            elif isinstance(indices[0], slice):
                start, stop, step = indices[0].start, indices[0].stop, indices[0].step
                if stop is None:
                    stop = self.shape[0]
                if start is None:
                    start = 0
                if step is None:
                    step = 1
                new_shape = ((stop - start) // step, self.shape[1])
            elif isinstance(indices[0], (list, np.ndarray)):
                if len(indices[0]) == 0: raise empty_lop_except
                new_shape = (len(indices[0]), self.shape[1])
            else:
                 raise idx_error_exception

            if indices[1] == Ellipsis:
                # all columns
                new_shape = self.shape
            elif isinstance(indices[1], int):
                # col F[:, i]
                new_shape = (new_shape[0], 1)
            elif isinstance(indices[1], slice):
                start, stop, step = indices[1].start, indices[1].stop, indices[1].step
                if stop is None:
                    stop = self.shape[1]
                if start is None:
                    start = 0
                if step is None:
                    step = 1
                new_shape = (new_shape[0], (stop - start) // step)
            elif isinstance(indices[1], (list, np.ndarray)):
                if len(indices[1]) == 0: raise empty_lop_except
                new_shape = (new_shape[0], len(indices[1]))
            else:
                 raise idx_error_exception
            return new_shape

    def concatenate(self, *ops, axis=0):
        """
        Returns the LazyLinearOp for the concatenation of self and op.

        Args:
            axis: axis of concatenation (0 for rows, 1 for columns).
        """
        from pyfaust import concatenate as cat
        nrows = self.shape[0]
        ncols = self.shape[1]
        out = self
        for op in ops:
            if axis == 0:
                out = out.vstack(op)
            elif axis == 1:
                out = out.hstack(op)
            else:
                raise ValueError('axis must be 0 or 1')
        return out

    def _vstack_slice(self, op, indices):
        rslice = indices[0]
        if rslice.step is not None:
            raise ValueError('Can\'t handle non-contiguous slice -- step > 1')
        if rslice.start == None:
            rslice = slice(0, rslice.stop, rslice.step)
        if rslice.stop == None:
            rslice = slice(rslice.start, self.shape[0] + op.shape[0], rslice.step)
        if rslice.stop > self.shape[0] + op.shape[0]:
            raise ValueError('Slice overflows the row dimension')
        if rslice.start >= 0 and rslice.stop <= self.shape[0]:
            # the slice is completly in self
            return lambda: self._slice(indices)
        elif rslice.start >= self.shape[0]:
            # the slice is completly in op
            return lambda: op._slice((slice(rslice.start - self.shape[0],
                                            rslice.stop - self.shape[0]) ,indices[1]))
        else:
            # the slice is overlapping self and op
            self_slice = self._slice((slice(rslice.start, self.shape[0]), indices[1]))
            op_slice = self._slice((slice(0, rslice.end - self.shape[0]), indices[1]))
            return lambda: self_slice.vstack(op_slice)

    def _vstack_mul_lambda(self, op, o):
        from scipy.sparse import issparse
        mul_mat = lambda : np.vstack((self @ o, op @ o))
        mul_vec = lambda : mul_mat().ravel()
        mul_mat_vec = lambda : mul_vec() if o.ndim == 1 else mul_mat()
        mul = lambda: mul_mat_vec() if isinstance(o, np.ndarray) \
                or issparse(o) else self.vstack(op) @ o
        return mul


    def vstack(self, op):
        """
        See pyfaust.lazylinop.vstack.
        """
        if self.shape[1] != op.shape[1]:
            raise ValueError('self and op numbers of columns must be the'
                             ' same')
        if not LazyLinearOp.isLazyLinearOp(op):
            op = LazyLinearOp.create_from_op(op)
        lambdas = {'@': lambda o: self._vstack_mul_lambda(op, o)(),
                   'H': lambda: self.H.hstack(op.H),
                   'T': lambda: self.T.hstack(op.T),
                   'slice': lambda indices: self._vstack_slice(op, indices)()
                  }
        new_op = LazyLinearOp(lambdas=lambdas,
                              shape=(self.shape[0] + op.shape[0], self.shape[1]),
                              root_obj=None)
        return new_op

    def _hstack_slice(self, op, indices):
        cslice = indices[1]
        if cslice.step is not None:
            raise ValueError('Can\'t handle non-contiguous slice -- step > 1')
        if cslice.stop > self.shape[1] + op.shape[1]:
            raise ValueError('Slice overflows the row dimension')
        if cslice.start >= 0 and cslice.stop <= self.shape[1]:
            # the slice is completly in self
            return lambda: self._slice(indices)
        elif cslice.start >= self.shape[1]:
            # the slice is completly in op
            return lambda: op._slice((indices[0], slice(cslice.start - self.shape[1],
                                            cslice.stop - self.shape[1])))
        else:
            # the slice is overlapping self and op
            self_slice = self._slice((indices[0], slice(cslice.start, self.shape[1])))
            op_slice = self._slice((indices[0], slice(0, cslice.end - self.shape[1])))
            return lambda: self_slice.vstack(op_slice)

    def _hstack_mul_lambda(self, op, o):
        from scipy.sparse import issparse
        if isinstance(o, np.ndarray) or issparse(o):
            if o.ndim == 1:
                return lambda: self @ o[:self.shape[1]] + op @ o[self.shape[1]:]
            else:
                return lambda: self @ o[:self.shape[1],:] + op @ o[self.shape[1]:, :]
        else:
            return lambda: \
                self @ o._slice((slice(0, self.shape[1]), slice(0,
                                                                o.shape[1]))) \
                    + op @ o._slice((slice(self.shape[1], o.shape[0]), slice(0, o.shape[1])))

    def hstack(self, op):
        """
        See pyfaust.lazylinop.hstack.
        """
        if self.shape[0] != op.shape[0]:
            raise ValueError('self and op numbers of rows must be the'
                             ' same')
        if not LazyLinearOp.isLazyLinearOp(op):
            op = LazyLinearOp.create_from_op(op)
        lambdas = {'@': lambda o: self._hstack_mul_lambda(op, o)(),
               'H': lambda: self.H.vstack(op.H),
               'T': lambda: self.T.vstack(op.T),
               'slice': lambda indices: self._hstack_slice(op, indices)()
              }
        new_op = LazyLinearOp(lambdas=lambdas,
                              shape=(self.shape[0], self.shape[1]
                                    + op.shape[1]),
                              root_obj=None)
        return new_op

    @property
    def real(self):
        """
        Returns the LazyLinearOp for real.
        """
        from scipy.sparse import issparse
        lambdas = {'@': lambda o: (self @ o.real).real + \
                   (self @ o.imag * 1j).real if isinstance(o, np.ndarray) \
                   or issparse(o) else real(self @ o),
                   'H': lambda: self.T.real,
                   'T': lambda: self.T.real,
                   'slice': lambda indices: self._slice(indices).real
                  }
        new_op = LazyLinearOp(lambdas=lambdas,
                           shape=tuple(self.shape),
                           root_obj=None)
        return new_op

    @property
    def imag(self):
        """
        Returns the LazyLinearOp for imag.
        """
        from scipy.sparse import issparse
        lambdas = {'@': lambda o: (self @ o.real).imag + \
                   1j * (self @ o.imag).imag if isinstance(o, np.ndarray) \
                   or issparse(o) else imag(self @ o),
                   'H': lambda: self.T.imag,
                   'T': lambda: self.T.imag,
                   'slice': lambda indices: self._slice(indices).imag
                  }
        new_op = LazyLinearOp(lambdas=lambdas,
                           shape=tuple(self.shape),
                           root_obj=None)
        return new_op


    @staticmethod
    def isLazyLinearOp(obj):
        """
        Returns True if obj is a LazyLinearOp, False otherwise.
        """
        return isinstance(obj, LazyLinearOp)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method == '__call__':
            if str(ufunc) == "<ufunc 'matmul'>" and len(inputs) >= 2 and \
               LazyLinearOp.isLazyLinearOp(inputs[1]):
                return inputs[1].__rmatmul__(inputs[0])
            elif str(ufunc) == "<ufunc 'multiply'>" and len(inputs) >= 2 and \
               LazyLinearOp.isLazyLinearOp(inputs[1]):
                return inputs[1].__rmul__(inputs[0])
            elif str(ufunc) == "<ufunc 'add'>" and len(inputs) >= 2 and \
                    LazyLinearOp.isLazyLinearOp(inputs[1]):
                return inputs[1].__radd__(inputs[0])
            elif str(ufunc) == "<ufunc 'subtract'>" and len(inputs) >= 2 and \
                    LazyLinearOp.isLazyLinearOp(inputs[1]):
                return inputs[1].__rsub__(inputs[0])
        elif method == 'reduce':
#            # not necessary numpy calls Faust.sum
#            if ufunc == "<ufunc 'add'>":
#                if len(inputs) == 1 and pyfaust.isLazyLinearOp(inputs[0]):
#                    #return inputs[0].sum(*inputs[1:], **kwargs)
#                else:
            return NotImplemented

    def __array__(self, *args, **kwargs):
        return self

    def __array_function__(self, func, types, args, kwargs):
        if func not in HANDLED_FUNCTIONS:
            return NotImplemented
        # Note: this allows subclasses that don't override
        # __array_function__ to handle Faust objects
        if not all(issubclass(t, LazyLinearOp) for t in types):
            return NotImplemented
        return HANDLED_FUNCTIONS[func](*args, **kwargs)

def LazyLinearOperator(shape, **kwargs):
    """
    Returns a LazyLinearOp defined by shape and at least a matvec or a matmat function.

    Args:
        shape: (tuple) dimensions (M, N).
        matvec: (callable) returns A * v (v a vector).
        rmatvec: (callable) returns A^H * v (v a vector of size N).
        matmat: (callable) returns A * V (V a dense matrix of dimensions (N, K)).
        rmatmat: (callable) returns A^H * V (V a dense matrix of dimensions (M, K)).
        dtype: data type of the matrix (can be None).

    Example:
        >>> # In this example we create a LazyLinearOp for the DFT using the fft from scipy
        >>> import numpy as np
        >>> from scipy.fft import fft, ifft
        >>> from pyfaust.lazylinop import LazyLinearOperator
        >>> lfft = LazyLinearOperator((8, 8), matmat=lambda x: fft(x, axis=0), rmatmat=lambda x: 8 * ifft(x, axis=0))
        >>> x = np.random.rand(8)
        >>> np.allclose(lfft * x, fft(x))
        True

    """
    matvec, rmatvec, matmat, rmatmat = [None for i in range(4)]
    def callable_err(k):
        return TypeError(k+' in kwargs must be a callable/function')
    for k in kwargs.keys():
        if k != 'dtype' and not callable(kwargs[k]):
            raise callable_err(k)
    if 'matvec' in kwargs.keys():
        matvec = kwargs['matvec']
    if 'rmatvec' in kwargs.keys():
        rmatvec = kwargs['rmatvec']
    if 'matmat' in kwargs.keys():
        matmat = kwargs['matmat']
    if 'rmatmat' in kwargs.keys():
        rmatmat = kwargs['rmatmat']
    if 'dtype' in kwargs.keys():
        dtype = kwargs['dtype']
    else:
        dtype = None

    if matvec is None and matmat is None:
        raise ValueError('At least a matvec or a matmat function must be'
                         ' passed in kwargs.')

    def _matmat(M, _matvec):
        if M.ndim == 1:
            return _matvec(M)

        out = np.empty((shape[0], M.shape[1]), dtype=dtype if dtype is not None
                      else M.dtype)
        for i in range(M.shape[1]):
            out[:, i] = _matvec(M[:,i])
        return out

    if matmat is None:
        matmat = lambda M: _matmat(M, matvec)

    if rmatmat is None and rmatvec is not None:
        rmatmat = lambda M: _matmat(M, rmatvec)

    return LazyLinearOp.create_from_funcs(matmat, rmatmat, shape, dtype=dtype)

def kron(A, B):
    """
    Returns the LazyLinearOp(Kron) for the Kronecker product A x B.

    Note: this specialization is particularly optimized for multiplying the
    operator by a vector.

    Args:
        A: LinearOperator (scaling factor),
        B: LinearOperator (block factor).

    Example:
        >>> from pyfaust.lazylinop import kron as lkron
        >>> import numpy as np
        >>> from pyfaust import rand
        >>> A = np.random.rand(100, 100)
        >>> B = np.random.rand(100, 100)
        >>> AxB = np.kron(A,B)
        >>> lAxB = lkron(A, B)
        >>> x = np.random.rand(AxB.shape[1], 1)
        >>> print(np.allclose(AxB@x, lAxB@x))
        True
        >>> from timeit import timeit
        >>> timeit(lambda: AxB @ x, number=10)
        0.4692082800902426
        >>> timeit(lambda: lAxB @ x, number=10)
        0.03464869409799576

    <b>See also:</b> numpy.kron.
    """
    def _kron(A, B, shape, op):
        from threading import Thread
        from multiprocessing import cpu_count
        from os import environ
        from pyfaust import isFaust
        #LazyLinearOp._sanitize_matmul(op) # TODO
        if hasattr(op, 'reshape') and hasattr(op, '__matmul__') and hasattr(op,
                                                                            '__getitem__'):
            if isFaust(A) or isFaust(B):
                parallel = False # e.g. for A, B Fausts in R^100x100 and op 128 columns
                # it was found that the sequential computation is faster
            else:
                parallel = True
            if 'KRON_PARALLEL' in environ:
                parallel = environ['KRON_PARALLEL'] == '1'
            nthreads = cpu_count() // 2
            if op.ndim == 1:
                op = op.reshape((op.size, 1))
                one_dim = True
            else:
                one_dim = False
            res = np.empty((shape[0], op.shape[1]))
            def out_col(j, ncols):
                for j in range(j, min(j + ncols, op.shape[1])):
                    op_mat = op[:, j].reshape((A.shape[1], B.shape[1]))
                    res[:, j] = (A @ op_mat @ B.T).reshape(shape[0])
            ncols = op.shape[1]
            if parallel:
                t = []
                cols_per_thread = ncols // nthreads
                if cols_per_thread * nthreads < ncols:
                    cols_per_thread += 1
                while len(t) < nthreads:
                    t.append(Thread(target=out_col, args=(cols_per_thread *
                                                          len(t),
                                                          cols_per_thread)))
                    t[-1].start()

                for j in range(nthreads):
                    t[j].join()
            else:
                out_col(0, ncols)
            if one_dim:
                res = res.ravel()
        else:
            raise TypeError('op must possess reshape, __matmul__ and'
                            ' __getitem__ attributes to be multiplied by a'
                            ' Kronecker LazyLinearOp (use toarray on the'
                            ' latter to multiply by the former)')
        return res
    shape = (A.shape[0] * B.shape[0], A.shape[1] * B.shape[1])
    return LazyLinearOperator(shape, matmat=lambda x: _kron(A, B, shape, x),
                              rmatmat=lambda x : _kron(A.T.conj(), B.T.conj(),
                                                       (shape[1], shape[0]), x))
