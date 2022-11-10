## @package pyfaust.lazylinop @brief The pyfaust module for lazy linear operators.

import numpy as np
from scipy.sparse.linalg import LinearOperator

HANDLED_FUNCTIONS = {}

class LazyLinearOp(LinearOperator):
    """
    This class implements a lazy linear operator. A LazyLinearOp is a
    specialization of a <a
    href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.LinearOperator.html">scipy.linalg.LinearOperator</a>.

    The evaluation of any defined operation is delayed until a multiplication
    by a dense matrix/vector, a call of LazyLinearOp.toarray or an explicit
    evaluation through LazyLinearOp.eval.

    To instantiate a LazyLinearOp look at pyfaust.lazylinop.asLazyLinearOp.

    Warning: This code is in a beta status.
    """
    def __init__(self, init_lambda, shape, root_obj):
        """
        Constructor. Not meant to be used directly.

        Args:
            init_lambda: starting operation.
            shape: the initial shape of the operator.
            root_obj: the initial object the operator is based on.

        <b>See also:</b> LazyLinearOp.create, pyfaust.lazylinop.asLazyLinearOp.
        """
        if not hasattr(root_obj, 'ndim'):
            raise TypeError('The starting object to initialize a'
                            ' LazyLinearOperator must possess a ndim'
                            ' attribute.')
        if root_obj.ndim != 2:
            raise ValueError('The starting object to initialize a LazyLinearOp '
                             'must have two dimensions, not: '+str(root_obj.ndim))

        self._lambda_stack = init_lambda
        self.shape = shape
        self._root_obj = root_obj
        self.dtype = None
        super(LazyLinearOp, self).__init__(self.dtype, self.shape)

    @staticmethod
    def create(obj):
        """
        Alias of pyfaust.lazylinop.asLazyLinearOp.

        Args:
            obj: cf. pyfaust.lazylinop.asLazyLinearOp.

        Returns:
            cf. pyfaust.lazylinop.asLazyLinearOp.

        Example:
            >>> from pyfaust.lazylinop import LazyLinearOp
            >>> import numpy as np
            >>> M = np.random.rand(10, 12)
            >>> lM = LazyLinearOp.create(M)
            >>> twolM = lM + lM
            >>> twolM
            <pyfaust.lazylinop.LazyLinearOp at 0x7fcd7d7750f0>
            >>> import pyfaust as pf
            >>> F = pf.rand(10, 12)
            >>> lF = LazyLinearOp.create(F)
            >>> twolF = lF + lF
            >>> twolF
            <pyfaust.lazylinop.LazyLinearOp at 0x7fcd7d774730>


        <b>See also:</b> pyfaust.rand.
        """
        return LazyLinearOp(lambda:obj, obj.shape, obj)

    def eval(self):
        """
        Evaluate the LazyLinearOp. All stacked virtual operations are evaluated.

        Returns:
            a linear operator object (whose type depends of the LazyLinearOp
            initialization through pyfaust.lazylinop.asLazyLinearOp and the operations made upon this object).

        Example:
            >>> import pyfaust as pf
            >>> from pyfaust.lazylinop import LazyLinearOp
            >>> import numpy as np
            >>> F = pf.rand(10, 12)
            >>> lF = LazyLinearOp.create(F)
            >>> twolF = 2 * lF
            >>> # twolF is a LazyLinearOp
            >>> # it is based on a Faust object
            >>> # so the evalution returns a Faust
            >>> twolF.eval()
            Faust size 10x12, density 2.03333, nnz_sum 244, 5 factor(s):
                - FACTOR 0 (double) SPARSE, size 10x10, density 0.5, nnz 50, addr: 0x562ec83e0a20
                - FACTOR 1 (double) SPARSE, size 10x10, density 0.5, nnz 50, addr: 0x562ec83be940
                - FACTOR 2 (double) SPARSE, size 10x12, density 0.333333, nnz 40, addr: 0x562ec83e2fa0
                - FACTOR 3 (double) SPARSE, size 12x11, density 0.454545, nnz 60, addr: 0x562ec82c66e0
                - FACTOR 4 (double) SPARSE, size 11x12, density 0.333333, nnz 44, addr: 0x562ec8330850

            >>> np.allclose(twolF.eval().toarray(), (2*F).toarray())
            True

        """
        return self._lambda_stack()

    def _checkattr(self, attr):
        if not hasattr(self._root_obj, attr):
            raise TypeError(attr+' is not supported by the root object of this'
                            ' LazyLinearOp')

    @staticmethod
    def _eval_if_lazy(o):
        return o.eval() if isLazyLinearOp(o) else o

    @property
    def ndim(self):
        """
        Returns the number of dimensions of the LazyLinearOp (it always 2).
        """
        return 2

    def transpose(self):
        """
        Returns the LazyLinearOp transpose.
        """
        self._checkattr('transpose')
        new_op = LazyLinearOp(init_lambda=lambda:
                                (self.eval()).transpose(),
                                shape=(self.shape[1], self.shape[0]),
                                root_obj=self._root_obj)
        return new_op

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
        new_op = LazyLinearOp(init_lambda=lambda:
                                (self.eval()).conj(),
                                shape=self.shape,
                                root_obj=self._root_obj)
        return new_op

    def conjugate(self):
        """
        Returns the LazyLinearOp conjugate.
        """
        return self.conj()

    def getH(self):
        """
        Returns the LazyLinearOp adjoint/transconjugate.
        """
        self._checkattr('getH')
        new_op = LazyLinearOp(init_lambda=lambda:
                                (self.eval()).getH(),
                                shape=(self.shape[1], self.shape[0]),
                                root_obj=self._root_obj)
        return new_op

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

    def __add__(self, op):
        """
        Returns the LazyLinearOp for self + op.

        Args:
            op: an object compatible with self for this binary operation.

        """
        self._checkattr('__add__')
        new_op = LazyLinearOp(init_lambda=lambda:
                                self.eval() + LazyLinearOp._eval_if_lazy(op),
                                shape=(tuple(self.shape)),
                                root_obj=self._root_obj)
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
        new_op = LazyLinearOp(init_lambda=lambda: self.eval() - LazyLinearOp._eval_if_lazy(op),
                                shape=(tuple(self.shape)),
                                root_obj=self._root_obj)
        return new_op

    def __rsub__(self, op):
        """
        Returns the LazyLinearOp for op - self.

        Args:
            op: an object compatible with self for this binary operation.

        """
        self._checkattr('__rsub__')
        new_op = LazyLinearOp(init_lambda=lambda:
                                LazyLinearOp._eval_if_lazy(op) -
                                self.eval(),
                                shape=(tuple(self.shape)),
                                root_obj=self._root_obj)
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


    def __truediv__(self, op):
        """
        Returns the LazyLinearOp for self / op.

        Args:
            op: an object compatible with self for this binary operation.

        """
        self._checkattr('__truediv__')
        new_op = LazyLinearOp(init_lambda=lambda:
                                self.eval() / LazyLinearOp._eval_if_lazy(op),
                                shape=(tuple(self.shape)),
                                root_obj=self._root_obj)
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

    def _sanitize_matmul(self, op):
        self._checkattr('__matmul__')
        if not hasattr(op, 'shape'):
            raise TypeError('op must have a shape attribute')
        if self.shape[1] != op.shape[0]:
            raise ValueError('dimensions must agree')

    def __matmul__(self, op):
        """
        Returns self @ op as a sparse matrix / dense array or as a LazyLinearOp.

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
            res = self.eval() @ op
        else:
            res = LazyLinearOp(init_lambda=lambda:
                                 self.eval() @ LazyLinearOp._eval_if_lazy(op),
                                 shape=(self.shape[0], op.shape[1]),
                                 root_obj=self._root_obj)
        return res

    def dot(self, op):
        """
        Alias of LazyLinearOp.__matmul__.
        """
        return self.__matmul__(op)

    def matvec(self, op):
        """
        Returns the LazyLinearOp for the multiplication self.matvec(op) or the np.ndarray
        resulting of the evaluation of self.matvec(op) if op is a np.ndarray.

        Args:
            op: must be a vector (row or column).

        <b>See also</b>: LazyLinearOp.__matmul__,<a
        href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.LinearOperator.matvec.html#scipy.sparse.linalg.LinearOperator.matvec">
        scipy.linalg.LinearOperator.matvec</a>
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
        Returns the LazyLinearOp for self^H @ v, where self^H is the conjugate transpose of A.
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
        Returns the LazyLinearOp for self^H @ v, where self^H is the conjugate transpose of A.
        """
        # LinearOperator need.
        return self.T.conj() @ op

    def __imatmul__(self, op):
        """
        Not implemented self @= op.
        """
        raise NotImplementedError(LazyLinearOp.__name__+".__imatmul__")
# can't do as follows, it recurses indefinitely because of self.eval
#        self._checkattr('__imatmul__')
#        self = LazyLinearOp(init_lambda=lambda:
#                              (self.eval()).\
#                              __imatmul__(LazyLinearOp._eval_if_lazy(op)),
#                              shape=(tuple(self.shape)),
#                              root_obj=self._root_obj)
#        return self

    def __rmatmul__(self, op):
        """
        Returns the LazyLinearOp for op @ self.

        Args:
            op: an object compatible with self for this binary operation.

        """
        self._checkattr('__rmatmul__')
        if not hasattr(op, 'shape'):
            raise TypeError('op must have a shape attribute')
        if self.shape[0] != op.shape[1]:
            raise ValueError('dimensions must agree')
        if isinstance(op, LazyLinearOp):
            res = LazyLinearOp(init_lambda=lambda:
                                 op.eval() @ self.eval(),
                                 shape=(self.shape[0], op.shape[1]),
                                 root_obj=self._root_obj)

        else:
            res = op @ self.eval()
        return res

    def __mul__(self, op):
        """
        Returns the LazyLinearOp for self * op.

        Args:
            op: an object compatible with self for this binary operation.

        """
        self._checkattr('__mul__')
        if isinstance(op, (float, int, complex)) or \
           op.ndim == 1 and op.size == self.shape[1] or \
           self.shape == op.shape or \
           op.shape[0] == 1 and op.shape[1] == self.shape[1] or \
           op.shape[1] == 1 and op.shape[0] == self.shape[0]:
            new_op = LazyLinearOp(init_lambda=lambda:
                                    self.eval() * LazyLinearOp._eval_if_lazy(op),
                                    shape=(tuple(self.shape)),
                                    root_obj=self._root_obj)
            return new_op
        else:
            raise ValueError('operands could not be broadcast together')

    def __rmul__(self, op):
        """
        Returns the LazyLinearOp for op * self.

        Args:
            op: an object compatible with self for this binary operation.

        """
        if isinstance(op, (float, int, complex)) or \
           op.ndim == 1 and op.size == self.shape[1] or \
           self.shape == op.shape or \
           op.shape[0] == 1 and op.shape[1] == self.shape[1] or \
           op.shape[1] == 1 and op.shape[0] == self.shape[0]:
            self._checkattr('__rmul__')
            new_op = LazyLinearOp(init_lambda=lambda:
                                    LazyLinearOp._eval_if_lazy(op) *
                                    self.eval(),
                                    shape=(tuple(self.shape)),
                                    root_obj=self._root_obj)
            return new_op
        else:
            raise ValueError('operands could not be broadcast together')

    def __imul__(self, op):
        """
        Not implemented self *= op.
        """
        raise NotImplementedError(LazyLinearOp.__name__+".__imul__")
#        # can't do as follows, it recurses indefinitely because of self.eval
#        self._checkattr('__imul__')
#        self = LazyLinearOp(init_lambda=lambda:
#                              (self.eval()).\
#                              __imul__(LazyLinearOp._eval_if_lazy(op)),
#                              shape=(tuple(self.shape)),
#                              root_obj=self._root_obj)
#        return self

    def toarray(self):
        """
        Returns the numpy array resulting from self evaluation.
        """
        ev_op = self.eval()
        if isinstance(ev_op, np.ndarray):
            return ev_op
        else:
            self._checkattr('toarray')
            return self.eval().toarray()

    def __getitem__(self, indices):
        """
        Returns the LazyLinearOp for indexing.

        Args:
            indices: array of length 1 or 2 which elements must be slice, integer or
            Ellipsis (...). Note that using Ellipsis for more than two indices is forbidden.

        """
        self._checkattr('__getitem__')
        if isinstance(indices, tuple) and len(indices) == 2 and isinstance(indices[0], int) and isinstance(indices[1], int):
            return self.eval().__getitem__(indices)
        else:
            return LazyLinearOp(init_lambda=lambda:
                                  (self.eval()).\
                                  __getitem__(indices),
                                  shape=self._newshape_getitem(indices),
                                  root_obj=self._root_obj)

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
        if axis == 0:
            for op in ops:
                nrows += op.shape[0]
        elif axis == 1:
            for op in ops:
                ncols += op.shape[1]
        new_shape = (nrows, ncols)
        new_op = LazyLinearOp(init_lambda=lambda:
                                cat((self.eval(),
                                     *[LazyLinearOp._eval_if_lazy(op) for op in
                                      ops]), axis=axis),
                                shape=(new_shape),
                                root_obj=self._root_obj)
        return new_op


    @property
    def real(self):
        """
        Returns the LazyLinearOp for real.
        """
        self._checkattr('real')
        new_op = LazyLinearOp(init_lambda=lambda:
                                (self.eval()).real,
                                shape=self.shape,
                                root_obj=self._root_obj)
        return new_op

    @property
    def imag(self):
        """
        Returns the LazyLinearOp for imag.
        """
        self._checkattr('imag')
        new_op = LazyLinearOp(init_lambda=lambda:
                                (self.eval()).imag,
                                shape=self.shape,
                                root_obj=self._root_obj)
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
               isLazyLinearOp(inputs[1]):
                return inputs[1].__rmatmul__(inputs[0])
            elif str(ufunc) == "<ufunc 'multiply'>" and len(inputs) >= 2 and \
               isLazyLinearOp(inputs[1]):
                return inputs[1].__rmul__(inputs[0])
            elif str(ufunc) == "<ufunc 'add'>" and len(inputs) >= 2 and \
                    isLazyLinearOp(inputs[1]):
                return inputs[1].__radd__(inputs[0])
            elif str(ufunc) == "<ufunc 'subtract'>" and len(inputs) >= 2 and \
                    isLazyLinearOp(inputs[1]):
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

def isLazyLinearOp(obj):
    """
    Returns True if obj is a LazyLinearOp, False otherwise.
    """
    return LazyLinearOp.isLazyLinearOp(obj)

def asLazyLinearOp(obj):
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
        >>> lM = asLazyLinearOp(M)
        >>> twolM = lM + lM
        >>> twolM
        <pyfaust.lazylinop.LazyLinearOp at 0x7fcd7d7750f0>
        >>> import pyfaust as pf
        >>> F = pf.rand(10, 12)
        >>> lF = asLazyLinearOp(F)
        >>> twolF = lF + lF
        >>> twolF
        <pyfaust.lazylinop.LazyLinearOp at 0x7fcd7d774730>


    <b>See also:</b> pyfaust.rand.
    """
    return LazyLinearOp2.create_from_op(obj)

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
    return LazyLinearOp2.create_from_op(obj)

def hstack(tup):
    """
    Concatenates lop1 and obj horizontally.

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
    Concatenates lop1 and obj vertically.

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

    <b>See also:</b> LazyLinearOpKron, numpy.kron.
    """
    return LazyLinearOpKron(A, B)

class LazyLinearOpKron(LazyLinearOp):
    """
    This class implements a specialization of LazyLinearOp dedicated to the
    Kronecker product of two linear operators.

    <b>See also:</b> pyfaust.lazylinop.kron.
    """

    def __init__(self, A, B):
        """
        Constructor for the A x B Kronecker product.
        """
        self.A = A
        self.B = B
        shape = (A.shape[0] * B.shape[0], A.shape[1] * B.shape[1])
        super(LazyLinearOpKron, self).__init__(lambda: A, shape, A)

    def conj(self):
        """
        Returns the LazyLinearOpKron conjugate.
        """
        return LazyLinearOpKron(asLazyLinearOp(self.A).conj(),
                                asLazyLinearOp(self.B).conj())

    def transpose(self):
        """
        Returns the LazyLinearOpKron transpose.
        """
        return LazyLinearOpKron(asLazyLinearOp(self.A).T, asLazyLinearOp(self.B).T)

    def getH(self):
        """
        Returns the LazyLinearOpKron adjoint/transconjugate.
        """
        return LazyLinearOpKron(asLazyLinearOp(self.A).getH(), asLazyLinearOp(self.B).getH())

    def __matmul__(self, op):
        """
        Returns the LazyLinearOpKron for the multiplication self @ op or if op is a np.ndarray it returns the np.ndarray (self @ op).

        Note: this specialization is particularly optimized for multiplying the
        operator by a vector.

        Args:
            op: an object compatible with self for this binary operation.

        <b>See also:</b> pyfaust.lazylinop.kron.
        """
        from threading import Thread
        from multiprocessing import cpu_count
        from os import environ
        from pyfaust import isFaust
        self._sanitize_matmul(op)
        force_eval = False # find no case where it is worth to be True
        if 'KRON_MATMUL_FORCE_EVAL' in environ:
            force_eval = environ['KRON_MATMUL_FORCE_EVAL'] == '1'
        if hasattr(op, 'reshape') and hasattr(op, '__matmul__') and hasattr(op,
                                                                            '__getitem__'):
            if force_eval:
                res = self.eval() @ op
            else:
                if isFaust(self.A) or isFaust(self.B):
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
                res = np.empty((self.shape[0], op.shape[1]))
                def out_col(j, ncols):
                    for j in range(j, min(j + ncols, op.shape[1])):
                        op_mat = op[:, j].reshape((self.A.shape[1], self.B.shape[1]))
                        res[:, j] = (LazyLinearOp._eval_if_lazy(self.A) @ op_mat @
                                  LazyLinearOp._eval_if_lazy(self.B).T).reshape(self.shape[0])
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
            res = LazyLinearOp(init_lambda=lambda:
                               self.eval() @ LazyLinearOp._eval_if_lazy(op),
                               shape=(self.shape[0], op.shape[1]),
                               root_obj=self._root_obj)
        return res

    def eval(self):
        A = self.A
        B = self.B
        if not isinstance(A, np.ndarray):
            A = A.toarray()
        if not isinstance(B, np.ndarray):
            B = B.toarray()
        return np.kron(A, B)

class LazyLinearOp2(LinearOperator):
    """
    This class implements a lazy linear operator. A LazyLinearOp is a
    specialization of a <a
    href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.LinearOperator.html">scipy.linalg.LinearOperator</a>.

    The evaluation of any defined operation is delayed until a multiplication
    by a dense matrix/vector, a call of LazyLinearOp.toarray or an explicit
    evaluation through LazyLinearOp.eval.

    To instantiate a LazyLinearOp look at pyfaust.lazylinop.asLazyLinearOp.

    Warning: This code is in a beta status.
    """
    def __init__(self, lambdas, shape, root_obj, dtype=None):
        """
        Constructor. Not meant to be used directly.

        Args:
            init_lambda: starting operation.
            shape: the initial shape of the operator.
            root_obj: the initial object the operator is based on.

        <b>See also:</b> LazyLinearOp.create, pyfaust.lazylinop.asLazyLinearOp.
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
        self._root_obj = root_obj #TODO: delete because it can't always be
        # defined (create_from_funcs and hybrid operations)
        self.dtype = dtype
        super(LazyLinearOp2, self).__init__(self.dtype, self.shape)

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
        Alias of pyfaust.lazylinop.asLazyLinearOp.

        Args:
            obj: cf. pyfaust.lazylinop.asLazyLinearOp.

        Returns:
            cf. pyfaust.lazylinop.asLazyLinearOp.

        Example:
            >>> from pyfaust.lazylinop import LazyLinearOp
            >>> import numpy as np
            >>> M = np.random.rand(10, 12)
            >>> lM = LazyLinearOp.create(M)
            >>> twolM = lM + lM
            >>> twolM
            <pyfaust.lazylinop.LazyLinearOp at 0x7fcd7d7750f0>
            >>> import pyfaust as pf
            >>> F = pf.rand(10, 12)
            >>> lF = LazyLinearOp.create(F)
            >>> twolF = lF + lF
            >>> twolF
            <pyfaust.lazylinop.LazyLinearOp at 0x7fcd7d774730>


        <b>See also:</b> pyfaust.rand.
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
            l['slice'] = None
        lop = LazyLinearOp2(lambdas, obj.shape, obj, dtype=obj.dtype)
        lopT = LazyLinearOp2(lambdasT, (obj.shape[1], obj.shape[0]), obj, dtype=obj.dtype)
        lopH = LazyLinearOp2(lambdasH, (obj.shape[1], obj.shape[0]), obj, dtype=obj.dtype)
        lopC = LazyLinearOp2(lambdasC, (obj.shape[0], obj.shape[1]), obj, dtype=obj.dtype)

        # TODO: refactor with create_from_funcs (in ctor)
        lambdas['T'] = lambda: lopT
        lambdas['H'] = lambda: lopH
        lambdas['slice'] = lambda indices: LazyLinearOp2._index_lambda(lop,
                                                                       indices)()
        lambdasT['T'] = lambda: lop
        lambdasT['H'] = lambda: lopC
        lambdasT['slice'] = lambda indices: LazyLinearOp2._index_lambda(lopT,
                                                                        indices)()
        lambdasH['T'] = lambda: lopC
        lambdasH['H'] = lambda: lop
        lambdasH['slice'] = lambda indices: LazyLinearOp2._index_lambda(lopH,
                                                                        indices)()
        return lop

    @staticmethod
    def create_from_scalar(s, shape=None):
        if shape is None:
            shape = (1, 1)
        a = np.ones(shape) * s
        return create_from_op(a)

    @staticmethod
    def create_from_funcs(matmat, rmatmat, shape, dtype=None):
        from scipy.sparse import eye as seye

        #MX = lambda X: matmat(np.eye(shape[1])) @ X
        MX = lambda X: matmat(X)
        MTX = lambda X: rmatmat(X.T).T
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

        lop = LazyLinearOp2(lambdas, shape, None)
        lopT = LazyLinearOp2(lambdasT, (shape[1], shape[0]), None)
        lopH = LazyLinearOp2(lambdasH, (shape[1], shape[0]), None)
        lopC = LazyLinearOp2(lambdasC, shape, None)

        lambdas['T'] = lambda: lopT
        lambdas['H'] = lambda: lopH
        lambdas['slice'] = lambda indices: LazyLinearOp2._index_lambda(lop,
                                                                       indices)()
        lambdasT['T'] = lambda: lop
        lambdasT['H'] = lambda: lopC
        lambdasT['slice'] = lambda indices: LazyLinearOp2._index_lambda(lopT,
                                                                        indices)()
        lambdasH['T'] = lambda: lopC
        lambdasH['H'] = lambda: lop
        lambdasH['slice'] = lambda indices: LazyLinearOp2._index_lambda(lopH,
                                                                        indices)()
        return lop

    def _checkattr(self, attr):
        if self._root_obj is not None and not hasattr(self._root_obj, attr):
            raise TypeError(attr+' is not supported by the root object of this'
                            ' LazyLinearOp')

    def _index_lambda(lop, indices):
        from scipy.sparse import eye as seye
        s = lambda: \
                LazyLinearOp2.create_from_op(seye(lop.shape[0],
                                                  format='csr')[indices[0]]) \
                @ lop @ LazyLinearOp2.create_from_op(seye(lop.shape[1], format='csr')[:, indices[1]])
        return s

    @property
    def ndim(self):
        """
        Returns the number of dimensions of the LazyLinearOp (it always 2).
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
        return self.T.H

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
        if not LazyLinearOp2.isLazyLinearOp(op):
            op = LazyLinearOp2.create_from_op(op)
        lambdas = {'@': lambda o: self @ o + op @ o,
                   'H': lambda: self.H + op.H,
                   'T': lambda: self.T + op.T,
                   'slice': lambda indices: self._slice(indices) + op._slice(indices)
                  }
        new_op = LazyLinearOp2(lambdas=lambdas,
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
#        self = LazyLinearOp2(init_lambda=lambda:
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
        if not LazyLinearOp2.isLazyLinearOp(op):
            op = LazyLinearOp2.create_from_op(op)
        lambdas = {'@': lambda o: self @ o - op @ o,
                   'H': lambda: self.H - op.H,
                   'T': lambda: self.T - op.T,
                   'slice': lambda indices: self._slice(indices) - op._slice(indices)
                  }
        new_op = LazyLinearOp2(lambdas=lambdas,
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
        if not LazyLinearOp2.isLazyLinearOp(op):
            op = LazyLinearOp2.create_from_op(op)
        lambdas = {'@': lambda o: op @ o - self @ o,
                   'H': lambda: op.H - self.H,
                   'T': lambda: op.T - self.T,
                   'slice': lambda indices: op._slice(indices) - self._slice(indices)
                  }
        new_op = LazyLinearOp2(lambdas=lambdas,
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
#        self = LazyLinearOp2(init_lambda=lambda:
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
#        self = LazyLinearOp2(init_lambda=lambda:
#                              (self.eval()).\
#                              __itruediv__(LazyLinearOp._eval_if_lazy(op)),
#                              shape=(tuple(self.shape)),
#                              root_obj=self._root_obj)
#        return self

    def _sanitize_matmul(self, op, swap=False):
        self._checkattr('__matmul__')
        if not hasattr(op, 'shape'):
            raise TypeError('op must have a shape attribute')
        if swap and op.shape[1] != self.shape[0] or not swap and self.shape[1] != op.shape[0]:
            raise ValueError('dimensions must agree')

    def __matmul__(self, op):
        """
        Returns self @ op as a sparse matrix / dense array or as a LazyLinearOp.

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
            if not LazyLinearOp2.isLazyLinearOp(op):
                op = LazyLinearOp2.create_from_op(op)
            lambdas = {'@': lambda o: self @ (op @ o),
                       'H': lambda: op.H @ self.H,
                       'T': lambda: op.T @ self.T,
                       'slice': lambda indices: self._slice((indices[0],
                                                            slice(0,
                                                                  self.shape[1])))\
                       @ op._slice((slice(0, op.shape[0]), indices[1]))
                      }
            res = LazyLinearOp2(lambdas=lambdas,
                               shape=(self.shape[0], op.shape[1]),
                               root_obj=None)
#            res = LazyLinearOp2.create_from_op(super(LazyLinearOp2,
#                                                     self).__matmul__(op))
        return res

    def dot(self, op):
        """
        Alias of LazyLinearOp.__matmul__.
        """
        return self.__matmul__(op)

    def matvec(self, op):
        """
        Returns the LazyLinearOp for the multiplication self.matvec(op) or the np.ndarray
        resulting of the evaluation of self.matvec(op) if op is a np.ndarray.

        Args:
            op: must be a vector (row or column).

        <b>See also</b>: LazyLinearOp.__matmul__,<a
        href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.LinearOperator.matvec.html#scipy.sparse.linalg.LinearOperator.matvec">
        scipy.linalg.LinearOperator.matvec</a>
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
        Returns the LazyLinearOp for self^H @ v, where self^H is the conjugate transpose of A.
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
        Returns the LazyLinearOp for self^H @ v, where self^H is the conjugate transpose of A.
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
        Returns the LazyLinearOp for op @ self.

        Args:
            op: an object compatible with self for this binary operation.

        """
        self._checkattr('__rmatmul__')
        from scipy.sparse import issparse
        self._sanitize_matmul(op, swap=True)
        if isinstance(op, np.ndarray) or issparse(op):
            res = op @ self.toarray()
        else:
            if not LazyLinearOp2.isLazyLinearOp(op):
                op = LazyLinearOp2.create_from_op(op)
            lambdas = {'@': lambda o: (op @ self) @ o,
                       'H': lambda: self.H @ op.H,
                       'T': lambda: self.T @ op.T,
                       'slice': lambda indices: (op @ self)._slice(indices)
                      }
            res = LazyLinearOp2(lambdas=lambdas,
                               shape=(op.shape[0], self.shape[1]),
                               root_obj=None)
        return res

    def __mul__(self, other):
        """
        Returns the LazyLinearOp for self * other.

        Args:
            other: a scalar or a vector.

        """
        self._checkattr('__mul__')
        from scipy.sparse import eye, diags
        if np.isscalar(other):
            S = eye(self.shape[1], format='csr') * other
            lop = LazyLinearOp2.create_from_op(S)
        elif hasattr(other, 'ndim') and other.ndim == 1 or other.ndim == 2 and other.shape[0] == 1:
            if other.size == 1:
                return self * self.item()
            else:
                D = diags(other.squeeze())
                lop = LazyLinearOp2.create_from_op(D)
        else:
            raise TypeError('Unsupported type for other')
        new_op = self @ lop
        return new_op

    def __rmul__(self, s):
        """
        Returns the LazyLinearOp for op * self.

        Args:
            s: a scalar.

        """
        # because s is a scalar, it is commutative
        # for vector broadcasting too
        return self * s


    def __imul__(self, op):
        """
        Not implemented self *= op.
        """
        raise NotImplementedError(LazyLinearOp.__name__+".__imul__")

    def toarray(self):
        """
        Returns the numpy array resulting from self evaluation.
        """
        #from scipy.sparse import eye
        #return self @ eye(self.shape[1], self.shape[1], format='csr')
        # don't use csr because of function based LazyLinearOp2
        # (e.g. scipy fft receives only numpy array)
        return self @ np.eye(self.shape[1])

    def __getitem__(self, indices):
        """
        Returns the LazyLinearOp for slicing/indexing.

        Args:
            indices: array of length 1 or 2 which elements must be slice, integer or
            Ellipsis (...). Note that using Ellipsis for more than two indices is forbidden.

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
        if self.shape[1] != op.shape[1]:
            raise ValueError('The number of columns of self and op are not the'
                             ' same')
        if not LazyLinearOp2.isLazyLinearOp(op):
            op = LazyLinearOp2.create_from_op(op)
        lambdas = {'@': lambda o: self._vstack_mul_lambda(op, o)(),
                   'H': lambda: self.H.hstack(op.H),
                   'T': lambda: self.T.hstack(op.T),
                   'slice': lambda indices: self._vstack_slice(op, indices)()
                  }
        new_op = LazyLinearOp2(lambdas=lambdas,
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
        if self.shape[0] != op.shape[0]:
            raise ValueError('The number of columns of self and op are not the'
                             ' same')
        if not LazyLinearOp2.isLazyLinearOp(op):
            op = LazyLinearOp2.create_from_op(op)
        lambdas = {'@': lambda o: self._hstack_mul_lambda(op, o)(),
               'H': lambda: self.H.vstack(op.H),
               'T': lambda: self.T.vstack(op.T),
               'slice': lambda indices: self._hstack_slice(op, indices)()
              }
        new_op = LazyLinearOp2(lambdas=lambdas,
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
        new_op = LazyLinearOp2(lambdas=lambdas,
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
                   or issparse(o) else real(self @ o),
                   'H': lambda: self.T.imag,
                   'T': lambda: self.T.imag,
                   'slice': lambda indices: self._slice(indices).imag
                  }
        new_op = LazyLinearOp2(lambdas=lambdas,
                           shape=tuple(self.shape),
                           root_obj=None)
        return new_op


    @staticmethod
    def isLazyLinearOp(obj):
        """
        Returns True if obj is a LazyLinearOp, False otherwise.
        """
        return isinstance(obj, LazyLinearOp2)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method == '__call__':
            if str(ufunc) == "<ufunc 'matmul'>" and len(inputs) >= 2 and \
               LazyLinearOp2.isLazyLinearOp(inputs[1]):
                return inputs[1].__rmatmul__(inputs[0])
            elif str(ufunc) == "<ufunc 'multiply'>" and len(inputs) >= 2 and \
               LazyLinearOp2.isLazyLinearOp(inputs[1]):
                return inputs[1].__rmul__(inputs[0])
            elif str(ufunc) == "<ufunc 'add'>" and len(inputs) >= 2 and \
                    LazyLinearOp2.isLazyLinearOp(inputs[1]):
                return inputs[1].__radd__(inputs[0])
            elif str(ufunc) == "<ufunc 'subtract'>" and len(inputs) >= 2 and \
                    LazyLinearOp2.isLazyLinearOp(inputs[1]):
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
    Returns a LazyLinearOp defined by shape and at least matvec.

    NOTE: At least a matvec or a matmat function must be passed in kwargs.

    Args:
        shape: (tuple) dimensions (M, N).
        matvec: (callable) returns A * v (v a vector).
        rmatvec: (callable) returns A^H * v (v a vector of size N).
        matmat: (callable) returns A * V (V a dense matrix of dimensions (N, K)).
        rmatmat: (callable) returns A^H * V (V a dense matrix of dimensions (M, K)).
        dtype: data type of the matrix (can be None).

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
        out = np.empty((shape[0], M.shape[1]), dtype=dtype if dtype is not None
                      else M.dtype)
        for i in range(M.shape[1]):
            out[:, i] = _matvec(M[:,i])
        return out

    if matmat is None:
        matmat = lambda M: _matmat(M, matvec)

    if rmatmat is None and rmatvec is not None:
        rmatmat = lambda M: _matmat(M, rmatvec)

    return LazyLinearOp2.create_from_funcs(matmat, rmatmat, shape, dtype=dtype)
