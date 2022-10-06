## @package pyfaust.lazylinop @brief The pyfaust module for lazy linear operators.

import numpy as np
from scipy.sparse.linalg import LinearOperator

HANDLED_FUNCTIONS = {}

class LazyLinearOp(LinearOperator):
    """
    This class implements a lazy linear operator.

    The evaluation of any defined operation is delayed.

    For creation and evaluation look at LazyLinearOp.create and
    LazyLinearOp.eval.
    """
    def __init__(self, init_lambda = None, shape=(0,0), root_obj=None):
        """
        Constructor. Not meant to be used directly.

        Args:
            init_lambda: starting operation

        <b>See also:</b> LazyLinearOp.create.
        """
        self._lambda_stack = init_lambda
        self.shape = shape
        self._root_obj = root_obj
        self.dtype = None
        super(LazyLinearOp, self).__init__(self.dtype, self.shape)

    @staticmethod
    def create(obj):
        """
        Creates a Lazy linear operator.

        NOTE: obj must support operations and attributes defined in this class.

        Args:
            obj: the root object on which the LazyLinearOp is based (it could
            be a numpy array, a scipy matrix, a Faust object).

        Returns:
            a LazyLinearOp instance based on obj.

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
            >>> lF = LazyLinearOp.create(M)
            >>> twolF = lF + lF
            >>> twolF
            <pyfaust.lazylinop.LazyLinearOp at 0x7fcd7d774730>
        """
        return LazyLinearOp(lambda:obj, obj.shape, obj)

    def eval(self):
        """
        Evaluate the LazyLinearOp. All stacked virtual operations are then
        evaluated.

        Returns:
            a "real" linear operator object (which the type depends of the
            LazyLinearOp initialization through LazyLinearOp.create and the
            operations made upon this object).

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
        return o.eval() if isinstance(o, LazyLinearOp) else o

    @property
    def ndim(self):
        """
        Returns the number of dimensions of the LazyLinearOp.
        """
        return 2

    def transpose(self):
        """
        Returns the LazyLinearOp transpose.
        """
        self._checkattr('transpose')
        new_op = self.__class__(init_lambda=lambda:
                                (self._lambda_stack()).transpose(),
                                shape=self.shape,
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
        new_op = self.__class__(init_lambda=lambda:
                                (self._lambda_stack()).conj(),
                                shape=(tuple(self.shape)),
                                root_obj=self._root_obj)
        return new_op

    def conjugate(self):
        """
        Returns the LazyLinearOp conjugate.
        """
        return self.conj()

    def getH(self):
        """
        Returns the LazyLinearOp adjoint.
        """
        self._checkattr('getH')
        new_op = self.__class__(init_lambda=lambda:
                                (self._lambda_stack()).getH(),
                                shape=(tuple(self.shape)),
                                root_obj=self._root_obj)
        return new_op

    @property
    def H(self):
        """
        The LazyLinearOp adjoint.
        """
        return self.getH()

    def _adjoint(self):
        """
        Returns the LazyLinearOp adjoint.
        """
        return self.H

    def __add__(self, op):
        """
        Returns the LazyLinearOp for self + op.
        """
        self._checkattr('__add__')
        new_op = self.__class__(init_lambda=lambda:
                                self._lambda_stack() + LazyLinearOp._eval_if_lazy(op),
                                shape=(tuple(self.shape)),
                                root_obj=self._root_obj)
        return new_op

    def __radd__(self, op):
        """
        Returns the LazyLinearOp for op + self.
        """
        return self.__add__(op)

    def __iadd__(self, op):
        """
        Not Implemented self += op.
        """
        raise NotImplementedError(self.__class__.__name__+".__iadd__")
# can't do as follows, it recurses indefinitely because on self.eval
#        self._checkattr('__iadd__')
#        self = self.__class__(init_lambda=lambda:
#                              (self._lambda_stack()).\
#                              __iadd__(LazyLinearOp._eval_if_lazy(op)),
#                              shape=(tuple(self.shape)),
#                              root_obj=self._root_obj)
#        return self


    def __sub__(self, op):
        """
        Returns the LazyLinearOp for self - op.
        """
        self._checkattr('__sub__')
        new_op = self.__class__(init_lambda=lambda: self._lambda_stack() - LazyLinearOp._eval_if_lazy(op),
                                shape=(tuple(self.shape)),
                                root_obj=self._root_obj)
        return new_op

    def __rsub__(self, op):
        """
        Returns the LazyLinearOp for op - self.
        """
        self._checkattr('__rsub__')
        new_op = self.__class__(init_lambda=lambda:
                                LazyLinearOp._eval_if_lazy(op) -
                                self._lambda_stack(),
                                shape=(tuple(self.shape)),
                                root_obj=self._root_obj)
        return new_op

    def __isub__(self, op):
        """
        Not implemented self -= op.
        """
        raise NotImplementedError(self.__class__.__name__+".__isub__")
# can't do as follows, it recurses indefinitely because on self.eval
#        self._checkattr('__isub__')
#        self = self.__class__(init_lambda=lambda:
#                              (self._lambda_stack()).\
#                              __isub__(LazyLinearOp._eval_if_lazy(op)),
#                              shape=(tuple(self.shape)),
#                              root_obj=self._root_obj)
#        return self


    def __truediv__(self, op):
        """
        Returns the LazyLinearOp for self / op.
        """
        self._checkattr('__truediv__')
        new_op = self.__class__(init_lambda=lambda:
                                self._lambda_stack() / LazyLinearOp._eval_if_lazy(op),
                                shape=(tuple(self.shape)),
                                root_obj=self._root_obj)
        return new_op

    def __itruediv__(self, op):
        """
        Not implemented self /= op.
        """
        raise NotImplementedError(self.__class__.__name__+".__itruediv__")
# can't do as follows, it recurses indefinitely because on self.eval
#
#        self._checkattr('__itruediv__')
#        self = self.__class__(init_lambda=lambda:
#                              (self._lambda_stack()).\
#                              __itruediv__(LazyLinearOp._eval_if_lazy(op)),
#                              shape=(tuple(self.shape)),
#                              root_obj=self._root_obj)
#        return self

    def __matmul__(self, op):
        """
        Returns the LazyLinearOp for self @ op or the np.ndarray resulting
        of the evaluation of self @ op if op is a np.ndarray.
        """
        self._checkattr('__matmul__')
        if not hasattr(op, 'shape'):
            raise TypeError('op must have a shape attribute')
        if self.shape[1] != op.shape[0]:
            raise ValueError('dimensions must agree')
        if isinstance(op, LazyLinearOp):
            res = self.__class__(init_lambda=lambda:
                                 self.eval() @ op.eval(),
                                 shape=(self.shape[0], op.shape[1]),
                                 root_obj=self._root_obj)
        else:
            res = self.eval() @ op
        return res

    def dot(self, op):
        """
        Returns the LazyLinearOp for self.dot(op) or the np.ndarray resulting
        of the evaluation of self @ op if op is a np.ndarray.
        """
        return self.__matmul__(op)

    def matvec(self, op):
        """
        Returns the LazyLinearOp for self.matvec(op) or the np.ndarray resulting
        of the evaluation of self @ op if op is a np.ndarray.

        Args:
            op: must be a vector (row or column).
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
        Returns self^H @ v, where self^H is the conjugate transpose of A.
        LinearOperator need.
        """
        return self.H @ op

    def _matmat(self, op):
        """
        See matmat.
        LinearOperator need.
        """
        if not hasattr(op, 'shape') or not hasattr(op, 'ndim'):
            raise TypeError('op must have shape and ndim attributes')
        if op.ndim > 2 or op.ndim == 0:
            raise ValueError('op.ndim must be 1 or 2')
        return self.__matmul__(op)

    def _rmatmat(self, op):
        """
        Returns self^H @ v, where self^H is the conjugate transpose of A.
        LinearOperator need.
        """
        return self.H @ op

    def __imatmul__(self, op):
        """
        Not implemented self @= op.
        """
        raise NotImplementedError(self.__class__.__name__+".__imatmul__")
# can't do as follows, it recurses indefinitely because on self.eval
#        self._checkattr('__imatmul__')
#        self = self.__class__(init_lambda=lambda:
#                              (self._lambda_stack()).\
#                              __imatmul__(LazyLinearOp._eval_if_lazy(op)),
#                              shape=(tuple(self.shape)),
#                              root_obj=self._root_obj)
#        return self

    def __rmatmul__(self, op):
        """
        Returns the LazyLinearOp for op @ self.
        """
        self._checkattr('__rmatmul__')
        if not hasattr(op, 'shape'):
            raise TypeError('op must have a shape attribute')
        if self.shape[0] != op.shape[1]:
            raise ValueError('dimensions must agree')
        if isinstance(op, LazyLinearOp):
            res = self.__class__(init_lambda=lambda:
                                 op.eval() @ self.eval(),
                                 shape=(self.shape[0], op.shape[1]),
                                 root_obj=self._root_obj)

        else:
            res = op @ self.eval()
        return res

    def __mul__(self, op):
        """
        Returns the LazyLinearOp for self * op.
        """
        self._checkattr('__mul__')
        if isinstance(op, (float, int, complex)) or \
           op.ndim == 1 and op.size == self.shape[1] or \
           self.shape == op.shape or \
           op.shape[0] == 1 and op.shape[1] == self.shape[1] or \
           op.shape[1] == 1 and op.shape[0] == self.shape[0]:
            new_op = self.__class__(init_lambda=lambda:
                                    self._lambda_stack() * LazyLinearOp._eval_if_lazy(op),
                                    shape=(tuple(self.shape)),
                                    root_obj=self._root_obj)
            return new_op
        else:
            raise ValueError('operands could not be broadcast together')

    def __rmul__(self, op):
        """
        Returns the LazyLinearOp for op * self.
        """
        if isinstance(op, (float, int, complex)) or \
           op.ndim == 1 and op.size == self.shape[1] or \
           self.shape == op.shape or \
           op.shape[0] == 1 and op.shape[1] == self.shape[1] or \
           op.shape[1] == 1 and op.shape[0] == self.shape[0]:
            self._checkattr('__rmul__')
            new_op = self.__class__(init_lambda=lambda:
                                    LazyLinearOp._eval_if_lazy(op) *
                                    self._lambda_stack(),
                                    shape=(tuple(self.shape)),
                                    root_obj=self._root_obj)
            return new_op
        else:
            raise ValueError('operands could not be broadcast together')

    def __imul__(self, op):
        """
        Not implemented self *= op.
        """
        raise NotImplementedError(self.__class__.__name__+".__imul__")
#        # can't do as follows, it recurses indefinitely because on self.eval
#        self._checkattr('__imul__')
#        self = self.__class__(init_lambda=lambda:
#                              (self._lambda_stack()).\
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
        """
        self._checkattr('__getitem__')
        if isinstance(indices, tuple) and len(indices) == 2 and isinstance(indices[0], int) and isinstance(indices[1], int):
            return self.eval().__getitem__(indices)
        else:
            return self.__class__(init_lambda=lambda:
                                  (self._lambda_stack()).\
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
            #self[i:j] a group of contiguous lines
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

    def concatenate(self, op, axis=0):
        """
        Returns the LazyLinearOp for the concatenation of self and op.

        Args:
            axis: axis of concatenation (0 for rows, 1 for columns).
        """
        from pyfaust import concatenate as cat
        new_shape = (self.shape[0] + op.shape[0] if axis == 0 else self.shape[0],
                     self.shape[1] + op.shape[1] if axis == 1 else self.shape[1])
        new_op = self.__class__(init_lambda=lambda:
                                cat((self._lambda_stack(),
                                     LazyLinearOp._eval_if_lazy(op)), axis=axis),
                                shape=(new_shape),
                                root_obj=self._root_obj)
        return new_op


    @property
    def real(self):
        """
        Returns the LazyLinearOp for real.
        """
        self._checkattr('real')
        new_op = self.__class__(init_lambda=lambda:
                                (self._lambda_stack()).real,
                                shape=self.shape,
                                root_obj=self._root_obj)
        return new_op

    @property
    def imag(self):
        """
        Returns the LazyLinearOp for imag.
        """
        self._checkattr('imag')
        new_op = self.__class__(init_lambda=lambda:
                                (self._lambda_stack()).imag,
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
            N = None
            fausts = []
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
        print("__array_function__")
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

    Args:
        obj: it can be a ndarray (ndim == 2), a Faust or a scipy matrix.
    """
    return LazyLinearOp.create(obj)

def hstack(lop1, obj):
    """
    Concatenates lop1 and obj horizontally.

    Args:
        lop1: a LazyLinearOp object.
        obj: any array / matrix / LazyLinearOp compatible in dimensions.
    """
    if isLazyLinearOp(lop1):
        return lop1.concatenate(obj, axis=1)
    else:
        raise TypeError('lop1 must be a LazyLinearOp')

def vstack(lop1, obj):
    """
    Concatenates lop1 and obj vertically.

    Args:
        lop1: a LazyLinearOp object.
        obj: any array / matrix / LazyLinearOp compatible in dimensions.
    """
    if isLazyLinearOp(lop1):
        return lop1.concatenate(obj, axis=0)
    else:
        raise TypeError('lop1 must be a LazyLinearOp')
