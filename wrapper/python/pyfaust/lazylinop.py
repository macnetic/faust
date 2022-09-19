# experimental block start
from abc import ABC, abstractmethod

class LazyLinearOp(ABC):

    @abstractmethod
    def __init__(self, init_lambda = None, shape=(0,0), root_obj=None):
        self.lambda_stack = init_lambda
        self._shape = shape
        self.root_obj = root_obj

    def eval(self):
        return self.lambda_stack()

    def _checkattr(self, attr):
        if not hasattr(self.root_obj, attr):
            raise TypeError(attr+' is not supported by the root object of this'
                            ' LazyLinearOp')

    @staticmethod
    def _eval_if_lazy(o):
        return o.eval() if isinstance(o, LazyLinearOp) else o

    @property
    def shape(self):
        return self._shape

    @property
    def ndim(self):
        return 2

    def transpose(self):
        self._checkattr('transpose')
        new_op = self.__class__(init_lambda=lambda:
                                (self.lambda_stack()).transpose(),
                                shape=(self.shape[1], self.shape[0]),
                                root_obj=self.root_obj)
        return new_op

    @property
    def T(self):
        return self.transpose()

    def conj(self):
        self._checkattr('conj')
        new_op = self.__class__(init_lambda=lambda:
                                (self.lambda_stack()).conj(),
                                shape=(tuple(self.shape)),
                                root_obj=self.root_obj)
        return new_op

    def conjugate(self):
        return self.conj()

    def getH(self):
        self._checkattr('getH')
        new_op = self.__class__(init_lambda=lambda:
                                (self.lambda_stack()).getH(),
                                shape=(tuple(self.shape)),
                                root_obj=self.root_obj)
        return new_op

    @property
    def H(self):
        return self.getH()

    def __add__(self, op):
        self._checkattr('__add__')
        new_op = self.__class__(init_lambda=lambda:
                                (self.lambda_stack()).\
                                __add__(LazyLinearOp._eval_if_lazy(op)),
                                shape=(tuple(self.shape)),
                                root_obj=self.root_obj)
        return new_op

    def __radd__(self, op):
        return self.__add__(op)

    def __iadd__(self, op):
        self._checkattr('__iadd__')
        self = self.__class__(init_lambda=lambda:
                              (self.lambda_stack()).\
                              __iadd__(LazyLinearOp._eval_if_lazy(op)),
                              shape=(tuple(self.shape)),
                              root_obj=self.root_obj)
        return self


    def __sub__(self, op):
        self._checkattr('__sub__')
        new_op = self.__class__(init_lambda=lambda:
                                (self.lambda_stack()).\
                                __sub__(LazyLinearOp._eval_if_lazy(op)),
                                shape=(tuple(self.shape)),
                                root_obj=self.root_obj)
        return new_op

    def __rsub__(self, op):
        self._checkattr('__rsub__')
        new_op = self.__class__(init_lambda=lambda:
                                (self.lambda_stack()).\
                                __rsub__(LazyLinearOp._eval_if_lazy(op)),
                                shape=(tuple(self.shape)),
                                root_obj=self.root_obj)
        return new_op

    def __isub__(self, op):
        self._checkattr('__isub__')
        self = self.__class__(init_lambda=lambda:
                              (self.lambda_stack()).\
                              __isub__(LazyLinearOp._eval_if_lazy(op)),
                              shape=(tuple(self.shape)),
                              root_obj=self.root_obj)
        return self


    def __truediv__(self, op):
        self._checkattr('__truediv__')
        new_op = self.__class__(init_lambda=lambda:
                                (self.lambda_stack()).\
                                __truediv__(LazyLinearOp._eval_if_lazy(op)),
                                shape=(tuple(self.shape)),
                                root_obj=self.root_obj)
        return new_op

    def __itruediv__(self, op):
        self._checkattr('__itruediv__')
        self = self.__class__(init_lambda=lambda:
                              (self.lambda_stack()).\
                              __itruediv__(LazyLinearOp._eval_if_lazy(op)),
                              shape=(tuple(self.shape)),
                              root_obj=self.root_obj)
        return self

    def __matmul__(self, op):
        self._checkattr('__matmul__')
        if not hasattr(op, 'shape'):
            raise TypeError('op must have a shape attribute')
        if self.shape[1] != op.shape[0]:
            raise ValueError('dimensions must agree')
        if isinstance(op, LazyLinearOp):
            res = self.__class__(init_lambda=lambda:
                                 self.eval() @ op.eval(),
                                 shape=(self.shape[0], op.shape[1]),
                                 root_obj=self.root_obj)
        else:
            res = self.eval() @ op
        return res

    def dot(self, op):
        return self.__matmul__(op)

    def matvec(self, op):
        if not hasattr(op, 'shape') or not hasattr(op, 'ndim'):
            raise TypeError('op must have shape and ndim attributes')
        if op.ndim > 2 or op.ndim == 0:
            raise ValueError('op.ndim must be 1 or 2')
        if op.ndim != 1 and op.shape[0] != 1 and op.shape[1] != 1:
            raise ValueError('op must be a vector -- attribute ndim to 1 or'
                             ' shape[0] or shape[1] to 1')
        return self.__matmul__(op)

    def __imatmul__(self, op):
        self._checkattr('__imatmul__')
        self = self.__class__(init_lambda=lambda:
                              (self.lambda_stack()).\
                              __imatmul__(LazyLinearOp._eval_if_lazy(op)),
                              shape=(tuple(self.shape)),
                              root_obj=self.root_obj)
        return self

    def __rmatmul__(self, op):
        self._checkattr('__rmatmul__')
        if not hasattr(op, 'shape'):
            raise TypeError('op must have a shape attribute')
        if self.shape[0] != op.shape[1]:
            raise ValueError('dimensions must agree')
        if isinstance(op, LazyLinearOp):
            res = self.__class__(init_lambda=lambda:
                                 op.eval() @ self.eval(),
                                 shape=(self.shape[0], op.shape[1]),
                                 root_obj=self.root_obj)

        else:
            res = op @ self.eval()
        return res

    def __mul__(self, op):
        self._checkattr('__mul__')
        new_op = self.__class__(init_lambda=lambda:
                                (self.lambda_stack()).\
                                __mul__(LazyLinearOp._eval_if_lazy(op)),
                                shape=(tuple(self.shape)),
                                root_obj=self.root_obj)
        return new_op

    def __rmul__(self, op):
        self._checkattr('__rmul__')
        new_op = self.__class__(init_lambda=lambda:
                                (self.lambda_stack()).\
                                __rmul__(LazyLinearOp._eval_if_lazy(op)),
                                shape=(tuple(self.shape)),
                                root_obj=self.root_obj)
        return new_op

    def __imul__(self, op):
        self._checkattr('__imul__')
        self = self.__class__(init_lambda=lambda:
                              (self.lambda_stack()).\
                              __imul__(LazyLinearOp._eval_if_lazy(op)),
                              shape=(tuple(self.shape)),
                              root_obj=self.root_obj)
        return self

    def toarray(self):
        self._checkattr('toarray')
        return self.eval().toarray()

    def __getitem__(self, indices):
        self._checkattr('__getitem__')
        if isinstance(indices, tuple) and len(indices) == 2 and isinstance(indices[0], int) and isinstance(indices[1], int):
            return self.eval().__getitem__(indices)
        else:
            self = self.__class__(init_lambda=lambda:
                                  (self.lambda_stack()).\
                                  __getitem__(indices),
                                  shape=(tuple(self.shape)),
                                  root_obj=self.root_obj)
            return self

    def concatenate(self, op, axis=0):
        from pyfaust import concatenate as cat
        new_op = self.__class__(init_lambda=lambda:
                                cat((self.lambda_stack(),
                                     LazyLinearOp._eval_if_lazy(op)), axis=axis),
                                shape=(tuple(self.shape)),
                                root_obj=self.root_obj)
        return new_op


class LazyFaust(LazyLinearOp):

    def __init__(self, init_lambda, shape, root_obj):
        super(LazyFaust, self).__init__(init_lambda = init_lambda, shape=shape,
                                        root_obj=root_obj)

    @staticmethod
    def create(F):
        return LazyFaust(lambda:F, F.shape, F)
# experimental block end
