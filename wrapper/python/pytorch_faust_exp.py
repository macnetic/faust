# this comment line matters for cmake
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

