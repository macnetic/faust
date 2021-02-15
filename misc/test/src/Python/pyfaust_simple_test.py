import sys
import unittest
from pyfaust import rand as frand, Faust, vstack, hstack, isFaust
from numpy.random import randint
import numpy as np
from scipy.sparse import csr_matrix
import tempfile
import os
import random

dev = 'cpu'

class PyfaustSimpleTest(unittest.TestCase):

    MIN_NUM_FACTORS = 1
    MAX_NUM_FACTORS = 4
    MAX_DIM_SIZE = 256
    MIN_DIM_SIZE = 3

    def setUp(self):
        """
        """
        nrows = randint(PyfaustSimpleTest.MIN_DIM_SIZE,
                        PyfaustSimpleTest.MAX_DIM_SIZE+1)
        ncols = randint(PyfaustSimpleTest.MIN_DIM_SIZE,
                        PyfaustSimpleTest.MAX_DIM_SIZE+1)
        nfacts = randint(PyfaustSimpleTest.MIN_NUM_FACTORS,
                         PyfaustSimpleTest.MAX_NUM_FACTORS+1)
        self.F = frand(nrows, ncols, num_factors=nfacts, dev=dev);
        self.nrows = nrows
        self.ncols = ncols
        self.nfacts = nfacts

    def test_toarray(self):
        print("Faust.toarray")
        # this test of toarray depends (and tests) Faust.factors() and
        # numfactors()
        factors = [ self.F.factors(i) for i in range(self.F.numfactors()) ]
        fullF = np.eye(self.ncols)
        for fac in reversed(factors):
            fullF = fac @ fullF
        self.assertTrue(np.allclose(self.F.toarray(), fullF))

    def test_nbytes(self):
        print("Faust.nbytes")
        # this test of nbytes depends (and tests) Faust.factors() and
        # numfactors()
        def sparse_fac_size(sfac):
            int_size = 4
            double_size = 8
            cplx_size = 16
            nnz = sfac.nnz
            nrows = sfac.shape[0]
            ncols = sfac.shape[1]
            size = nnz*int_size+(nrows+1)*int_size
            if sfac.dtype == np.complex:
                size += cplx_size*nnz
            elif sfac.dtype == np.float:
                size += double_size*nnz
            return size
        Fsize = 0
        for i in range(0, self.F.numfactors()):
            fac = self.F.factors(i)
            if isinstance(fac, np.ndarray):
                Fsize += fac.nbytes
            elif isinstance(fac,csr_matrix):
                Fsize += sparse_fac_size(fac)
        self.assertEqual(Fsize, self.F.nbytes)


    def test_shape(self):
        print("Faust.shape")
        self.assertEqual(self.nrows, self.F.shape[0])
        self.assertEqual(self.ncols, self.F.shape[1])

    def test_ndim(self):
        print("Faust.ndim")
        self.assertEqual(self.F.ndim, 2)

    def test_size(self):
        print("Faust.size")
        self.assertEqual(self.F.size, self.nrows*self.ncols)

    def test_device(self):
        print("Faust.device")
        self.assertEqual(self.F.device, dev)

    def test_transpose(self):
        print("Faust.transpose")
        self.assertTrue(np.allclose(self.F.T.toarray().T, self.F.toarray()))
        self.assertTrue(np.allclose(self.F.transpose().toarray().T, self.F.toarray()))

    def test_conj(self):
        print("Faust.conj")
        self.assertTrue(np.allclose(self.F.conj().toarray(), self.F.toarray().conj()))
        self.assertTrue(np.allclose(self.F.conjugate().toarray(),
                                    self.F.toarray().conjugate()))

    def test_transconj(self):
        print("Faust.H")
        self.assertTrue(np.allclose(self.F.H.toarray().conj().T, self.F.toarray()))
        self.assertTrue(np.allclose(self.F.getH().toarray().conj().T, self.F.toarray()))

    def test_pruneout(self):
        print("Faust.pruneout")
        self.assertLessEqual(self.F.pruneout().nbytes, self.F.nbytes)
        self.assertTrue(np.allclose(self.F.pruneout().toarray(),
                                    self.F.toarray()))

    def test_add(self):
        print("Faust.__add__,__radd__")
        G = frand(self.nrows, self.ncols)
        self.assertTrue(np.allclose((self.F+G).toarray(),
                                    self.F.toarray()+G.toarray()))
        self.assertTrue(np.allclose((self.F+G.toarray()).toarray(),
                                    self.F.toarray()+G.toarray()))

    def test_sub(self):
        print("Faust.__sub__,__rsub__")
        G = frand(self.nrows, self.ncols)
        self.assertTrue(np.allclose((self.F-G).toarray(),
                                    self.F.toarray()-G.toarray()))
        self.assertTrue(np.allclose((self.F-G.toarray()).toarray(),
                                    self.F.toarray()-G.toarray()))

    def test_div(self):
        print("Faust.__truediv__")
        self.assertTrue(np.allclose((self.F/2).toarray(), self.F.toarray()/2))

    def test_matmul(self):
        print("Faust.__matmul__, dot, __rmatmul__")
        G = frand(self.ncols, self.nrows)
        self.assertTrue(np.allclose((self.F@G).toarray(),
                                    self.F.toarray()@G.toarray()))
        self.assertTrue(np.allclose((self.F@G.toarray()),
                                    self.F.toarray()@G.toarray()))
        self.assertTrue(np.allclose((self.F.dot(G)).toarray(),
                                    self.F.toarray().dot(G.toarray())))

    def test_concatenate(self):
        print("Faust.concatenate, pyfaust.vstack, pyfaust.hstack")
        G = frand(self.nrows, self.ncols)
        self.assertTrue(np.allclose((self.F.concatenate(G)).toarray(),
                                    np.concatenate((self.F.toarray(),
                                                    G.toarray()))))
        self.assertTrue(np.allclose((self.F.concatenate(G.toarray())).toarray(),
                                    np.concatenate((self.F.toarray(),
                                                    G.toarray()))))
        self.assertTrue(np.allclose(vstack((self.F, G.toarray())).toarray(),
                                    np.vstack((self.F.toarray(),
                                                    G.toarray()))))
        self.assertTrue(np.allclose(hstack((self.F, G.toarray())).toarray(),
                                    np.hstack((self.F.toarray(),
                                                    G.toarray()))))

        def test_isFaust(self):
            print("test pyfaust.isFaust")
            self.assertTrue(isFaust(self.F))

    def test_nnz_sum(self):
        print("Faust.nnz_sum")
        nnz_sum = 0
        for i in range(0, self.F.numfactors()):
            nnz_sum += self.F.factors(i).nnz
        self.assertEqual(nnz_sum, self.F.nnz_sum())

    def test_density(self):
        print("Faust.density")
        self.assertEqual(self.F.density(), self.F.nnz_sum()/self.F.size)

    def test_rcg(self):
        print("Faust.rcg")
        self.assertEqual(self.F.density(), self.F.nnz_sum()/self.F.size)
        self.assertEqual(self.F.rcg(), 1/self.F.density())

    def test_norm(self):
        print("Faust.norm")
        for nt in ['fro', 1, 2, np.inf]:
            print(self.F.norm(nt),
                                        np.linalg.norm(self.F.toarray(), nt))
            self.assertTrue(np.allclose(self.F.norm(nt),
                                        np.linalg.norm(self.F.toarray(), nt),
                                        rtol=1e-2))

    def test_normalize(self):
        print("Faust.normalize")
        FA = self.F.toarray()
        for nt in ['fro', 1, 2, np.inf]:
            NF = self.F.normalize(nt)
            NFA = NF.toarray()
            for j in range(NFA.shape[1]):
                n = np.linalg.norm(FA[:,j],  2
                                   if nt ==
                                   'fro' else
                                   nt,)
                self.assertTrue(n == 0 or np.allclose(FA[:,j]/n, NFA[:,j], rtol=1e-3))

    def test_numfactors(self):
        print("Faust.numfactors")
        self.assertEqual(self.nfacts, self.F.numfactors())
        self.assertEqual(self.nfacts, self.F.__len__())

    def test_factors(self):
        print("Faust.factors")
        Fc = Faust([self.F.factors(i) for i in range(self.F.numfactors())])
        self.assertTrue(np.allclose(Fc.toarray(), self.F.toarray()))
        i = randint(0,self.F.numfactors())
        j = randint(0,self.F.numfactors())
        if self.F.numfactors() > 1:
            if i > j:
                tmp = i
                i = j
                j = tmp
            elif i == j:
                if j < self.F.numfactors()-1:
                    j += 1
                elif i > 0:
                    i -= 1
                else:
                    return # irrelevant test factors(i) already tested above
            Fp = self.F.factors(range(i,j+1))
            print(Fp)
            for k in range(i,j+1):
                #self.assertTrue(np.allclose(self.F.factors(k), Fp.factors(k)))
                self._assertAlmostEqual(self.F.factors(k), Fp.factors(k-i))


    def _assertAlmostEqual(self, a, b):
        if not isinstance(a, np.ndarray):
            a = a.toarray()
        if not isinstance(b, np.ndarray):
            b = b.toarray()
        self.assertTrue(np.allclose(a, b))

    def test_left_right(self):
        print("Faust.right, Faust.left")
        i = randint(0, self.F.numfactors())
        left = self.F.left(i)
        for k in range(0, i+1):
            if isFaust(left):
                fac = left.factors(k)
            else:
                fac = left
            if not isinstance(fac, np.ndarray):
                fac = fac.toarray()
            Ffac = self.F.factors(k)
            if not isinstance(Ffac, np.ndarray):
                Ffac = Ffac.toarray()
            self.assertTrue(np.allclose(fac, Ffac))

    def test_save(self):
        print("Faust.save")
        tmp_dir = tempfile.gettempdir()+os.sep
        rand_suffix = random.Random().randint(1,1000)
        test_file = tmp_dir+"A"+str(rand_suffix)+".mat"
        self.F.save(test_file)
        Fs = Faust(test_file)
        self._assertAlmostEqual(Fs, self.F)

    def test_astype(self):
        print("test Faust.astype, Faust.dtype")
        if self.F.dtype == np.float:
            self.assertEqual(self.F.astype(np.complex).dtype, np.complex)
        else:
            self.assertEqual(self.F.astype(np.float).dtype, np.float)

    def test_pinv(self):
        print("Faust.pinv")
        self._assertAlmostEqual(self.F.pinv(), np.linalg.pinv(self.F.toarray()))

    def test_issparse(self):
        print("Faust.issparse")
        self.assertEqual(self.F.issparse(), np.all([isinstance(self.F.factors(i),
                                                        csr_matrix) for i in
                                             range(0, self.F.numfactors())]))

    def test_swap_cols(self):
        print("Faust.swap_cols")
        j1 = randint(0, self.F.shape[1])
        j2 = randint(0, self.F.shape[1])
        sF = self.F.swap_cols(j1, j2)
        Fa = self.F.toarray()
        sFa = sF.toarray()
        self._assertAlmostEqual(sFa[:,j1], Fa[:,j2])
        self.assertAlmostEqual(sF.norm(), self.F.norm())

    def test_swap_rows(self):
        print("Faust.swap_rows")
        i1 = randint(0, self.F.shape[0])
        i2 = randint(0, self.F.shape[0])
        sF = self.F.swap_rows(i1, i2)
        Fa = self.F.toarray()
        sFa = sF.toarray()
        self._assertAlmostEqual(sFa[i1,:], Fa[i2,:])
        self.assertAlmostEqual(sF.norm(), self.F.norm())

    def test_optimize_memory(self):
        print("Faust.optimize_memory")
        self.assertLessEqual(self.F.optimize_memory().nbytes, self.F.nbytes)


if __name__ == "__main__":
    if(len(sys.argv)> 1):
        dev = sys.argv[1]
        if dev != 'cpu' and not dev.startswith('gpu'):
            raise ValueError("dev argument must be cpu or gpu.")
        del sys.argv[1] # deleted to avoid interfering with unittest
    if(len(sys.argv) > 1):
        #ENOTE: test only a single test if name passed on command line
        singleton = unittest.TestSuite()
        singleton.addTest(PyfaustSimpleTest(sys.argv[1]))
        unittest.TextTestRunner().run(singleton)
    else:
        unittest.main()




