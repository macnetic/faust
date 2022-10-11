import unittest
import pyfaust as pf
from pyfaust.lazylinop import LazyLinearOp, vstack, hstack
import numpy.linalg as LA
import numpy as np

class TestLazyLinearOpFaust(unittest.TestCase):

    def setUp(self):
        self.lop = LazyLinearOp.create(pf.rand(10, 15))
        self.lopA = self.lop.toarray()
        self.lop2 = LazyLinearOp.create(pf.rand(10, 15))
        self.lop2A = self.lop2.toarray()
        self.lop3 = LazyLinearOp.create(pf.rand(15, 10))
        self.lop3A = self.lop3.toarray()

    def test_shape(self):
        self.assertEqual(self.lop.shape, self.lop.eval().shape)

    def test_ndim(self):
        self.assertEqual(self.lop.ndim, 2)

    def test_transp(self):
        lopT = self.lop.T
        self.assertAlmostEqual(LA.norm(lopT.toarray()-self.lopA.T), 0)
        lopT = self.lop.transpose()
        self.assertAlmostEqual(LA.norm(lopT.toarray()-self.lopA.T), 0)
        self.assertEqual(self.lop.shape[0], lopT.shape[1])
        self.assertEqual(self.lop.shape[1], lopT.shape[0])

    def test_conj(self):
        lopC = self.lop.conj()
        self.assertAlmostEqual(LA.norm(lopC.toarray()-self.lopA.conj()), 0)

    def test_adjoint(self):
        lopH = self.lop.H
        self.assertAlmostEqual(LA.norm(lopH.toarray()-self.lopA.conj().T), 0)
        lopH = self.lop.getH()
        self.assertAlmostEqual(LA.norm(lopH.toarray()-self.lopA.conj().T), 0)
        self.assertEqual(self.lop.shape[0], lopH.shape[1])
        self.assertEqual(self.lop.shape[1], lopH.shape[0])

    def test_add(self):
        ladd = self.lop + self.lop2
        self.assertAlmostEqual(LA.norm(ladd.toarray()-(self.lopA+self.lop2A)),
                               0)
        M = np.random.rand(*self.lop.shape)
        ladd2 = self.lop + M
        self.assertTrue(isinstance(ladd2, LazyLinearOp))
        self.assertAlmostEqual(LA.norm(ladd2.toarray()-(self.lopA + M)),
                               0)
    def test_iadd(self):
        self.assertRaises(NotImplementedError, self.lop.__iadd__, self.lop2)

    def test_radd(self):
        M = np.random.rand(*self.lop.shape)
        ladd2 = M + self.lop
        self.assertTrue(isinstance(ladd2, LazyLinearOp))
        self.assertAlmostEqual(LA.norm(ladd2.toarray()-(M + self.lopA)), 0)

    def test_sub(self):
        lsub = self.lop - self.lop2
        self.assertAlmostEqual(LA.norm(lsub.toarray() - (self.lopA - self.lop2A)),
                               0)
        M = np.random.rand(*self.lop.shape)
        lsub2 = self.lop - M
        self.assertTrue(isinstance(lsub2, LazyLinearOp))
        self.assertAlmostEqual(LA.norm(lsub2.toarray() - (self.lopA - M)), 0)

    def test_rsub(self):
        M = np.random.rand(*self.lop.shape)
        lsub2 = M - self.lop
        self.assertTrue(isinstance(lsub2, LazyLinearOp))
        self.assertAlmostEqual(LA.norm(lsub2.toarray()-(M - self.lopA)), 0)

    def test_isub(self):
        self.assertRaises(NotImplementedError, self.lop.__isub__, self.lop2)

    def test_matmul_dot_matvec(self):
        lmul = self.lop @ self.lop3
        self.assertAlmostEqual(LA.norm(lmul.toarray() - (self.lopA @ self.lop3A)),
                               0)
        lmul = self.lop.dot(self.lop3)
        self.assertAlmostEqual(LA.norm(lmul.toarray() - (self.lopA @ self.lop3A)),
                               0)
        M = np.random.rand(self.lop.shape[1], 15)
        lmul2 = self.lop @ M
        self.assertTrue(isinstance(lmul2, np.ndarray))
        self.assertAlmostEqual(LA.norm(lmul2 - (self.lopA @ M)),
                               0)
        lmul2 = self.lop.matvec(M[:,0])
        self.assertTrue(isinstance(lmul2, np.ndarray))
        self.assertAlmostEqual(LA.norm(lmul2 - (self.lopA @ M[:,0])),
                               0)

    def test_rmatmul(self):
        M = np.random.rand(15, self.lop.shape[0])
        lmul2 = M @ self.lop
        self.assertTrue(isinstance(lmul2, np.ndarray))
        self.assertAlmostEqual(LA.norm(lmul2 - (M @ self.lopA)),
                               0)

    def test_imatmul(self):
        self.assertRaises(NotImplementedError, self.lop.__imatmul__, self.lop2)

    def test_mul(self):
        v = np.random.rand(self.lop.shape[1])
        lmul2 = self.lop * v
        self.assertTrue(isinstance(lmul2, LazyLinearOp))
        self.assertAlmostEqual(LA.norm(lmul2.toarray() - (self.lopA * v)),
                               0)
        v = np.random.rand(1, self.lop.shape[1])
        lmul2 = self.lop * v
        self.assertTrue(isinstance(lmul2, LazyLinearOp))
        self.assertAlmostEqual(LA.norm(lmul2.toarray() - (self.lopA * v)),
                               0)
        s = np.random.rand(1, 1)[0, 0]
        lmul2 = self.lop * s
        self.assertTrue(isinstance(lmul2, LazyLinearOp))
        self.assertAlmostEqual(LA.norm(lmul2.toarray() - (self.lopA * s)),
                               0)

    def test_rmul(self):
        v = np.random.rand(self.lop.shape[1])
        lmul2 = v * self.lop
        self.assertTrue(isinstance(lmul2, LazyLinearOp))
        self.assertAlmostEqual(LA.norm(lmul2.toarray() - (v * self.lopA)),
                               0)

        v = np.random.rand(1, self.lop.shape[1])
        lmul2 = v * self.lop
        self.assertTrue(isinstance(lmul2, LazyLinearOp))
        self.assertAlmostEqual(LA.norm(lmul2.toarray() - (v * self.lopA)),
                               0)

        s = np.random.rand(1, 1)[0, 0]
        lmul2 = s * self.lop
        self.assertTrue(isinstance(lmul2, LazyLinearOp))
        self.assertAlmostEqual(LA.norm(lmul2.toarray() - (s * self.lopA)),
                               0)

    def test_concatenate(self):
        lcat = self.lop.concatenate(self.lop2, axis=0)
        self.assertAlmostEqual(LA.norm(lcat.toarray() - np.vstack((self.lopA,
                                                                   self.lop2A))),
                               0)
        self.assertEqual(lcat.shape[0], self.lop.shape[0] + self.lop2.shape[0])
        lcat = self.lop.concatenate(self.lop2, axis=1)
        self.assertAlmostEqual(LA.norm(lcat.toarray() - np.hstack((self.lopA,
                                                                   self.lop2A))),
                               0)
        self.assertEqual(lcat.shape[1], self.lop.shape[1] + self.lop2.shape[1])
        # auto concat
        lcat = self.lop.concatenate(self.lop, axis=0)
        self.assertAlmostEqual(LA.norm(lcat.toarray() - np.vstack((self.lopA,
                                                                   self.lopA))),
                               0)
        self.assertEqual(lcat.shape[0], self.lop.shape[0] + self.lop.shape[0])
        # using hstack and vstack
        lcat = vstack((self.lop, self.lop2, self.lop))
        self.assertAlmostEqual(LA.norm(lcat.toarray() - np.vstack((self.lopA,
                                                                   self.lop2A,
                                                                  self.lopA))),
                               0)
        self.assertEqual(lcat.shape[0], 2 * self.lop.shape[0] + self.lop2.shape[0])
        lcat = hstack((self.lop, self.lop2, self.lopA))
        self.assertAlmostEqual(LA.norm(lcat.toarray() - np.hstack((self.lopA,
                                                                   self.lop2A,
                                                                  self.lopA))),
                               0)
        self.assertEqual(lcat.shape[1], 2 * self.lop.shape[1] + self.lop2.shape[1])



    def test_chain_ops(self):
        lchain = self.lop + self.lop2
        lchain = lchain @ self.lop3
        lchain = 2 * lchain
        v = np.random.rand(lchain.shape[1])
        lchain = lchain * v
        lchain = lchain.concatenate(self.lop3, axis=0)
        mat_ref = np.vstack(((2 * (self.lopA + self.lop2A) @ self.lop3A) * v,
                             self.lop3A))
        self.assertAlmostEqual(LA.norm(lchain.toarray() - mat_ref),
                               0)

    def test_get_item(self):
        n1 = self.lop.shape[0]//2
        n2 = self.lop.shape[1]//2
        lslice = self.lop[3:n1, 3:n2]
        lsliceA = self.lopA[3:n1, 3:n2]
        self.assertAlmostEqual(LA.norm(lslice.toarray()-lsliceA), 0)

    def test_real(self):
        cF = pf.rand(10, 15, field='complex')
        lcF = LazyLinearOp.create(cF)
        lcF = lcF.real
        self.assertAlmostEqual(LA.norm(lcF.toarray()-cF.real.toarray()), 0)

    def test_imag(self):
        cF = pf.rand(10, 15, field='complex')
        lcF = LazyLinearOp.create(cF)
        lcF = lcF.imag
        self.assertAlmostEqual(LA.norm(lcF.toarray()-cF.imag.toarray()), 0)

    def test_aslazylinop(self):
        from pyfaust.lazylinop import asLazyLinearOp
        cF = pf.rand(10, 15, field='complex')
        lcF = asLazyLinearOp(cF)
        self.assertTrue(pf.lazylinop.isLazyLinearOp(lcF))
        self.assertEqual(cF.shape, lcF.shape)

class TestLazyLinearOpFaustKron(TestLazyLinearOpFaust):

    def setUp(self):
        from pyfaust.lazylinop import kron as lkron
        lop_A = LazyLinearOp.create(pf.rand(10, 15))
        lop_B = LazyLinearOp.create(pf.rand(10, 15))
        self.lop = lkron(lop_A, lop_B)
        self.lopA = self.lop.toarray()
        self.lop2 = LazyLinearOp.create(pf.rand(*self.lop.shape))
        self.lop2A = self.lop2.toarray()
        self.lop3 = LazyLinearOp.create(pf.rand(self.lop.shape[1], 10))
        self.lop3A = self.lop3.toarray()

if '__main__' == __name__:
    unittest.main()
