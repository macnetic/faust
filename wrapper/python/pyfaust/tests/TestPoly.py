import unittest
from pyfaust.poly import basis, poly, expm_multiply
from pyfaust import isFaust
import numpy as np
from numpy.linalg import norm
from scipy.sparse import csr_matrix, random, vstack as spvstack, eye as speye
from scipy.sparse.linalg import expm_multiply as scipy_expm_multiply


class TestPoly(unittest.TestCase):

    def __init__(self, methodName='runTest', dev='cpu', dtype='double'):
        super(TestPoly, self).__init__(methodName)
        self.dev = dev
        if dtype == 'real':  # backward compat
            dtype = 'double'
        self.dtype = dtype

    def setUp(self):
        self.d = 50
        self.density = 0.02
        self.L = random(self.d, self.d, .02, format='csr', dtype=self.dtype)
        self.L @= self.L.H
        self.K = 5

    def verif_basis(self, d, F, K, L):
        # assert the dimensions are consistent to L
        self.assertEqual(F.shape[0], (K+1)*L.shape[0])
        self.assertEqual(F.shape[1], L.shape[0])
        # assert the 0-degree polynomial matrix is the identity
        last_fac = F.factors(F.numfactors()-1).toarray()
        Id = np.eye(d)
        self.assertTrue(np.allclose(Id, last_fac))
        if K >= 1:
            # assert the 1-degree polynomial matrix is in the form [Id ; L]
            deg1_fac = F.factors(F.numfactors()-2).toarray()
            self.assertTrue(np.allclose(deg1_fac[:d, :], Id))
            self.assertTrue(np.allclose(deg1_fac[d:2*d, :], L.toarray()))
            if K >= 2:
                # assert the 2-degree polynomial matrix is
                # in the form [Id ; [-Id, L]]
                I2d = np.eye(2*d)
                deg2_fac = F.factors(F.numfactors()-3).toarray()
                self.assertTrue(np.allclose(deg2_fac[:2*d, :], I2d))
                self.assertTrue(np.allclose(deg2_fac[2*d:, :d], -Id))
                self.assertTrue(np.allclose(deg2_fac[2*d:, d:], 2*L.toarray()))
                if K >= 3:
                    # assert the n-degree polynomial matrix is in the form
                    # [I_nd ; [0 , -Id, 2L]]
                    for n in range(3, K):
                        Ind = np.eye(n*d)
                        degn_fac = F.factors(F.numfactors()-n-1).toarray()
                        self.assertTrue(np.allclose(degn_fac[:n*d, :], Ind))
                        self.assertTrue(np.allclose(degn_fac[n*d:, -2*d:-d],
                                                    -Id))
                        self.assertTrue(np.allclose(degn_fac[n*d:, -d:],
                                                    2*L.toarray()))
                        zero_part = degn_fac[n*d:, :-2*d]
                        self.assertTrue(np.linalg.norm(zero_part) == 0)

    def test_basis_py(self):
        print("Test basis(impl='py')")
        d = self.d
        L = self.L
        K = self.K
        F = basis(L, K, 'chebyshev', dev=self.dev, impl='py')
        self.verif_basis(d, F, K, L)

    def test_basis(self):
        print("Test basis(impl='native')")
        d = self.d
        L = self.L
        K = self.K
        F = basis(L, K, 'chebyshev', dev=self.dev, impl='native')
        self.verif_basis(d, F, K, L)

    def test_basisT0(self):
        print("Test basis(T0)")
        d = self.d
        L = self.L
        K = self.K
        density = self.density
        T0 = random(d, 2, density, format='csr', dtype=self.dtype)
        F = basis(L, K, 'chebyshev', dev=self.dev, T0=T0)
        print(F)
        # assert the dimensions are consistent to L and TO
        self.assertEqual(F.shape[0], (K+1)*L.shape[0])
        self.assertEqual(F.shape[1], T0.shape[1])
        # assert the 0-degree polynomial matrix is T0
        last_fac = F.factors(F.numfactors()-1).toarray()
        self.assertTrue(np.allclose(T0.toarray(), last_fac))

    def test_poly(self):
        print("Test poly()")
        self._test_poly_impl('native')

    def test_poly_py(self):
        print("Test poly(impl='py')")
        self._test_poly_impl('py')

    def _test_poly_impl(self, impl):
        from pyfaust import Faust
        d = self.d
        L = self.L
        K = self.K
        for L in [self.L, Faust(self.L)]:
            if impl == 'native' and isFaust(L):
                # native impl doesn't handle Faust L
                self.assertRaises(TypeError, basis, L, K, 'chebyshev',
                                  dev=self.dev, impl=impl)
                continue
            F = basis(L, K, 'chebyshev',
                      dev=self.dev, impl=impl).astype(self.dtype)
            self.assertEqual(F.shape[0], (K+1) * L.shape[0])
            coeffs = np.random.rand(K+1).astype(self.dtype)
            G = poly(coeffs, F, impl=impl)
            # Test polynomial as Faust
            poly_ref = np.zeros((d, d))
            for i, c in enumerate(coeffs[:]):
                poly_ref += c * F[d*i:(i+1)*d, :]
            self.assertAlmostEqual((G-poly_ref).norm(), 0)
            # Test polynomial as array
            GM = poly(coeffs, F.toarray(), impl=impl)
            self.assertTrue(isinstance(GM, np.ndarray))
            err = norm(GM - poly_ref.toarray())/norm(poly_ref.toarray())
            self.assertLessEqual(err, 1e-6)
            # Test polynomial-vector product
            x = np.random.rand(F.shape[1], 1).astype(L.dtype)
            # Three ways to do (not all as efficient as each other)
            Fx1 = poly(coeffs, F, dev=self.dev, impl=impl)@x
            Fx2 = poly(coeffs, F@x, dev=self.dev, impl=impl)
            Fx3 = poly(coeffs, F, X=x, dev=self.dev, impl=impl)
            err = norm(Fx1-Fx2)/norm(Fx1)
            self.assertLessEqual(err, 1e-6)
            self.assertTrue(np.allclose(Fx1, Fx3))
            # Test polynomial-matrix product
            X = np.random.rand(F.shape[1], 18).astype(L.dtype)
            FX1 = poly(coeffs, F, dev=self.dev, impl=impl)@X
            FX2 = poly(coeffs, F@X, dev=self.dev, impl=impl)
            FX3 = poly(coeffs, F, X=X, dev=self.dev, impl=impl)
            err = norm(FX1-FX2)/norm(FX1)
            self.assertLessEqual(err, 1e-6)
            self.assertTrue(np.allclose(FX2, FX3))
            # Test creating the polynomial basis on the fly
            G2 = poly(coeffs, 'chebyshev', L, impl=impl)
            self.assertAlmostEqual((G-G2).norm(), 0)
            GX = poly(coeffs, 'chebyshev', L, X=X, dev=self.dev, impl=impl)
            err = norm(FX1-GX)/norm(FX1)
            self.assertLessEqual(err, 1e-6)
            # Test polynomial-matrix product with arbitrary T0
            F_ = basis(L, K, 'chebyshev', dev=self.dev, T0=csr_matrix(X),
                       impl=impl)
            GT0eqX = poly(coeffs, F_, dev=self.dev, impl=impl).toarray()
            self.assertTrue(np.allclose(GT0eqX, FX1))

    def test_expm_multiply(self):
        print("Test expm_multiply()")
        L = self.L
        L = L@L.T
        # test expm_multiply on a vector
        x = np.random.rand(L.shape[1]).astype(L.dtype)
        pts_args = {'start': -.5, 'stop': -0.1, 'num': 3, 'endpoint': True}
        t = np.linspace(**pts_args).astype(L.dtype)
        y = expm_multiply(L, x, t)
        y_ref = scipy_expm_multiply(L, x, **pts_args)
        self.assertTrue(norm(y-y_ref)/norm(y_ref) < 1e-2)
        # test expm_multiply on a matrix
        X = np.random.rand(L.shape[1], 32).astype(L.dtype)
        pts_args = {'start': -.5, 'stop': -0.1, 'num': 3, 'endpoint': True}
        t = np.linspace(**pts_args)
        y = expm_multiply(L, X, t)
        y_ref = scipy_expm_multiply(L, X, **pts_args)
        self.assertTrue(norm(y-y_ref)/norm(y_ref) < 1e-2)
        # test expm_multiply with (non-default) tradeoff=='memory'
        X = np.random.rand(L.shape[1], 32).astype(L.dtype)
        pts_args = {'start': -.5, 'stop': -0.1, 'num': 3, 'endpoint': True}
        t = np.linspace(**pts_args)
        y = expm_multiply(L, X, t, tradeoff='memory')
        y_ref = scipy_expm_multiply(L, X, **pts_args)
        self.assertTrue(norm(y-y_ref)/norm(y_ref) < 1e-2)
        # test expm_multiply with (non-default) group_coeffs=False, poly_meth=2
        X = np.random.rand(L.shape[1], 32).astype(L.dtype)
        pts_args = {'start': -.5, 'stop': -0.1, 'num': 3, 'endpoint': True}
        t = np.linspace(**pts_args)
        y = expm_multiply(L, X, t, group_coeffs=False, poly_meth=2)
        y_ref = scipy_expm_multiply(L, X, **pts_args)
        self.assertTrue(norm(y-y_ref)/norm(y_ref) < 1e-2)
        # error cases
        self.assertRaisesRegex(TypeError, 'A must be a csr_matrix',
                               expm_multiply, L.toarray(), X, t)
        self.assertRaisesRegex(ValueError, 'A must be symmetric '
                               'positive definite.',
                               expm_multiply,
                               spvstack((L,
                                         speye(1, L.shape[1]))).tocsr(),
                               X, t)
        self.assertRaisesRegex(ValueError,
                               "tradeoff must be 'memory' or 'time'",
                               expm_multiply, L, X, t,
                               tradeoff='anything')
        self.assertRaisesRegex(ValueError, 'poly_meth must be 1 or 2',
                               expm_multiply, L, X, t, poly_meth=3)
        self.assertRaisesRegex(ValueError, 'group_coeffs must be a bool',
                               expm_multiply, L, X, t, group_coeffs='anything')
        self.assertRaisesRegex(ValueError,
                               "group_coeffs can't be True if poly_meth == 1.",
                               expm_multiply, L, X, t, poly_meth=1,
                               group_coeffs=True)
        t_non_neg_err = 'pyfaust.poly.expm_multiply handles only negative '
        'time points.'
        self.assertRaisesRegex(ValueError, t_non_neg_err,
                               expm_multiply, L, X, - t, poly_meth=2,
                               group_coeffs=True)
        self.assertRaisesRegex(ValueError, t_non_neg_err,
                               expm_multiply, L, X, - t)

    def test_poly_cat(self):
        print("Test poly._cat")
        from pyfaust.poly import _cat
        from pyfaust import rand as frand
        from scipy.sparse import random as srand
        type_err = "Inconsistent cat"
        axis_err = "axis must be 0 or 1"
        # only test errors because, functional cases
        # are tested indirectly in other tests
        # errors related to dimensions are not tested
        # but it depends on pyfaust/scipy core cat (assumed already tested)
        F1 = frand(10, 10, dtype=self.dtype)
        S1 = srand(10, 10, density=.5, dtype=self.dtype)
        for axis in [0, 1]:
            # 1. try to cat a csr_matrix/Faust
            self.assertRaisesRegex(TypeError, type_err, _cat, (F1, S1), axis)
            # 2. try to cat a csr_matrix/np.ndarray
            self.assertRaisesRegex(TypeError, type_err, _cat, (F1.toarray(),
                                                               S1), axis)
            # 3. try to cat a Faust/np.ndarray
            self.assertRaisesRegex(TypeError, type_err, _cat, (F1,
                                                               S1.toarray()),
                                   axis)
        # 4. test bad axis (2)
        self.assertRaisesRegex(ValueError, axis_err, _cat, (F1,
                                                            S1.toarray()),
                               2)

    def test_invm_multiply(self):
        print("Test invm_multiply()")
        from scipy.sparse import random
        try:
            from pyfaust.poly import invm_multiply
            from numpy.linalg import inv
            np.random.seed(42)  # for reproducibility
            A = random(64, 64, .1, format='csr')
            A = A@A.T
            B = np.random.rand(A.shape[1], 2)
            # tradeoff=='time' is already tested in docstrings
            A_inv_B = invm_multiply(A, B, rel_err=1e-2, max_K=2048,
                                    tradeoff='memory')
            A_inv_B_ref = inv(A.toarray())@B
            self.assertTrue(norm(A_inv_B-A_inv_B_ref)/norm(A_inv_B) <= 0.0273)
            # test error cases
            self.assertRaisesRegex(ValueError, 'poly_meth must be 1 or 2',
                                   invm_multiply, A, B, poly_meth=3)
            self.assertRaisesRegex(ValueError,
                                   "tradeoff must be 'memory' or 'time'",
                                   invm_multiply, A, B, tradeoff='anything')
            # test singular A error case
            sA = csr_matrix(A)
            sA[-2:, :] = 0
            self.assertRaisesRegex(Exception,
                                   "a <= 0 error: A is a singular matrix or its"
                                   " spectrum contains negative values.",
                                   invm_multiply, A, B)
        except ImportError:
            print("invm_multiply is only avaiable in experimental packages.")
