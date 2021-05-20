import unittest
from pyfaust.poly import basis
from numpy.random import randint
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import aslinearoperator
import tempfile
import os
import random

dev = 'cpu'
field = 'real'


class TestPoly(unittest.TestCase):

    def __init__(self, methodName='runTest', dev='cpu', field='real'):
        super(TestPoly, self).__init__(methodName)
        self.dev = dev
        self.field = field

    def setUp(self):
        pass

    def test_basis(self):
        print("Test basis()")
        from scipy.sparse import random
        d = 50
        if self.field == 'complex':
            dtype = 'complex'
        else:
            dtype = 'double'
        L = random(d, d, .02, format='csr', dtype=dtype)
        L = L@L.T
        K = 5
        F = basis(L, K, 'chebyshev', dev=self.dev)
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
                # assert the 2-degree polynomial matrix is in the form [Id ; [-Id, L]]
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
                        self.assertTrue(np.allclose(degn_fac[n*d:, -2*d:-d], -Id))
                        self.assertTrue(np.allclose(degn_fac[n*d:, -d:],
                                                    2*L.toarray()))
                        zero_part = degn_fac[n*d:, :-2*d]
                        self.assertTrue(np.linalg.norm(zero_part) == 0)




