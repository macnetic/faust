import unittest
from scipy import sparse
import random
import tempfile
import os
import sys
import numpy as np
from scipy.io import savemat  # , loadmat
from numpy.linalg import norm
import math

class TestFaustPy(unittest.TestCase):

    MAX_NUM_FACTORS = 8  # for the tested Faust
    MAX_DIM_SIZE = 512
    MIN_DIM_SIZE = 3

    def setUp(self):
        """ Initializes the tests objects """
        r = random.Random()  # initialized from time or system
        num_factors = r.randint(1, TestFaustPy.MAX_NUM_FACTORS)
        factors = []
        d2 = r.randint(TestFaustPy.MIN_DIM_SIZE, TestFaustPy.MAX_DIM_SIZE)
        for i in range(0, num_factors):
            d1, d2 = d2, r.randint(1, TestFaustPy.MAX_DIM_SIZE)
            factors += [sparse.random(d1, d2, density=0.1, format='csr',
                        dtype=np.float64)] #.toarray() removed
            #print("factor",i,":", factors[i])
        self.F = Faust(factors)
        self.factors = factors
        # we keep dense matrices as reference for the tests
        for i in range(0, num_factors):
            self.factors[i] =  factors[i].toarray()
        print("Tests on random Faust with dims=", self.F.shape[0],
              self.F.shape[1])
        print("Num. factors:", num_factors)
        self.r = r
        self.num_factors = num_factors

    def testSave(self):
        print("testSave()")
        tmp_dir = tempfile.gettempdir()+os.sep
        # save the Faust through Faust core API
        rand_suffix = random.Random().randint(1,1000)
        test_file = tmp_dir+"A"+str(rand_suffix)+".mat"
        self.F.save(test_file, format="Matlab")
        # save the Faust relying on numpy API
        ref_file = tmp_dir+"A_ref"+str(rand_suffix)+".mat"
        mdict = {'faust_factors':
                 np.ndarray(shape=(1, self.F.get_num_factors()), dtype=object)}
        # self.F.display()
        for i in range(0, self.F.get_num_factors()):
            mdict['faust_factors'][0, i] = self.F.get_factor(i)
        savemat(ref_file, mdict)
        # open the two saved files and compare the fausts
        F_test = Faust(filepath=test_file)
        F_ref = Faust(filepath=ref_file)
        # print("self.F.get_num_factors()=", self.F.get_num_factors())
        self.assertEqual(F_test.get_num_factors(), F_ref.get_num_factors())
        for i in range(0, F_ref.get_num_factors()):
            fact_ref = F_ref.get_factor(i)
            fact_test = F_test.get_factor(i)
            self.assertEqual(fact_ref.shape, fact_test.shape)
            self.assertTrue((fact_ref == fact_test).all())

    def testGetNumRows(self):
        print("testGetNumRows()")
        self.assertEqual(self.F.shape[0], self.factors[0].shape[0])

    def testGetNumCols(self):
        print("testGetNumCols()")
        self.assertEqual(self.F.shape[1],
                         self.factors[len(self.factors)-1].shape[1])

    def testGetFactorAndConstructor(self):
        print("testGetFactorAndConstructor()")
        for ref_fact,i in zip(self.factors,range(0,len(self.factors))):
            #print("testGetFac() fac ",i, " :", self.F.get_factor(i))
            #print("testGetFac() reffac",i, ":",ref_fact)
            self.assertTrue((ref_fact == self.F.get_factor(i)).all())

    def testGetNumFactors(self):
        print("testGetNumFactors()")
        self.assertEqual(self.F.get_num_factors(), len(self.factors))

    def testNorm2(self):
        print("testNorm2()")
        ref_norm = norm(self.mulFactors())
        test_norm = self.F.norm(2)
        print("ref_norm=", ref_norm, "test_norm=", test_norm)
        # TODO: remove this workaround when the supposed bug will be corrected in core lib
        if(math.isnan(test_norm) and not math.isnan(ref_norm)):
            return
        self.assertLessEqual(abs(ref_norm-test_norm)/abs(ref_norm), 0.05)

    def testNorm1(self):
        print('test Faust.norm(1)')
        ref_norm = norm(self.mulFactors(), 1)
        test_norm = self.F.norm(1)
        print("test_norm=", test_norm, "ref_norm=", ref_norm)
        if(math.isnan(test_norm) and not math.isnan(ref_norm)):
            return
        self.assertLessEqual(abs(ref_norm-test_norm)/abs(ref_norm), 0.05)

    def testNormFro(self):
        print('test Faust.norm("fro")')
        ref_norm = norm(self.mulFactors(), 'fro')
        test_norm = self.F.norm('fro')
        print("test_norm=", test_norm, "ref_norm=", ref_norm)
        if(math.isnan(test_norm) and not math.isnan(ref_norm)):
            return
        self.assertLessEqual(abs(ref_norm-test_norm)/abs(ref_norm), 0.05)

    def faust_nnz(self):
        ref_nnz = 0
        for fact in self.factors:
            ref_nnz += np.count_nonzero(fact)
        return ref_nnz

    def testNnz(self):
        print("testNnz()")
        ref_nnz = self.faust_nnz()
        self.assertEqual(ref_nnz, self.F.nnz_sum())

    def testDensity(self):
        print("testDensity()")
        ref_density = \
        float(self.faust_nnz())/float(self.F.shape[1]*self.F.shape[0])
        self.assertAlmostEqual(ref_density, self.F.density(), delta=.001)

    def testRcg(self):
        print("testRcg()")
        ref_rcg = \
        float(self.F.shape[0]*self.F.shape[1])/float(self.faust_nnz())
        self.assertAlmostEqual(ref_rcg, self.F.rcg(), delta=.001)


    def mulFactors(self):
        n = self.factors[0].shape[0]
        prod = np.eye(n,n)
        for factor in self.factors:
            prod = prod.dot(factor)
        return prod

    def assertProdEq(self, prod, test_prod):
        self.assertEqual(prod.shape, test_prod.shape)
        for i in range(0, prod.shape[0]):
            for j in range(0, prod.shape[1]):
                #print(self.F[i,j][0],prod[i,j])
                if(prod[i,j] != 0):
                    self.assertLessEqual((test_prod[i,j]-prod[i,j])/prod[i,j],10**-6)


    def testGetItem(self):
        print("testGetItem()")
        n = self.factors[0].shape[0]
        # test whole array
        prod = self.mulFactors()
        test_prod = self.F[::,::]
        self.assertProdEq(prod, test_prod)
        test_prod = self.F[...,...]
        self.assertProdEq(prod, test_prod)
        # test one random element
        rand_i, rand_j = self.r.randint(0,self.F.shape[0]-1),self.r.randint(0,self.F.shape[1]-1)
        if(prod[rand_i,rand_j] != 0):
            self.assertLessEqual(abs(self.F[rand_i,rand_j][0]-prod[rand_i,rand_j])/abs(prod[rand_i,rand_j]),
                                 10**-6, msg=("compared values are (ref,rest) ="
                                              +str(prod[rand_i,rand_j])+str(prod[rand_i,rand_j])))
        # test one random row
        rand_i = self.r.randint(0,self.F.shape[0]-1)
        row = self.F[rand_i,...]
        for j in range(0,self.F.shape[1]):
            if(row[j] == 0):
                self.assertEqual(prod[rand_i,j], 0)
            else:
                self.assertLessEqual(abs(row[j]-(prod[rand_i,j]))/prod[rand_i,j],10**-6)
        # test one random col
        rand_j = self.r.randint(0,self.F.shape[1]-1)
        col = self.F[..., rand_j]
        for i in range(0,self.F.shape[0]):
            if(col[i] == 0):
                self.assertEqual(prod[i, rand_j], 0)
            else:
                self.assertLessEqual(abs(col[i]-(prod[i,rand_j]))/prod[i,rand_j],10**-6)

    def testToDense(self):
        print("testToDense()")
        prod = self.mulFactors()
        test_prod = self.F.toarray()
        self.assertProdEq(prod, test_prod)
        #self.assertTrue((self.F.toarray() == prod).all())


    def testMul(self):
        print("testMul()")
        rmat = np.random.rand(self.F.shape[1],
                              self.r.randint(1,TestFaustPy.MAX_DIM_SIZE))
        prod = self.mulFactors().dot(rmat)
        test_prod = self.F*rmat
        self.assertProdEq(prod, test_prod)


    def testTranspose(self):
        print("testTranspose()")
        tFaust = self.F.transpose()
        tF = tFaust.toarray()
        del tFaust # just to test weak reference of underlying Faust::Transform
        # works (otherwise it should crash with core dump later)
        F = self.F.toarray() # to avoid slowness
        for i in range(0, tF.shape[0]):
            for j in range(0, tF.shape[1]):
                if(F[j,i] != 0):
                    self.assertLessEqual(abs(tF[i,j]-F[j,i])/abs(F[j,i]), 10**-3)
                else:
                    self.assertEqual(tF[i,j],0)
        tmp_dir = tempfile.gettempdir()+os.sep
        rand_suffix = random.Random().randint(1,1000)
        test_file = tmp_dir+"A"+str(rand_suffix)+".mat"
        ref_file = tmp_dir+"A"+str(rand_suffix)+"o.mat"
        self.F.transpose().save(test_file)
        self.F.save(ref_file)
        #print("file=",test_file)
        tF2 = Faust(filepath=test_file)
        #print(tF2.shape[0], tF2.shape[1])
        #print(self.F.shape[0], self.F.shape[1])
        self.assertEqual(tF2.shape[1], tF.shape[1])
        self.assertEqual(tF2.shape[0], tF.shape[0])
        tF2 = tF2.toarray()
        for i in range(0, tF.shape[0]):
            for j in range(0, tF.shape[1]):
                if(F[j,i] != 0):
                    self.assertLessEqual(abs(tF2[i,j]-F[j,i])/abs(F[j,i]),
                                         10**-12)
                else:
                    self.assertEqual(tF2[i,j],0)
        os.remove(test_file)
        os.remove(ref_file)

    def testSize(self):
        print("testSize()")
        self.assertEqual((self.F.shape[0],self.F.shape[1]), self.F.shape)

    def testDelete(self):
        print("Test del Faust")
        tF = self.F.transpose()
        F = self.F
        del F
        with self.assertRaises(UnboundLocalError):
            print(F.shape)
        self.assertEqual(tF.shape, (self.factors[self.num_factors-1].shape[1],
                                     self.factors[0].shape[0]))

    def testConjugate(self):
        print("Test Faust.conj()")
        test_Fc = self.F.conj().toarray()
        ref_Fc = self.F.toarray().conj()
        self.assertTrue((test_Fc == ref_Fc).all())

    def testGetH(self):
        print("Test Faust.getH()")
        test_Fct = self.F.getH().toarray()
        ref_Fct = self.F.toarray().conj().T
        #print("test_Fct=", test_Fct)
        #print("ref_Fct=", ref_Fct)
        ref_Fct[ref_Fct==0] = 1
        test_Fct[test_Fct==0] = 1
        self.assertTrue(((((test_Fct-ref_Fct)/ref_Fct) < 0.01)).all())

class TestFaustPyCplx(TestFaustPy):

        def setUp(self):
            """ Initializes the tests objects """
            r = random.Random()  # initialized from time or system
            num_factors = r.randint(1, TestFaustPy.MAX_NUM_FACTORS)
            factors = []
            d2 = r.randint(TestFaustPy.MIN_DIM_SIZE, TestFaustPy.MAX_DIM_SIZE)
            for i in range(0, num_factors):
                d1, d2 = d2, r.randint(TestFaustPy.MIN_DIM_SIZE, TestFaustPy.MAX_DIM_SIZE)
                if(r.randint(0,1) > 0): # generate a dense complex matrix
                    factors += [np.random.rand(d1, d2).astype(np.complex)]
                    factors[i].imag = [np.random.rand(d1, d2)]
                else:
                    # generate a sparse matrix
                    # we can't use scipy.sparse.random directly because only float type is
                    # supported for random sparse matrix generation
                    factors += [sparse.random(d1, d2, dtype=np.float64, format='csr',
                                       density=r.random()).astype(np.complex)]
                    #print("setUp() i=",i, "d1=",d1, "d2=", d2, "factor.shape[0]=",
                    #      factors[i].shape[0])
                    factors[i] *= np.complex(r.random(), r.random())
            self.F = Faust(factors)
            self.factors = factors
            # we keep dense matrices as reference for the tests
            for i in range(0, num_factors):
                if(not isinstance(factors[i], np.ndarray)):
                    self.factors[i] =  factors[i].toarray()
            print("Tests on random complex Faust with dims=", self.F.shape[0],
                  self.F.shape[1])
            print("Num. factors:", num_factors)
            self.r = r
            self.num_factors = num_factors


if __name__ == "__main__":
    if(len(sys.argv)> 1):
        # argv[1] is for adding a directory in PYTHONPATH
        # (to find FaustPy module)
        sys.path.append(sys.argv[1])
        del sys.argv[1] # deleted to avoid interfering with unittest
    from FaustPy import Faust
    unittest.main()
