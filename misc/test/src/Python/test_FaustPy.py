import unittest
from scipy import sparse
import random
import tempfile
import os
import sys
import numpy as np
from scipy.io import savemat,loadmat
from numpy.linalg import norm
from scipy.sparse import spdiags, bsr_matrix
import math
from os.path import dirname

class TestFaustPy(unittest.TestCase):

    MAX_NUM_FACTORS = 4  # for the tested Faust
    MAX_DIM_SIZE = 256
    MIN_DIM_SIZE = 3

    def setUp(self):
        """ Initializes the tests objects """
        r = random.Random()  # initialized from time or system
        num_factors = r.randint(1, TestFaustPy.MAX_NUM_FACTORS)
        factors = []
        d2 = r.randint(TestFaustPy.MIN_DIM_SIZE, TestFaustPy.MAX_DIM_SIZE)
        for i in range(0, num_factors):
            d1, d2 = d2, r.randint(1, TestFaustPy.MAX_DIM_SIZE)
            factors += [sparse.random(d1, d2, density=0.5, format='csr',
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
                 np.ndarray(shape=(1, self.F.numfactors()), dtype=object)}
        # self.F.display()
        for i in range(0, self.F.numfactors()):
            mdict['faust_factors'][0, i] = self.F.factors(i)
        savemat(ref_file, mdict)
        # open the two saved files and compare the fausts
        F_test = Faust(filepath=test_file)
        F_ref = Faust(filepath=ref_file)
        # print("self.F.numfactors()=", self.F.numfactors())
        self.assertEqual(F_test.numfactors(), F_ref.numfactors())
        for i in range(0, F_ref.numfactors()):
            fact_ref = F_ref.factors(i)
            fact_test = F_test.factors(i)
            if(not isinstance(fact_ref, np.ndarray)):
                # fact_ref is a sparse matrix
                fact_ref = fact_ref.toarray()
            if(not isinstance(fact_test, np.ndarray)):
                # fact_test is a sparse matrix
                fact_test = fact_test.toarray()
            self.assertEqual(fact_test.shape, fact_ref.shape)
            self.assertTrue((fact_ref == fact_test).all())

    def testLoadNative(self):
        print("test load native")
        import pyfaust as pf
        F = pf.rand(32, 295, density=.5)
        F.save('rand_faust.mat')
        F_ = pf.Faust.load_native('rand_faust.mat')
        self.assertLessEqual((F-F_).norm()/F.norm(), 1e-6)

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
            #print("testGetFac() fac ",i, " :", self.F.factors(i))
            #print("testGetFac() reffac",i, ":",ref_fact)
            test_fact = self.F.factors(i)
            if(not isinstance(test_fact, np.ndarray)):
                # test_fact is a sparse matrix
                test_fact = test_fact.toarray()
            #print(ref_fact.shape, test_fact.shape)
            test_fact = np.asfortranarray(test_fact)
            self.assertTrue((ref_fact == test_fact).all())
        # and factors() on transpose Faust ?
        tF = self.F.transpose()
        for i in range(0,len(self.factors)):
            #print("testGetFac() fac ",i, " :", self.F.factors(i))
            #print("testGetFac() reffac",i, ":",ref_fact)
            test_fact = tF.factors(i)
            ref_fact = self.F.factors(len(self.factors)-i-1).transpose()
            if(not isinstance(test_fact, np.ndarray)):
                # test_fact is a sparse matrix
                test_fact = test_fact.toarray()
            #print(ref_fact.shape, test_fact.shape)
            test_fact = np.asfortranarray(test_fact)
            self.assertTrue((ref_fact == test_fact).all())

    def testGetNumFactors(self):
        print("testGetNumFactors()")
        self.assertEqual(self.F.numfactors(), len(self.factors))

    def testNormalize(self):
        print("test Faust.normalize()")
        test_args = [
            [],
            [2],
            ['fro'],
            [float('Inf')],
            [np.inf],
            [1],
            [2,0],
            ['fro', 0],
            [np.inf, 0],
            [1, 1],
            [2,1],
            ['fro', 1],
            [np.inf, 1],
            {'axis':0,'ord':1},
            {'axis':0, 'ord':2},
            {'axis':0, 'ord':np.inf},
            {'axis':1,'ord':1},
            {'axis':1, 'ord':2},
            {'axis':1, 'ord':np.inf}
        ]
        F = self.F
        for args in test_args:
            axis = 1 #default
            ord = 'fro' # default
            print('signature: ', args)
            if(isinstance(args,dict)):
                test_NF = F.normalize(**args)
                if('axis' in args.keys()):
                    axis = args['axis']
                if('ord' in args.keys()):
                    ord = args['ord']
            else:
                test_NF = F.normalize(*args)
                if(len(args) > 0):
                   ord = args[0]
                   if(len(args) > 1):
                        axis = args[1]
            ref_full_NF = F.toarray()
            # print("axis=", axis, "ord=", ord)
            i = self.r.randint(0, F.shape[axis]-1)
            #for i in range(0,F.shape[axis]):
            # test only one random column/row (speeder)
            if(axis == 0 and norm(ref_full_NF[i:i+1,:], ord) != 0):
                vec_ref_full_NF = ref_full_NF[i,:] = ref_full_NF[i,:]/norm(ref_full_NF[i:i+1,:], ord)
            elif(axis == 1 and norm(ref_full_NF[:,i:i+1], ord) != 0):
                vec_ref_full_NF = ref_full_NF[:,i] = \
                        ref_full_NF[:,i]/norm(ref_full_NF[:,i:i+1], ord)
            else:
                continue
            if(F.dtype == np.complex):
                places=1 # accuracy is less good with complex
            else:
                places=3
            if(axis==0):
                n1 = norm(ref_full_NF[i,:])
                n2 = test_NF[i,:].norm()
                vec_test_NF = test_NF[i,:]
            else: #axis == 1
                n1 = norm(ref_full_NF[:,i])
                n2 = test_NF[:,i].norm()
                vec_test_NF = test_NF[:,i]
            self.assertAlmostEqual(n1,n2,
                places=places,
                       msg="\nref_full_F=\n"+str(F.toarray())+"\nvec_ref_full_NF=\n"+str(vec_ref_full_NF)+"\nvec_test_NF=\n"+ \
                       str(vec_test_NF.toarray())+"\nF=\n"+ \
                       str(str(F.save('/tmp/normalize_test.mat'))+", i="+str(i)))


    def testNormInf(self):
        print("testNormInf()")
        ref_norm = norm(self.mulFactors(), np.inf)
        test_norm = self.F.norm(np.inf)
        self.F.save("norminf_test.mat");
        self.F.save('/tmp/F_inf_test.mat')
        print("ref_norm=", ref_norm, "test_norm=", test_norm,
              "test_toarray_norm=", norm(self.F.toarray(),np.inf))
        # TODO: remove this workaround when the supposed bug will be corrected in core lib
        if(math.isnan(test_norm) and not math.isnan(ref_norm)):
            return
        self.assertLessEqual(abs(ref_norm-test_norm)/abs(ref_norm), 0.05)

    def testNorm2(self):
        print("testNorm2()")
        ref_norm = norm(self.mulFactors(),2)
        test_norm = self.F.norm(2)
        print("ref_norm=", ref_norm, "test_norm=", test_norm)
        # TODO: remove this workaround when the supposed bug will be corrected in core lib
        if(math.isnan(test_norm) and not math.isnan(ref_norm)):
            return
        if(ref_norm == 0.0):
            self.assertLessEqual(test_norm, 0.0005)
        self.assertLessEqual(abs(ref_norm-test_norm)/ref_norm, 0.05)

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
        # mul. factors in same order than core C++ to in
        # Faust::Transform::product (calling multiply for that)
        n = self.factors[-1].shape[1]
        prod = np.eye(n,n)
        for i in range(len(self.factors)-1,-1,-1):
            prod = self.factors[i].dot(prod)
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
        test_prod = self.F[::,::].toarray()
        self.assertProdEq(prod, test_prod)
        # reminder the Faust is of minimal size 3 for the both dims (MIN_DIM_SIZE)
        # test_prod = self.F[...,...].toarray() # forbidden (only one index can
        # be an ellipsis)
        # self.assertProdEq(prod, test_prod)
        # test one random element
        rand_i, rand_j = self.r.randint(0,self.F.shape[0]-2),self.r.randint(0,self.F.shape[1]-2)
        if(prod[rand_i,rand_j] != 0):
            self.assertLessEqual(abs(self.F[rand_i,rand_j]-prod[rand_i,rand_j])/abs(prod[rand_i,rand_j]),
                                 10**-6, msg=("compared values are (ref,rest) ="
                                              +str(prod[rand_i,rand_j])+str(prod[rand_i,rand_j])))
        # test one random row
        rand_i = self.r.randint(0,self.F.shape[0]-2)
        row = self.F[rand_i,...].toarray()
        for j in range(0,self.F.shape[1]):
            if(row[0,j] == 0):
                self.assertEqual(prod[rand_i,j], 0)
            else:
                self.assertLessEqual(abs(row[0,j]-(prod[rand_i,j]))/prod[rand_i,j],10**-6)
        # test one random col
        rand_j = self.r.randint(0,self.F.shape[1]-2)
        col = self.F[..., rand_j].toarray()
        for i in range(0,self.F.shape[0]):
            if(col[i,0] == 0):
                self.assertEqual(prod[i, rand_j], 0)
            else:
                self.assertLessEqual(abs(col[i,0]-(prod[i,rand_j]))/prod[i,rand_j],10**-6)
        # test that the sliced Faust is consistent with the sliced full matrix
        # of the Faust
        rand_ii, rand_jj = \
        self.r.randint(rand_i+1,self.F.shape[0]-1),self.r.randint(rand_j+1,self.F.shape[1]-1)
        sF = self.F[rand_i:rand_ii,rand_j:rand_jj]
        n1 = norm(sF.toarray()-self.F.toarray()[rand_i:rand_ii,rand_j:rand_jj])
        n2 = norm(self.F.toarray()[rand_i:rand_ii,rand_j:rand_jj])
        if(n2 == 0): # avoid nan error
            self.assertLessEqual(n1,0.0005)
        else:
            print("slices:", rand_i, rand_ii, rand_j, rand_jj)
            self.F.save('test_sF.mat')
            self.assertLessEqual(n1/n2,0.005)
        # test a second slice on sF only if the shape allows it # min shape (1,1)
        if(sF.shape[0] >= 2 and sF.shape[1] >= 2):
            rand_i, rand_j = \
            self.r.randint(0,sF.shape[0]-2),self.r.randint(0,sF.shape[1]-2)
            rand_ii, rand_jj = \
                    self.r.randint(rand_i+1,sF.shape[0]-1),self.r.randint(rand_j+1,sF.shape[1]-1)
            sF2 = sF[rand_i:rand_ii,rand_j:rand_jj]
            n1 = norm(sF2.toarray()-sF.toarray()[rand_i:rand_ii,rand_j:rand_jj])
            n2 = norm(sF.toarray()[rand_i:rand_ii,rand_j:rand_jj])
            if(n2 == 0): # avoid nan error
                self.assertLessEqual(n1,0.0005)
            else:
                self.assertLessEqual(n1/n2,0.005)
        rand_i, rand_j = \
        self.r.randint(0,self.F.shape[0]-2),self.r.randint(0,self.F.shape[1]-2)
        rand_ii, rand_jj = \
        self.r.randint(rand_i+1,self.F.shape[0]-1),self.r.randint(rand_j+1,self.F.shape[1]-1)
        print("test fancy indexing on rows")
        F = self.F
        num_inds = self.r.randint(1,F.shape[0])
        row_ids = [ self.r.randint(0,F.shape[0]-1) for i in
                   range(0,num_inds)]
        F_rows = F[row_ids,:]
        self.assertTrue((F_rows.toarray() == F.toarray()[row_ids,:]).all())
        print("test fancy indexing on cols")
        num_inds = self.r.randint(1,F.shape[1])
        col_ids = [ int(self.r.randint(0,F.shape[1]-1)) for i in
                   range(0,num_inds)]
        F_cols = F[:,col_ids]
        n1 = norm(F_cols.toarray() - F.toarray()[:,col_ids])
        n2 = norm(F.toarray()[:,col_ids])
        if(n2 == 0): # avoid nan error
            self.assertLessEqual(n1,0.0005)
        else:
            self.assertLessEqual(n1/n2, 0.005)
        print("test fancy indexing on rows with slice on cols")
        num_inds = self.r.randint(1,F.shape[0]-1)
        col_slice_start = self.r.randint(0,F.shape[1]-1)
        col_slice_stop = self.r.randint(col_slice_start+1,F.shape[1])
        row_ids = [ self.r.randint(0,F.shape[0]-1) for i in
                   range(0,num_inds)]
        F_rows = F[row_ids,col_slice_start:col_slice_stop]
        n1 = \
        norm(F_rows.toarray()-F.toarray()[row_ids,col_slice_start:col_slice_stop])
        n2 = norm(F.toarray()[row_ids,col_slice_start:col_slice_stop])
        if(n2 == 0): # avoid nan error
            self.assertLessEqual(n1,0.0005)
        else:
            self.assertLessEqual(n1/n2, 0.005)
        print("test fancy indexing on cols with slice on rows")
        num_inds = self.r.randint(1,F.shape[1])
        col_ids = [ self.r.randint(0,F.shape[1]-1) for i in
                   range(0,num_inds)]
        row_slice_start = self.r.randint(0,F.shape[0]-1)
        row_slice_stop = self.r.randint(row_slice_start+1, F.shape[0])
        F_cols = F[row_slice_start:row_slice_stop,col_ids]
        n1 = norm(F_cols.toarray() -
                  F.toarray()[row_slice_start:row_slice_stop,col_ids])
        n2 = norm(F.toarray()[row_slice_start:row_slice_stop,col_ids])
        if(n2 == 0): # avoid nan error
            self.assertLessEqual(n1,0.0005)
        else:
            self.assertLessEqual(n1/n2, 0.005)
        print("test slicing with positive step equal or greater than one")
        # rows
        rstep = self.r.randint(1, F.shape[0])
        assert(rstep >= 1)
        sF = F[rand_i:rand_ii:rstep]
        n1 = norm(sF.toarray()-F.toarray()[rand_i:rand_ii:rstep])
        n2 = norm(F.toarray()[rand_i:rand_ii:rstep])
        if(n2 == 0): # avoid nan error
            self.assertLessEqual(n1,0.0005)
        else:
            self.assertLessEqual(n1/n2, 0.005)
        # columns
        cstep = self.r.randint(1, F.shape[1])
        assert(cstep >= 1)
        sF = F[:,rand_j:rand_jj:cstep]
        n1 = norm(sF.toarray()-F.toarray()[:,rand_j:rand_jj:cstep])
        n2 = norm(F.toarray()[:,rand_j:rand_jj:cstep])
        if(n2 == 0): # avoid nan error
            self.assertLessEqual(n1,0.0005)
        else:
            self.assertLessEqual(n1/n2, 0.005)
        # rows and cols
        sF = F[rand_i:rand_ii:rstep,rand_j:rand_jj:cstep]
        n1 = \
        norm(sF.toarray()-F.toarray()[rand_i:rand_ii:rstep,rand_j:rand_jj:cstep])
        n2 = norm(F.toarray()[rand_i:rand_ii:rstep,rand_j:rand_jj:cstep])
        if(n2 == 0): # avoid nan error
            self.assertLessEqual(n1,0.0005)
        else:
            self.assertLessEqual(n1/n2, 0.005)
        print("test slicing with negative step")
        # rows
        sF = F[rand_ii:rand_i:-rstep]
        n1 = norm(sF.toarray()-F.toarray()[rand_ii:rand_i:-rstep])
        n2 = norm(F.toarray()[rand_ii:rand_i:-rstep])
        if(n2 == 0): # avoid nan error
            self.assertLessEqual(n1,0.0005)
        else:
            self.assertLessEqual(n1/n2, 0.005)
        # cols
        sF = F[:,rand_jj:rand_j:-cstep]
        n1 = norm(sF.toarray()-F.toarray()[:,rand_jj:rand_j:-cstep])
        n2 = norm(F.toarray()[:,rand_jj:rand_j:-cstep])
        if(n2 == 0): # avoid nan error
            self.assertLessEqual(n1,0.0005)
        else:
            self.assertLessEqual(n1/n2, 0.005)
        # rows and cols
        sF = F[rand_ii:rand_i:-rstep,rand_jj:rand_j:-cstep]
        n1 = \
        norm(sF.toarray()-F.toarray()[rand_ii:rand_i:-rstep,rand_jj:rand_j:-cstep])
        n2 = norm(F.toarray()[rand_ii:rand_i:-rstep,rand_jj:rand_j:-cstep])
        if(n2 == 0): # avoid nan error
            self.assertLessEqual(n1,0.0005)
        else:
            self.assertLessEqual(n1/n2, 0.005)
        print("test slicing with start and stop indices being negative (relative "
              "position to the end: -1 for the last elt etc.")
        # assuming F.shape[0] > rand_ii > rand_i >= 0, likewise for rand_j,
        # rand_jj
        # reminder: slice(None) means 0:-1:1
        # test multiple confs in a loop
        for slice_row, slice_col in \
        [(slice(rand_i-F.shape[0],rand_ii-F.shape[0]),slice(None)), # conf1
         (slice(None),slice(rand_j-F.shape[1],rand_jj-F.shape[1])), # conf2
         # conf 3
         (slice(rand_i-F.shape[0],rand_ii-F.shape[0]),
          slice(rand_j-F.shape[1],rand_jj-F.shape[1])),
         #conf 4
         (slice(rand_i-F.shape[0],rand_ii-F.shape[0], rstep),slice(None)),
         #conf 5
         (slice(None),slice(rand_j-F.shape[1],rand_jj-F.shape[1],cstep))
        ]:
            sF = F[slice_row, slice_col]
            n1 = norm(sF.toarray()-F.toarray()[slice_row,slice_col])
            n2 = norm(F.toarray()[slice_row,slice_col])
            if(n2 == 0): # avoid nan error
                self.assertLessEqual(n1,0.0005)
            else:
                self.assertLessEqual(n1/n2, 0.005)

    def testLazySlicingIndexing(self):
        print("testLazySlicingIndexing()")
        from pyfaust import rand as frand
        for F in [frand(512, 512, fac_type='dense'),
                  frand(512, 512, fac_type='sparse')] + \
                 [frand(512, 512, fac_type='dense', num_factors=1),
                  frand(512, 512, fac_type='sparse', num_factors=1)]:
            # test that two slicings in a row make a consistent Faust
            self.assertTrue(np.allclose(F[:15, :42][:12, :13].toarray(),
                                        F.toarray()[:15, :42][:12, :13]))
            # the same but with a transpose between
            self.assertTrue(np.allclose(F[:15, :42].T[:12, :13].toarray(),
                                        F.toarray()[:15, :42].T[:12, :13]))
            # test that two indexings in a row make a consistent Faust
            I = list(range(0, 15, 3))
            J = list(range(0, 42, 5))
            I2 = list(range(0, len(I), 2))
            J2 = list(range(0, len(J), 2))
            self.assertTrue(np.allclose(F[:,J].toarray(),
                                        F.toarray()[:,J]))
            self.assertTrue(np.allclose(F[I][:,J][I2][:,J2].toarray(),
                                        F.toarray()[I][:,J][I2][:,J2]))
            # the same with a transpose between
            self.assertTrue(np.allclose(F[I][:,J].T[J2][:,I2].toarray(),
                                        F.toarray()[I][:,J].T[J2][:,I2]))
            # test that a slicing followed by an indexing works
            self.assertTrue(np.allclose(F[:15, :42][I][:,J].toarray(),
                                        F.toarray()[:15, :42][I][:,J]))
            # the same with a transpose between
            self.assertTrue(np.allclose(F[:15, :42].T[J][:,I].toarray(),
                                        F.toarray()[:15, :42].T[J][:,I]))
            # now do the same but the indexing first
            self.assertTrue(np.allclose(F[I][:,J][:len(I)//2, :len(J)//2].toarray(),
                                        F.toarray()[I][:,J][:len(I)//2, :len(J)//2]))
            self.assertTrue(np.allclose(F[I][:,J].T[:len(I)//2, :len(J)//2].toarray(),
                                        F.toarray()[I][:,J].T[:len(I)//2, :len(J)//2]))
            # other cases
            self.assertTrue(np.allclose(F[I,15:42].toarray(),
                                        F.toarray()[I,15:42]))
            self.assertTrue(np.allclose(F[15:42, J].toarray(),
                                        F.toarray()[15:42, J]))
            self.assertTrue(np.allclose(F.T[J,15:42].toarray(),
                                        F.T.toarray()[J,15:42]))
            self.assertTrue(np.allclose(F.T[15:42, I].toarray(),
                                        F.T.toarray()[15:42, I]))
            self.assertTrue(np.allclose(F[I,::2].toarray(),
                                        F.toarray()[I,::2]))
            self.assertTrue(np.allclose(F[::2, J].toarray(),
                                        F.toarray()[::2, J]))
            self.assertTrue(np.allclose(F.T[J,::2].toarray(),
                                        F.T.toarray()[J,::2]))
            self.assertTrue(np.allclose(F.T[::2, I].toarray(),
                                        F.T.toarray()[::2, I]))



    def testToDense(self):
        print("testToDense()")
        prod = self.mulFactors()
        test_prod = self.F.toarray()
        self.assertProdEq(prod, test_prod)
        #self.assertTrue((self.F.toarray() == prod).all())


    def testPlus(self):
        from pyfaust import rand as frand
        from numpy.random import rand
        print("testPlus()")
        print("addition of a Faust-scalar")
        scals = [ self.r.random()*100,
                 np.complex(self.r.random()*100,self.r.random()*100),
                 self.r.randint(1,100)]
        F = self.F
        for s in scals:
            print("scalar=", s)
            test_F = F+s
            ref_full_F = self.mulFactors()+s
            self.assertAlmostEqual(norm(test_F.toarray()-ref_full_F), 0, places=3)
        print("addition Faust-Faust")
        fausts = \
        [ frand(F.shape[0], F.shape[0])@Faust(rand(F.shape[0],F.shape[1])),
         frand(F.shape[0],F.shape[0], density=.5,
                           field='complex')@Faust(rand(F.shape[0],F.shape[1]))]
        for i in range(0,len(fausts)):
            F2 = fausts[i]
            self.assertAlmostEqual(norm((F+F2).toarray()-
                                        (F.toarray()+F2.toarray())), 0,
                                   places=2)

    def testMinus(self):
        from pyfaust import rand as frand
        from numpy.random import rand
        print("testMinus()")
        print("subtraction of a Faust-scalar")
        scals = [ self.r.random()*100,
                 np.complex(self.r.random()*100,self.r.random()*100),
                 self.r.randint(1,100)]
        F = self.F
        for s in scals:
            print("scalar=", s)
            test_F = F+s
            ref_full_F = self.mulFactors()+s
            self.assertAlmostEqual(norm(test_F.toarray()-ref_full_F), 0,
                                   places=3)
        print("subtraction Faust-Faust")
        fausts = \
        [ frand(F.shape[0], F.shape[0])@Faust(rand(F.shape[0],F.shape[1])),
         frand(F.shape[0], F.shape[0], density=.5,
                           field='complex')@Faust(rand(F.shape[0],F.shape[1]))]
        for i in range(0,len(fausts)):
            F2 = fausts[i]
            self.assertAlmostEqual(norm((F-F2).toarray()-
                                        (F.toarray()-F2.toarray())), 0,
                                   places=2)

    def testScalDiv(self):
        print("test div by a scalar: real and complex.")
        scals = [ self.r.random()*100,
                 np.complex(self.r.random()*100,self.r.random()*100),
                 self.r.randint(1,100)]
        F = self.F
        for s in scals:
            test_F = F/s
            ref_full_F = self.mulFactors()/s
            self.assertAlmostEqual(norm(test_F.toarray()-ref_full_F), 0, places=3)

    def test_prod_opt(self):
        from pyfaust import FaustMulMode
        print("test GREEDY and DYNPROG prod opt methods and optimize_time.")
        GREEDY = 4
        self.F.save('/tmp/F.mat')
        H = self.F.clone()
        H.m_faust.set_FM_mul_mode(GREEDY) # FaustMulMode.GREEDY replaced by GREEDY local variable because GREEDY is not a visible opt. method anymore
        G = self.F.clone()
        G.m_faust.set_FM_mul_mode(FaustMulMode.DYNPROG)
        self.assertTrue(np.allclose(self.F.toarray(), H.toarray()))
        self.assertTrue(np.allclose(self.F.toarray(), G.toarray()))
        M = np.random.rand(self.F.shape[1], self.F.shape[0]).astype(self.F.dtype)
        self.assertTrue(np.allclose(self.F@M, H@M))
        self.assertTrue(np.allclose(self.F@M, G@M))
        S = sparse.random(self.F.shape[1], self.F.shape[0], .2, format='csr').astype(self.F.dtype)
        self.assertTrue(np.allclose(self.F@S, H@S))
        self.assertTrue(np.allclose(self.F@S, G@S))
        # test any method chosen by optimize_time
        I = self.F.optimize_time()
        self.assertTrue(np.allclose(self.F.toarray(), I.toarray()))
        self.assertTrue(np.allclose(self.F@M, I@M))
        self.assertTrue(np.allclose(self.F@S, I@S))
        # using the F@M benchmark
        J = self.F.optimize_time(mat=M)
        self.assertTrue(np.allclose(self.F.toarray(), J.toarray()))
        self.assertTrue(np.allclose(self.F@M, J@M))
        self.assertTrue(np.allclose(self.F@S, J@S))
        # using the F@S benchmark
        K = self.F.optimize_time(mat=S)
        self.assertTrue(np.allclose(self.F.toarray(), K.toarray()))
        self.assertTrue(np.allclose(self.F@M, K@M))
        self.assertTrue(np.allclose(self.F@S, K@S))

    def testMul(self):
        print("testMul()")
        print("test mul by a full real matrix")
        rmat = np.random.rand(self.F.shape[1],
                              self.r.randint(1,TestFaustPy.MAX_DIM_SIZE))
        prod = self.mulFactors().dot(rmat)
        test_prod = self.F.dot(rmat)
        self.assertProdEq(prod, test_prod)
        print("test mul by a full complex matrix")
        j = np.complex(0,1)
        rand = np.random.rand
        cmat = rand(rmat.shape[0], rmat.shape[1]) + j*rand(rmat.shape[0],
                                                           rmat.shape[1])
        prod = self.mulFactors().dot(cmat)
        test_prod = self.F.dot(cmat)
        self.assertProdEq(prod, test_prod)
        print("test mul by a real scalar")
        import random
        r = random.random()*100
        test_prod = self.F*r
        ref_prod = self.mulFactors()*r
        self.assertLess(norm(test_prod.toarray()-ref_prod)/norm(ref_prod),
                        1**-5)
        self.assertLess(norm((self.F.T*r).toarray()-self.F.toarray().T*r)/norm(self.F.toarray().T*r),1**-5)
        self.assertLess(norm((self.F.H*r).toarray()-self.F.toarray().T.conj()*r)/norm(self.F.toarray().T.conj()*r),1**-5)

        print("test mul by a complex scalar")
        c = np.complex(r, random.random()*100)
        test_prod = self.F*c
        ref_prod = self.mulFactors()*c
        self.assertLess(norm(test_prod.toarray()-ref_prod)/norm(ref_prod),
                        1**-5)
        print("test mul of two Fausts")
        from pyfaust import rand as frand, Faust
        F = self.F
        r_fausts = [ frand(F.shape[1], F.shape[1], self.r.randint(1,10)),
                    frand(F.shape[1], F.shape[1], self.r.randint(1,10), field='complex')]
        for i in range(0,len(r_fausts)):
            rF = r_fausts[i]
            assert(isinstance(rF, Faust))
            test_prod = F@rF
            ref_prod = self.mulFactors().dot(rF.toarray())
#            print('test_prod=', test_prod)
#            print('ref_prof=', ref_prod.shape)
#            print("test_prod=", test_prod.toarray())
#            print("ref_prof=", ref_prod)
            self.assertLess(norm(test_prod.toarray()-ref_prod)/norm(ref_prod),
                            1**-5)

        from scipy.sparse import dia_matrix, csr_matrix
        print("test mul of a Faust by a dia_matrix")
        D = dia_matrix((np.random.rand(1,F.shape[1]),np.array([0])),
                       shape=(F.shape[1],F.shape[1]))
        test_prod = F@D
        self.assertTrue(np.allclose(test_prod, F.toarray().dot(D.toarray())))
        print("test mul of a Faust by a complex dia_matrix")
        D = \
        dia_matrix((np.random.rand(1,F.shape[1])+np.random.rand(1,F.shape[1])*np.complex(0,1),np.array([0])),
                       shape=(F.shape[1],F.shape[1]))
        test_prod = F@D
        self.assertTrue(np.allclose(test_prod, F.toarray().dot(D.toarray())))
        Mr = csr_matrix(rand(F.shape[1],10))
        Mc = csr_matrix(rand(F.shape[1],10)+np.complex(0,1)*rand(F.shape[1],10))
        print("test mul Faust-csr_matrix")
        self.assertTrue(np.allclose(F@Mr, F.toarray().dot(Mr.toarray())))
        print("test mul Faust-complex csr_matrix")
        self.assertTrue(np.allclose(F@Mc, F.toarray().dot(Mc.toarray())))

    def testConcatenate(self):
        print("testConcatenate()")
        from pyfaust import rand
        from numpy.linalg import norm
        F = self.F
        FAUST,SPARSE,FULL=0,1,2
        for typeH in range(0,3):
            for cat_axis in [0,1]:
                if(isinstance(self,TestFaustPyCplx)):
                    field = 'complex'
                else:
                    field = 'real'
                G = rand(F.shape[(cat_axis+1)%2], F.shape[(cat_axis+1)%2], self.r.randint(1,TestFaustPy.MAX_NUM_FACTORS), field=field)
                # add one random factor to get a random number of rows to test
                # vertcat and a random number of cols to test horzcat
                if cat_axis == 0:
                    M = sparse.csr_matrix(np.random.rand(
                        self.r.randint(1,TestFaustPy.MAX_DIM_SIZE),
                        F.shape[(cat_axis+1)%2]).astype(F.dtype))
                    H = Faust([M]+[G.factors(i) for i in
                           range(0,G.numfactors())])
                else:
                    M = sparse.csr_matrix(np.random.rand(
                        F.shape[(cat_axis+1)%2],
                        self.r.randint(1,TestFaustPy.MAX_DIM_SIZE)).astype(F.dtype))
                    H = Faust([G.factors(i) for i in
                           range(0,G.numfactors())]+[M])
                if(typeH == FAUST):
                    H_ = H
                elif(typeH == SPARSE):
                    from scipy.sparse import csr_matrix
                    H_ = csr_matrix(H.toarray())
                else: # typeH == FULL
                    H_ = H.toarray()
                print("testConcatenate() F.shape, H.shape", F.shape, H.shape)
                C = F.concatenate(H_,axis=cat_axis)
                ref_C = np.concatenate((F.toarray(),
                                        H.toarray()),
                                       axis=cat_axis)
                self.assertEqual(C.shape, ref_C.shape)
                self.assertLessEqual(norm(C.toarray()-ref_C)/norm(ref_C),
                            10**-5)
        # test random number of Fausts concatenation
        from pyfaust import rand as frand, concatenate as cat, isFaust
        for cat_axis in [0, 1]:
            n = self.r.randint(3, 18)
            fausts = []
            field_names = ['real', 'complex']
            fac_types = ['sparse', 'dense', 'mixed']
            for i in range(n):
                fac_type_id = self.r.randint(0,2)
                field_id = self.r.randint(0,1)
                is_faust  = bool(self.r.randint(0,1)) or i == 0
                nrows = self.r.randint(2, 128) if cat_axis == 0 else F.shape[0]
                ncols = self.r.randint(2, 128) if cat_axis == 1 else F.shape[1]
                if is_faust:
                    fausts += [frand(nrows, ncols,
                                     fac_type=fac_types[fac_type_id],
                                     field=field_names[field_id])]
                else:
                    fausts += [frand(nrows, ncols, fac_type=fac_types[fac_type_id],
                                     field=field_names[field_id],
                                     num_factors=1).factors(0)]
            Fc = cat(tuple(fausts), axis=cat_axis)
            arrays = []
            for F in fausts:
                if not isinstance(F, np.ndarray):
                    arrays += [F.toarray()]
                else:
                    arrays += [F]
            Mc = np.concatenate(arrays, axis=cat_axis)
            self.assertTrue(np.allclose(Fc.toarray(), Mc))





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

    def testShape(self):
        print("testShape()")
        self.assertEqual((self.F.shape[0],self.F.shape[1]), self.F.shape)

    def testSize(self):
        print("testSize()")
        test_size = self.F.size
        ref_size = self.mulFactors().size
        self.assertEqual(test_size,ref_size)

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

    def test_real(self):
        print("test Faust.real")
        rF = self.F.real
        rF_ref = self.mulFactors().real
        self.assertTrue(np.allclose(rF.toarray(), rF_ref))

    def test_imag(self):
        print("test Faust.imag")
        iF = self.F.imag
        iF_ref = self.mulFactors().imag
        self.assertTrue(np.allclose(iF.toarray(), iF_ref))

    def test_left(self):
        print("Test Faust.left()")
        for F in [self.F, self.F.T, self.F.H]:
            for i in range(len(F)):
                lFt = F.left(i)
                if i == 0:
                    lFt = Faust([lFt])
                lFr = Faust([F.factors(j) for j in range(i+1)])
                self.assertEqual(len(lFr), len(lFt))
                for fid in range(i+1):
                    a = lFt.factors(fid)
                    b = lFr.factors(fid)
                    if(not isinstance(a, np.ndarray)):
                        a = a.toarray()
                    if(not isinstance(b, np.ndarray)):
                        b = b.toarray()
                    self.assertTrue(np.allclose(a,b))

    def test_right(self):
        print("Test Faust.right()")
        for F in [self.F, self.F.T, self.F.H]:
            for i in range(len(F)):
                rFt = F.right(i)
                if i == len(F)-1:
                    rFt = Faust([rFt])
                rFr = Faust([F.factors(j) for j in range(i, len(F))])
                self.assertEqual(len(rFr), len(rFt))
                for fid in range(len(rFr)):
                    a = rFt.factors(fid)
                    b = rFr.factors(fid)
                    if(not isinstance(a, np.ndarray)):
                        a = a.toarray()
                    if(not isinstance(b, np.ndarray)):
                        b = b.toarray()
                    self.assertTrue(np.allclose(a,b))

    def test_circulant(self):
        v = np.random.rand(1024)
        from scipy.linalg import circulant
        from pyfaust import circ
        self.assertTrue(np.allclose(circulant(v), circ(v).toarray()))

    def test_anticirculant(self):
        v = np.random.rand(1024)
        from scipy.linalg import circulant
        from pyfaust import anticirc
        P = np.zeros((len(v), len(v)))
        I = np.arange(len(v)-1, -1, -1)
        J = np.arange(0, len(v))
        P[I, J] = 1
        self.assertTrue(np.allclose(circulant(v), anticirc(v).toarray()@P))

    def test_toeplitz(self):
        r = np.random.rand(1024)
        c = np.random.rand(2048)
        from scipy.linalg import toeplitz
        from pyfaust import toeplitz as ftoeplitz
        self.assertTrue(np.allclose(toeplitz(c, r), ftoeplitz(c, r).toarray()))

    def test_bsr_get_fact(self):
        print("test_bsr_get_fact")
        from pyfaust import rand_bsr, Faust
        from numpy import allclose
        nfacs = 5
        F = rand_bsr(15, 15, 3, 5, nfacs)
        for i in range(nfacs):
            bf = F.factors(i)
            G = Faust(bf)
            bf2  = G.factors(0)
            self.assertTrue(allclose(bf.toarray(), bf2.toarray()))

    def test_palm4msa_constraints_consistency_checking(self):
        print("Test ParamsPalm4MSA.are_constraints_consistent")
        from pyfaust.factparams import ParamsPalm4MSA, ConstraintList, StoppingCriterion
        from pyfaust.proj import splin, normcol
        ##### error test 1: last constraint number of columns vs matrix number of columns
        M = np.random.rand(500, 32).astype(self.F.dtype)
        cons = ConstraintList('splin', 5, 500, 32, 'normcol', 1.0, 32, 33)
        projs = [splin((500,32), 5), normcol((32,33), 1.0)] # the same with projs
        stop_crit = StoppingCriterion(num_its=200)
        # WARNING: for assertRaisesRegex it must be on one line or only the first line
        # is tested as regex (is that a bug?)
        err_msg = "ParamsPalm4MSA error: the matrix M to factorize must have the same number of columns as the last constraint. They are respectively: 32 and 33."
        param = ParamsPalm4MSA(cons, stop_crit)
        self.assertRaisesRegex(ValueError, err_msg, param.are_constraints_consistent, M)
        param = ParamsPalm4MSA(projs, stop_crit)
        self.assertRaisesRegex(ValueError, err_msg, param.are_constraints_consistent, M)
        ##### error test 2: first constraint number of rows vs matrix number of rows
        cons = ConstraintList('splin', 5, 501, 32, 'normcol', 1.0, 32, 32)
        projs = [splin((501, 32), 5), normcol((32, 32), 1.0)] # the same with projs
        stop_crit = StoppingCriterion(num_its=200)
        # cf. WARNING above
        err_msg = "ParamsPalm4MSA error: the matrix M to factorize must have the same number of rows as the first constraint. They are respectively: 500 and 501."
        param = ParamsPalm4MSA(cons, stop_crit)
        self.assertRaisesRegex(ValueError, err_msg, param.are_constraints_consistent, M)
        param = ParamsPalm4MSA(projs, stop_crit)
        self.assertRaisesRegex(ValueError, err_msg, param.are_constraints_consistent, M)
        ##### error test 3: the constraints dimensions must follow a well-defined matrix product
        cons = ConstraintList('splin', 5, 500, 33, 'normcol', 1.0, 32, 32)
        projs = [splin((500, 33), 5), normcol((32, 32), 1.0)] # the same with projs
        stop_crit = StoppingCriterion(num_its=200)
        # cf. WARNING above
        err_msg = "The 1-index constraint number of rows \(which is 32\) must be equal to the 0-index constraint number of columns \(which is 33\)."
        param = ParamsPalm4MSA(cons, stop_crit)
        self.assertRaisesRegex(ValueError, err_msg, param.are_constraints_consistent, M)
        param = ParamsPalm4MSA(projs, stop_crit)
        self.assertRaisesRegex(ValueError, err_msg, param.are_constraints_consistent, M)
        ##### sucess test: the constraints dimensions must follow a well-defined matrix product
        cons = ConstraintList('splin', 5, 500, 32, 'normcol', 1.0, 32, 32)
        projs = [splin((500, 32), 5), normcol((32, 32), 1.0)] # the same with projs
        stop_crit = StoppingCriterion(num_its=200)
        # cf. WARNING above
        param = ParamsPalm4MSA(cons, stop_crit)
        self.assertTrue(param.are_constraints_consistent(M))
        param = ParamsPalm4MSA(projs, stop_crit)
        self.assertTrue(param.are_constraints_consistent(M))

    def test_hierarchical_constraints_consistency_checking(self):
        print("ParamsHierarchical.are_constraints_consistent")
        from pyfaust.factparams import (ParamsHierarchical, ConstraintList,
                                        StoppingCriterion)
        from pyfaust.proj import splin, sp, normcol
        M = np.random.rand(500, 32)
        stop_crit1 = StoppingCriterion(num_its=200)
        stop_crit2 = StoppingCriterion(num_its=200)
        ### 1st test: it must work, constraints/projs are valid
        fact_cons = ConstraintList('splin', 5, 500, 32, 'sp', 96, 32, 32,
                                       'sp', 96, 32, 32)
        res_cons = ConstraintList('normcol', 1, 32, 32, 'sp', 666, 32, 32,
                                      'sp', 333, 32, 32)
        param = ParamsHierarchical(fact_cons, res_cons, stop_crit1,
                                       stop_crit2)
        self.assertTrue(param.are_constraints_consistent(M))
        # test the same with projectors
        fact_projs = [splin((500, 32), 5), sp((32,32), 96), sp((32,32), 96)]
        res_projs = [normcol((32,32), 1), sp((32,32), 666), sp((32,32), 333)]
        param = ParamsHierarchical(fact_projs, res_projs, stop_crit1,
                                       stop_crit2)
        self.assertTrue(param.are_constraints_consistent(M))
        ### 2nd test: erroneous end constraint sizes (which must be equal to M size)
        fact_cons = ConstraintList('splin', 5, 501, 32, 'sp', 96, 32, 32,
                                   'sp', 96, 32, 32)
        res_cons = ConstraintList('normcol', 1, 32, 32, 'sp', 666, 32, 32,
                                  'sp', 333, 32, 32)

        param = ParamsHierarchical(fact_cons, res_cons, stop_crit1,
                                   stop_crit2)
        err_msg = "The number of rows of the 0-index factor constraint \(501\) must be equal to the number of rows of the matrix to factorize \(500\) \(is_fact_side_left=False\)."
        self.assertRaisesRegex(ValueError, err_msg, param.are_constraints_consistent, M)
        fact_projs = [splin((501, 32), 5), sp((32,32), 96), sp((32,32), 96)]
        res_projs = [normcol((32,32), 1), sp((32,32), 666), sp((32,32), 333)]
        param = ParamsHierarchical(fact_cons, res_cons, stop_crit1,
                                   stop_crit2)
        self.assertRaisesRegex(ValueError, err_msg, param.are_constraints_consistent, M)
        ### 3rd test: constraint intermediary inconsistent sizes
        fact_cons = ConstraintList('splin', 5, 500, 31, 'sp', 96, 32, 32,
                                   'sp', 96, 32, 32)
        res_cons = ConstraintList('normcol', 1, 32, 32, 'sp', 666, 32, 32,
                                  'sp', 333, 32, 32)

        param = ParamsHierarchical(fact_cons, res_cons, stop_crit1,
                                   stop_crit2)
        err_msg = "The number of columns \(31\) of the 0-index factor constraint must be equal to the number of rows \(32\) of the 0-index residual constraint \(is_fact_side_left=False\)."
        param = ParamsHierarchical(fact_cons, res_cons, stop_crit1,
                                   stop_crit2)
        self.assertRaisesRegex(ValueError, err_msg, param.are_constraints_consistent, M)
        fact_cons = ConstraintList('splin', 5, 500, 31, 'sp', 96, 32, 32,
                                   'sp', 96, 32, 32)
        res_cons = ConstraintList('normcol', 1, 31, 32, 'sp', 666, 32, 32,
                                  'sp', 333, 32, 32)

        param = ParamsHierarchical(fact_cons, res_cons, stop_crit1,
                                   stop_crit2)
        err_msg = "The number of rows \(32\) of the 1-index factor constraint must be equal to the number of columns \(31\) of the 0-index factor constraint \(is_fact_side_left=False\)."
        param = ParamsHierarchical(fact_cons, res_cons, stop_crit1,
                                   stop_crit2)
        self.assertRaisesRegex(ValueError, err_msg, param.are_constraints_consistent, M)
        fact_cons = ConstraintList('splin', 5, 500, 32, 'sp', 96, 32, 32,
                                   'sp', 96, 32, 32)
        res_cons = ConstraintList('normcol', 1, 32, 32, 'sp', 666, 64, 32,
                                  'sp', 333, 32, 32)

        param = ParamsHierarchical(fact_cons, res_cons, stop_crit1,
                                   stop_crit2)
        err_msg = "The number of columns \(32\) of the 1-index factor constraint must be equal to the number of rows \(64\) of the 1-index residual constraint \(is_fact_side_left=False\)."
        param = ParamsHierarchical(fact_cons, res_cons, stop_crit1,
                                   stop_crit2)
        self.assertRaisesRegex(ValueError, err_msg, param.are_constraints_consistent, M)

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
                                       density=min(r.random()+.5,1)).astype(np.complex)]
                    #print("setUp() i=",i, "d1=",d1, "d2=", d2, "factor.shape[0]=",
                    #      factors[i].shape[0])
                    factors[i] *= np.complex(r.random(), r.random())*5 # 5 is
                    # arbitrary
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

class TestFaustFactory(unittest.TestCase):

    def testFactPalm4MSA(self):
        from pyfaust.fact import palm4msa
        print("Test pyfaust.fact.palm4msa()")
        from pyfaust.factparams import ConstraintReal,\
                ConstraintInt, ConstraintName
        from pyfaust.factparams import ParamsPalm4MSA, StoppingCriterion
        #num_facts = 2
        #is_update_way_R2L = False
        #init_lambda = 1.0
#        init_facts = list()
#        init_facts.append(np.zeros([500,32]))
#        init_facts.append(np.eye(32))
        #M = np.random.rand(500, 32)
        M = \
                loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/config_compared_palm2.mat")['data']
        # default step_size
        cons1 = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5)
        cons2 = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32,
                                 32, 1.0)
        stop_crit = StoppingCriterion(num_its=200)
        param = ParamsPalm4MSA([cons1, cons2], stop_crit, init_facts=None,
                               is_update_way_R2L=False, init_lambda=1.0,
                               is_verbose=False, constant_step_size=False)
        F, _lambda = palm4msa(M, param, ret_lambda=True)
        #F.display()
        #print("normF", F.norm("fro"))
        self.assertEqual(F.shape, M.shape)
        E = F.toarray()-M
        #print("err.:",norm(F.toarray(), "fro"),  norm(E,"fro"), norm (M,"fro"))
        print("err:", norm(E,"fro")/norm(M,"fro"))
        print("_lambda:", _lambda)
        # matrix to factorize and reference relative error come from
        # misc/test/src/C++/test_palm4MSA.cpp
        self.assertAlmostEqual(norm(E,"fro")/norm(M,"fro"), 0.270954109668, places=4)

    def testParamsPalm4MSA(self):
        # call Palm4MSA specifying params
        from os import dup2, pipe # for
        from pyfaust.fact import palm4msa
        from pyfaust.factparams import (ParamsPalm4MSA, ConstraintList,
                                        StoppingCriterion,
                                        ConstraintInt,
                                        ConstraintReal, ParamsFact)
        import numpy as np
        from tempfile import gettempdir
        from os.path import join
        M = np.random.rand(500, 32) 
        cons = ConstraintList('splin', 5, 500, 32, 'normcol', 1.0, 32, 32) 
        # or alternatively using pyfaust.proj
        # from pyfaust.proj import splin, normcol
        # cons = [ splin((500,32), 5), normcol((32,32), 1.0)]
        stop_crit = StoppingCriterion(num_its=200)
        param = ParamsPalm4MSA(cons, stop_crit)
        param.is_verbose = True
        param.grad_calc_opt_mode = 1 
        param.factor_format = 'dense'
        param.packing_RL = False
        tmp_dir = gettempdir()
        tmp_file = join(tmp_dir, "verbose_output_of_palm4msa_test")
        print("tmp_file:", tmp_file)
        f = open(tmp_file, 'w')
        dup2(1,2)
        dup2(f.fileno(), 1)
        F = palm4msa(M, param)
        print()
        f.close()
        dup2(2,1)
        # retrieve the params effectively used from C++ core output
        # reconstruct a ParamsPalm4MSA from the values found
        param_test = ParamsPalm4MSA(cons, stop_crit)
        param_test.constraints = []
        for line in open(tmp_file, 'r').readlines():
            print(line, end='')
            if(line.startswith('NFACTS')):
                param_test.num_facts = int(line.split(':')[-1].strip())
            if(line.startswith('VERBOSE')):
                param_test.is_verbose = bool(int(line.split(':')[-1].strip()))
            if(line.startswith('UPDATEWAY')):
                param_test.is_update_way_R2L = \
                bool(int(line.split(':')[-1].strip()))
            if(line.startswith('INIT_LAMBDA')):
                param_test.init_lambda = float(line.split(':')[-1].strip())
            if(line.startswith('ISCONSTANTSTEPSIZE')):
                param_test.constant_step_size = \
                bool(int(line.split(':')[-1].strip()))
            if(line.startswith('step_size')):
                param_test.step_size = float(line.split(':')[-1].strip())
            if(line.startswith('gradCalcOptMode')):
                param_test.grad_calc_opt_mode = int(line.split(':')[-1].strip())
            if(line.startswith('factors format')):
                param_test.factor_format = int(line.split(':')[-1].strip())
                param_test.factor_format = \
                ParamsFact.factor_format_int2str(param_test.factor_format)
            if(line.startswith('packing_RL')):
                param_test.packing_RL = int(line.split(':')[-1].strip()) != 0
                print(param_test.packing_RL)
            if(line.startswith('errorThreshold')):
                param_test.stop_crit.tol = float(line.split(':')[-1].strip())
                param_test.stop_crit._is_criterion_error = True
            if(line.startswith('nb_it')):
                param_test.stop_crit.num_its = int(line.split(':')[-1].strip())
            if(line.startswith('maxIteration')):
                param_test.stop_crit.maxiter = \
                bool(int(line.split(':')[-1].strip()))
            if(line.startswith('type_cont')):
                colon_fields = line.split(':')
                const_type_name = colon_fields[1].strip()
                if(const_type_name.startswith('FAUST_INT CONSTRAINT_NAME_SPLIN')):
                    nrows = int(colon_fields[2].split(' ')[0].strip())
                    ncols = int(colon_fields[3].split(' ')[0].strip())
                    cons_val = int(colon_fields[-1].strip())
                    cons = ConstraintInt('splin', nrows, ncols, cons_val)
                    param_test.constraints += [cons]
                if(const_type_name.startswith('FAUST_REAL CONSTRAINT_NAME_NORMCOL')):
                    nrows = int(colon_fields[2].split(' ')[0].strip())
                    ncols = int(colon_fields[3].split(' ')[0].strip())
                    cons_val = float(colon_fields[-1].strip())
                    cons = ConstraintReal('normcol', nrows, ncols, cons_val)
                    param_test.constraints += [cons]
        # compare original param instance and the reconstructed one
        self.assertEqual(param_test.num_facts, param.num_facts)
        self.assertEqual(param_test.is_verbose, param.is_verbose)
        self.assertEqual(param_test.is_update_way_R2L, param.is_update_way_R2L)
        self.assertEqual(param_test.init_lambda, param.init_lambda)
        self.assertEqual(param_test.constant_step_size,
                         param.constant_step_size)
        self.assertEqual(param_test.step_size, param.step_size)
        self.assertEqual(param_test.grad_calc_opt_mode, param.grad_calc_opt_mode)
#        self.assertEqual(param_test.factor_format, param.factor_format)
#        self.assertEqual(param_test.packing_RL, param.packing_RL)
        self.assertEqual(param_test.stop_crit.maxiter, param.stop_crit.maxiter)
        self.assertEqual(param_test.stop_crit._is_criterion_error,
                         param.stop_crit._is_criterion_error)
        self.assertEqual(param_test.stop_crit.num_its,
                         param.stop_crit.num_its)
        self.assertEqual(len(param_test.constraints), len(param.constraints))
        for i in range(len(param_test.constraints)):
            self.assertEqual(param_test.constraints[i]._num_rows,
                             (param.constraints[i]._num_rows))
            self.assertEqual(param_test.constraints[i]._num_cols,
                             (param.constraints[i]._num_cols))
            self.assertEqual(param_test.constraints[i].name,
                             param_test.constraints[i].name)
            self.assertEqual(param_test.constraints[i]._cons_value,
                             param_test.constraints[i]._cons_value)
        os.remove(tmp_file)

    def testFactPalm4MSA2020(self):
        from pyfaust.fact import palm4msa
        print("Test pyfaust.fact.palm4msa2020()")
        from pyfaust.factparams import ConstraintReal,\
                ConstraintInt, ConstraintName
        from pyfaust.factparams import ParamsPalm4MSA, StoppingCriterion
        #num_facts = 2
        #is_update_way_R2L = False
        #init_lambda = 1.0
#        init_facts = list()
#        init_facts.append(np.zeros([500,32]))
#        init_facts.append(np.eye(32))
        #M = np.random.rand(500, 32)
        M = \
                loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/config_compared_palm2.mat")['data']
        # default step_size
        cons1 = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5)
        cons2 = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32,
                                 32, 1.0)
        stop_crit = StoppingCriterion(num_its=200)
        param = ParamsPalm4MSA([cons1, cons2], stop_crit, init_facts=None,
                               is_update_way_R2L=False, init_lambda=1.0,
                               is_verbose=False, constant_step_size=False)
        F, _lambda = palm4msa(M, param, ret_lambda=True, backend=2020)
        #F.display()
        #print("normF", F.norm("fro"))
        self.assertEqual(F.shape, M.shape)
        E = F.toarray()-M
        #print("err.:",norm(F.toarray(), "fro"),  norm(E,"fro"), norm (M,"fro"))
        print("err:", norm(E,"fro")/norm(M,"fro"))
        print("_lambda:", _lambda)
        # matrix to factorize and reference relative error come from
        # misc/test/src/C++/test_palm4MSA.cpp
        self.assertAlmostEqual(norm(E,"fro")/norm(M,"fro"), 0.270954109668, places=4)

    def testFactHierarch(self):
        print("Test pyfaust.fact.hierarchical()")
        from pyfaust.fact import hierarchical
        from pyfaust.factparams import ParamsHierarchical, StoppingCriterion
        from pyfaust.factparams import ConstraintReal, ConstraintInt,\
                ConstraintName
        num_facts = 4
        is_update_way_R2L = False
        init_lambda = 1.0
        #M = np.random.rand(500, 32)
        M = \
        loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/matrix_hierarchical_fact.mat")['matrix']
        # default step_size
        fact0_cons = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5)
        fact1_cons = ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 96)
        fact2_cons = ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 96)
        res0_cons = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32, 32, 1)
        res1_cons =  ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 666)
        res2_cons =  ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 333)
        stop_crit1 = StoppingCriterion(num_its=200)
        stop_crit2 = StoppingCriterion(num_its=200)
        param = ParamsHierarchical([fact0_cons, fact1_cons, fact2_cons],
                                       [res0_cons, res1_cons, res2_cons],
                                       stop_crit1, stop_crit2,
                                       is_verbose=False)
        F = hierarchical(M, param, backend=2020)
        self.assertEqual(F.shape, M.shape)
        E = F.toarray()-M
        #print("err.:",norm(F.toarray(), "fro"),  norm(E,"fro"), norm (M,"fro"),
        print("err: ", norm(E,"fro")/norm(M,"fro"))
        # matrix to factorize and reference relative error come from
        # misc/test/src/C++/hierarchicalFactorization.cpp
        self.assertAlmostEqual(norm(E,"fro")/norm(M,"fro"),
                               0.267959, places=5)

    def testFactHierarchCplx(self):
        print("Test pyfaust.fact.hierarchicalCplx()")
        from pyfaust.fact import hierarchical
        from pyfaust.factparams import ParamsHierarchical, StoppingCriterion
        from pyfaust.factparams import ConstraintReal,\
                ConstraintInt, ConstraintName
        num_facts = 4
        is_update_way_R2L = False
        init_lambda = 1.0
        #M = np.random.rand(500, 32)
        M = \
        loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/matrix_hierarchical_fact.mat")['matrix']
        M = M + np.complex(0,1)*M
        # default step_size
        fact0_cons = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5)
        fact1_cons = ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 96)
        fact2_cons = ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 96)
        res0_cons = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32, 32, 1)
        res1_cons =  ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 666)
        res2_cons =  ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 333)
        stop_crit1 = StoppingCriterion(num_its=200)
        stop_crit2 = StoppingCriterion(num_its=200)
        param = ParamsHierarchical([fact0_cons, fact1_cons, fact2_cons],
                                       [res0_cons, res1_cons, res2_cons],
                                       stop_crit1, stop_crit2,
                                       is_verbose=False)
        F = hierarchical(M, param)
        self.assertEqual(F.shape, M.shape)
        #print(F.toarray())
        E = F.toarray()-M
        #print("err.:",norm(F.toarray(), "fro"),  norm(E,"fro"), norm (M,"fro"),
        print("err: ", norm(E,"fro")/norm(M,"fro"))
        # matrix to factorize and reference relative error come from
        # misc/test/src/C++/hierarchicalFactorization.cpp
        self.assertAlmostEqual(norm(E,"fro")/norm(M,"fro"), 0.273 , places=3)

    def testFactPalm4MSACplx(self):
        print("Test pyfaust.fact.palm4msaCplx()")
        from pyfaust.fact import palm4msa
        from pyfaust.factparams import ParamsPalm4MSA,StoppingCriterion
        from pyfaust.factparams import ConstraintReal,\
                ConstraintInt, ConstraintName
        num_facts = 2
        is_update_way_R2L = False
        init_lambda = 1.0
#        init_facts = list()
#        init_facts.append(np.zeros([500,32]))
#        init_facts.append(np.eye(32))
        #M = np.random.rand(500, 32)
        M = \
        loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/config_compared_palm2.mat")['data']
        M = M + np.complex(0,1)*M
        # default step_size
        cons1 = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5)
        cons2 = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32,
                                 32, 1.0)
        stop_crit = StoppingCriterion(num_its=200)
        param = ParamsPalm4MSA([cons1, cons2], stop_crit, init_facts=None,
                               is_verbose=False, constant_step_size=False,
                               is_update_way_R2L=False, init_lambda=1.0)
        F,_lambda = palm4msa(M, param, ret_lambda=True)
        #F.display()
        #print("normF", F.norm("fro"))
        self.assertEqual(F.shape, M.shape)
        #print(F.toarray())
        E = F.toarray()-M
        #print("err.:",norm(F.toarray(), "fro"),  norm(E,"fro"), norm (M,"fro"))
        print("err:", norm(E,"fro")/norm(M,"fro"))
        print("lambda:", _lambda)
        # matrix to factorize and reference relative error come from
        # misc/test/src/C++/test_palm4MSA.cpp
        self.assertAlmostEqual(norm(E,"fro")/norm(M,"fro"), 0.272814, places=4)

    def testHadamard(self):
        print("Test pyfaust.wht()")
        from pyfaust import wht
        pow2_exp = random.Random().randint(1,10)
        n = 2**pow2_exp
        H = wht(n, False)
        fH = H.toarray()
        self.assertEqual(np.count_nonzero(fH), fH.size)
        for i in range(0,n-1):
            for j in range(i+1,n):
                self.assertTrue((fH[i,::].dot(fH[j,::].T) == 0).all())
        assert(np.allclose(wht(n).toarray(),
               wht(n, False).normalize().toarray()))

    def testFourier(self):
        print("Test pyfaust.dft()")
        from pyfaust import dft
        from numpy.fft import fft
        pow2_exp = random.Random().randint(1,10)
        n = 2**pow2_exp
        F = dft(n, False)
        fF = F.toarray()
        ref_fft = fft(np.eye(n))
        self.assertAlmostEqual(norm(ref_fft-fF)/norm(ref_fft),0)
        assert(np.allclose(dft(n).toarray(),
                           dft(n, False).normalize().toarray()))

    def testFourierDiagOpt(self):
        print("Test pyfaust.dft(diag_opt=True)")
        from pyfaust import dft
        pow2_exp = random.Random().randint(1,10)
        n = 2**pow2_exp
        for normed in [True, False]:
            F = dft(n, normed=normed, diag_opt=False)
            oF = dft(n, normed=normed, diag_opt=True)
            fF = oF.toarray()
            ref_fft = F.toarray()
            self.assertAlmostEqual(norm(ref_fft-fF)/norm(ref_fft),0)

    def testRandButterfly(self):
        print("Test pyfaust.rand_butterfly")
        from pyfaust import wht, rand_butterfly
        H = wht(32).toarray()
        for dtype in ['float32', 'double', 'complex']:
            F = rand_butterfly(32, dtype=dtype)
            self.assertTrue(not np.allclose(F.toarray(), H))
            ref_I, ref_J = np.nonzero(H)
            I, J, = np.nonzero(F.toarray())
            self.assertTrue((I == ref_I).all())
            self.assertTrue((J == ref_J).all())
            self.assertTrue(F.dtype == np.dtype(dtype))

    def testDCT(self):
        print("Test pyfaust.dct()")
        from pyfaust import dct
        from scipy.fft import dct as sdct
        from numpy.random import rand
        n = 512
        DCT = dct(n, normed=False)
        x = rand(n)
        y1 = DCT@x
        y2 = sdct(x)
        self.assertTrue(np.allclose(y1, y2))
        DCTn = dct(n, normed=True)
        self.assertTrue(np.allclose(DCTn@x, DCT.normalize()@x))

    def testDST(self):
        print("Test pyfaust.dst()")
        from pyfaust import dst
        from scipy.fft import dst as sdst
        from numpy.random import rand
        n = 512
        DST = dst(n, normed=False)
        x = rand(n)
        y1 = DST@x
        y2 = sdst(x)
        self.assertTrue(np.allclose(y1, y2))
        DSTn = dst(n, normed=True)
        self.assertTrue(np.allclose(DSTn@x, DST.normalize()@x))

    def testFGFTGivens(self):
        print("Test fact.fgft_givens()")
        print(sys.path)
        import pyfaust.fact
        L = loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat")['Lap']
        L = L.astype(np.float64)
        J = \
        int(loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat")['J'])
        D, F = pyfaust.fact.eigtj(L, J, nGivens_per_fac=1, verbosity=0, enable_large_Faust=True)
        D = spdiags(D, [0], L.shape[0], L.shape[0])
        print("Lap norm:", norm(L, 'fro'))
        err = norm((F@D.toarray())@F.T.toarray()-L,"fro")/norm(L,"fro")
        print("err: ", err)
        # the error reference is from the C++ test,
        # misc/test/src/C++/GivensFGFT.cpp.in
        self.assertAlmostEqual(err, 0.0326529, places=7)
        # check nnz for the F.numfactors()-1 first factors
        for fac in [F.factors(i) for i in range(0,F.numfactors())]:
            nnz = np.count_nonzero(fac.toarray())
            self.assertEqual(nnz, 2+L.shape[0])

    def testFGFTGivensParallel(self):
        print("Test pyfaust.fact.fgft_givens() -- parallel")
        from pyfaust.fact import eigtj
        L = loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat")['Lap']
        L = L.astype(np.float64)
        J = \
                int(loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat")['J'])
        t = int(L.shape[0]/2)
        D, F = eigtj(L, J, nGivens_per_fac=t, verbosity=0,
                           enable_large_Faust=True)
        D = spdiags(D, [0], L.shape[0], L.shape[0])
        print("Lap norm:", norm(L, 'fro'))
        err = norm((F@D.toarray())@F.T.toarray()-L,"fro")/norm(L,"fro")
        print("err: ", err)
        # the error reference is from the C++ test,
        # misc/test/src/C++/GivensFGFTParallel.cpp.in (_double version)
        self.assertAlmostEqual(err, 0.0398154, places=7)
        D2, F2 = eigtj(L, J, nGivens_per_fac=t, verbosity=0,
                       enable_large_Faust=True)
        D2 = spdiags(D, [0], L.shape[0], L.shape[0])
        print("Lap norm:", norm(L, 'fro'))
        err2 = norm((F2@D.toarray())@F2.T.toarray()-L,"fro")/norm(L,"fro")
        print("err2: ", err2)
        # the error reference is from the C++ test,
        # misc/test/src/C++/GivensFGFTParallel.cpp.in
        self.assertEqual(err, err2)
        # check nnz for the last factors (avoiding errors due to pivot image
        # near to zero in first factors)
        for fac in [F.factors(i) for i in range(F.numfactors()-10,F.numfactors())]:
            nnz = np.count_nonzero(fac.toarray())
            self.assertEqual(nnz, 2*t+L.shape[0])

    def testFGFTGivensSparse(self):
        print("Test fact.fgft_givens_sparse()")
        print(sys.path)
        from scipy.sparse import csr_matrix
        import pyfaust.fact
        L = loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat")['Lap']
        L = L.astype(np.float64)
        J = \
        int(loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat")['J'])
        D, F = pyfaust.fact.eigtj(csr_matrix(L), J, nGivens_per_fac=0,
                                        verbosity=0,enable_large_Faust=True)
        D = spdiags(D, [0], L.shape[0], L.shape[0])
        print("Lap norm:", norm(L, 'fro'))
        err = norm((F@D.toarray())@F.T.toarray()-L,"fro")/norm(L,"fro")
        print("err: ", err)
        # the error reference is from the C++ test,
        # misc/test/src/C++/GivensFGFT.cpp.in
        self.assertAlmostEqual(err, 0.0326529, places=7)
        # check nnz for the F.numfactors()-1 first factors
        for fac in [F.factors(i) for i in range(0,F.numfactors())]:
            nnz = np.count_nonzero(fac.toarray())
            self.assertEqual(nnz, 2+L.shape[0])

    def testFGFTGivensParallelSparse(self):
        print("Test pyfaust.fact.fgft_givens_sparse() -- parallel")
        from pyfaust.fact import eigtj
        from scipy.sparse import csr_matrix
        L = loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat")['Lap']
        L = L.astype(np.float64)
        J = \
                int(loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat")['J'])
        t = int(L.shape[0]/2)
        D, F = eigtj(csr_matrix(L), J, nGivens_per_fac=t, verbosity=0, enable_large_Faust=True)
        D = spdiags(D, [0], L.shape[0], L.shape[0])
        print("Lap norm:", norm(L, 'fro'))
        err = norm((F@D.toarray())@F.T.toarray()-L,"fro")/norm(L,"fro")
        print("err: ", err)
        # the error reference is from the C++ test,
        # misc/test/src/C++/GivensFGFTParallel.cpp.in (_double version)
        self.assertAlmostEqual(err, 0.0398154, places=7)
        D2, F2 = eigtj(csr_matrix(L), J, nGivens_per_fac=t, verbosity=0,
                       enable_large_Faust=True)
        D2 = spdiags(D, [0], L.shape[0], L.shape[0])
        print("Lap norm:", norm(L, 'fro'))
        err2 = norm((F2@D.toarray())@F2.T.toarray()-L,"fro")/norm(L,"fro")
        print("err2: ", err2)
        # the error reference is from the C++ test,
        # misc/test/src/C++/GivensFGFTParallel.cpp.in
        self.assertEqual(err, err2)
        # check nnz for the last factors
        for fac in [F2.factors(i) for i in range(F2.numfactors()-10, F2.numfactors())]:
            nnz = np.count_nonzero(fac.toarray())
            self.assertEqual(nnz, 2*t+L.shape[0])

    def testeigtj_abserr(self):
        from numpy.random import rand
        from pyfaust.fact import eigtj
        from scipy.io import loadmat
        import sys
        err = .01
        M = rand(128,128)
        M = M.dot(M.T)
        D,U = eigtj(M, tol=err, relerr=False)
        self.assertAlmostEqual(norm(M-U@np.diag(D)@U.T), err, places=3 )
        L = loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat")['Lap']
        L = L.astype(np.float64)
        M_cplx = L*np.complex(1,0) + L*np.complex(0,1)
        M_cplx = M_cplx.dot(np.matrix(M_cplx).H)
        err=.1
        D,U = eigtj(M_cplx, tol=err, relerr=False, verbosity=0)
        self.assertLessEqual(norm(M_cplx-U@np.diag(D)@U.H), err)

    def testeigtj_relerr(self):
        from numpy.random import rand
        from pyfaust.fact import eigtj
        err = .01
        M = rand(128,128)
        M = M.dot(M.T)
        D,U = eigtj(M, tol=err)
        self.assertLessEqual(norm(M-U@np.diag(D)@U.T)/norm(M), err)
        M_cplx = M + rand(128,128)*np.complex(0,1)
        M_cplx = M_cplx.dot(np.matrix(M_cplx).H)
        D,U = eigtj(M_cplx, tol=err)
        self.assertLessEqual(norm(M_cplx-U@np.diag(D)@U.H)/norm(M_cplx), err)

    def testFactPalm4MSA_fgft(self):
        print("Test pyfaust.fact._palm4msa_fgft()")
        from pyfaust.factparams import ConstraintReal,\
                ConstraintInt, ConstraintName, StoppingCriterion
        #from pyfaust.factparams import ParamsPalm4MSAFGFT, StoppingCriterion
        from numpy import diag,copy
        from pyfaust.fact import _palm4msa_fgft
        from pyfaust.factparams import ParamsPalm4MSAFGFT

        L = \
        loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/ref_test_PALM4SMA_FFT2")['data']
        init_D = \
        loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/ref_test_PALM4SMA_FFT2")['p_init_D']
        init_D = copy(diag(init_D))
        init_facts1 = \
                loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/ref_test_PALM4SMA_FFT2")['p_init_facts1']
        init_facts2 = \
                loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/ref_test_PALM4SMA_FFT2")['p_init_facts2']
        init_facts = [ init_facts1, init_facts2 ]
        L = L.astype(np.float64)
        # other params should be set from file also but do it manually
        cons1 = ConstraintInt(ConstraintName(ConstraintName.SP), 128, 128,
                              12288)
        cons2 = ConstraintInt(ConstraintName(ConstraintName.SP), 128,
                                 128, 384)
        stop_crit = StoppingCriterion(num_its=100)
        param = ParamsPalm4MSAFGFT([cons1, cons2], stop_crit,
                                   init_facts=init_facts,
                                   init_D=init_D,
                                   is_update_way_R2L=False, init_lambda=128,
                                   is_verbose=True, step_size=1e-6)
        F, D, _lambda = _palm4msa_fgft(L, param, ret_lambda=True)
        print("Lap norm:", norm(L, 'fro'))
        print("out lambda:", _lambda)
        D = diag(D)
        err = norm((F.toarray()@D)@F.T.toarray()-L,"fro")/norm(L,"fro")
        print("err: ", err)
        # the error reference is from the C++ test,
        # misc/test/src/C++/test_Palm4MSAFFT.cpp.in
        self.assertAlmostEqual(err, 1.39352e-5, places=5)

    def testFactHierarchFGFT(self):
        print("Test pyfaust.fact.fgft_palm()")
        import pyfaust.fact
        from pyfaust.factparams import ParamsHierarchical, StoppingCriterion
        from pyfaust.factparams import ConstraintReal, ConstraintInt,\
                ConstraintName
        from numpy import diag, copy
        num_facts = 4
        is_update_way_R2L = False
        init_lambda = 1.0
        #M = np.random.rand(500, 32)
        U = \
        loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/HierarchicalFactFFT_test_U_L_params.mat")['U']
        Lap = \
        loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/HierarchicalFactFFT_test_U_L_params.mat")['Lap'].astype(np.float)
        init_D = \
        loadmat(dirname(sys.argv[0])+"/../../../../misc/data/mat/HierarchicalFactFFT_test_U_L_params.mat")['init_D']
        params_struct = \
        loadmat(dirname(sys.argv[0])+'/../../../../misc/data/mat/HierarchicalFactFFT_test_U_L_params.mat')['params']
        nfacts = params_struct['nfacts'][0,0][0,0] #useless
        niter1 = params_struct['niter1'][0,0][0,0]
        niter2 = params_struct['niter2'][0,0][0,0]
        verbose = params_struct['verbose'][0,0][0,0]==1
        is_update_way_R2L = params_struct['update_way'][0,0][0,0]==1
        init_lambda = params_struct['init_lambda'][0,0][0,0]
        # default step_size
        stepsize = params_struct['stepsize'][0,0][0,0]
        factside = params_struct['fact_side'][0,0][0,0] == 1 # ignored
        # for convenience I set the constraints manually and don't take them
        # from mat file, but they are the same
        fact0_cons = ConstraintInt(ConstraintName(ConstraintName.SP), 128,
                                   128, 12288)
        fact1_cons = ConstraintInt(ConstraintName(ConstraintName.SP), 128, 128,
                                  6144)
        fact2_cons = ConstraintInt(ConstraintName(ConstraintName.SP), 128, 128,
                                  3072)
        res0_cons = ConstraintInt(ConstraintName(ConstraintName.SP), 128,
                                   128, 384)
        res1_cons =  ConstraintInt(ConstraintName(ConstraintName.SP), 128, 128,
                                  384)
        res2_cons =  ConstraintInt(ConstraintName(ConstraintName.SP), 128, 128,
                                   384)
        stop_crit1 = StoppingCriterion(num_its=niter1)
        stop_crit2 = StoppingCriterion(num_its=niter2)
        param = ParamsHierarchical([fact0_cons, fact1_cons, fact2_cons],
                                       [res0_cons, res1_cons, res2_cons],
                                       stop_crit1, stop_crit2,
                                       is_verbose=verbose,
                                       init_lambda=init_lambda,
                                       constant_step_size=False)
        diag_init_D = copy(diag(init_D))
        print("norm(init_D):", norm(init_D))
        F,D, _lambda = pyfaust.fact.fgft_palm(U, Lap, param,
                                                           diag_init_D,
                                                           ret_lambda=True)
        print("out_lambda:", _lambda)
        self.assertEqual(F.shape, U.shape)
        D = diag(D)
        err = norm((F.toarray()@D)@F.T.toarray()-Lap,"fro")/norm(Lap,"fro")
        # matrix to factorize and reference relative error come from
        # misc/test/src/C++/hierarchicalFactorizationFFT.cpp
        self.assertAlmostEqual(err, 0.08480, places=4)

    def test_splin(self):
        from pyfaust.proj import splin
        from random import randint
        from numpy.random import rand
        from numpy import count_nonzero
        from numpy.linalg import norm
        min_n, min_m = 5, 5
        m = randint(min_m, 128)
        n = randint(min_n, 128)
        M = rand(m,n)
        k = randint(1,n)
        p = splin((m,n),k, normalized=True)
        Mp = p(M)
        for i in range(0,m):
            # np.savez('M.npz', M, Mp, k)
            self.assertLessEqual(count_nonzero(Mp[i,:]), k)
        self.assertAlmostEqual(norm(Mp), 1)

    def test_spcol(self):
        from pyfaust.proj import spcol
        from random import randint
        from numpy.random import rand
        from numpy import count_nonzero
        from numpy.linalg import norm
        min_n, min_m = 5, 5
        m = randint(min_m, 128)
        n = randint(min_n, 128)
        M = rand(m,n)
        k = randint(1,m)
        p = spcol((m,n),k, normalized=True)
        Mp = p(M)
        for i in range(0,n):
            # np.savez('M.npz', M, Mp, k)
            self.assertLessEqual(count_nonzero(Mp[:,i]), k)
        self.assertAlmostEqual(norm(Mp), 1)

    def test_splincol(self):
        from pyfaust.proj import splincol
        from random import randint
        from numpy.random import rand
        from numpy import count_nonzero
        from numpy.linalg import norm
        min_n, min_m = 5, 5
        m = randint(min_m, 128)
        n = randint(min_n, 128)
        M = rand(m,n)
        k = randint(1,m)
        p = splincol((m,n),k, normalized=True)
        Mp = p(M)
        # TODO: define sparsity assertions to verify on Mp
        self.assertAlmostEqual(norm(Mp), 1)

    def test_sp(self):
        from pyfaust.proj import sp
        from random import randint
        from numpy.random import rand
        from numpy import count_nonzero
        from numpy.linalg import norm
        min_n, min_m = 5, 5
        m = randint(min_m, 128)
        n = randint(min_n, 128)
        M = rand(m,n)
        k = randint(1,m*n)
        p = sp((m,n),k, normalized=True)
        Mp = p(M)
        for i in range(0,n):
            # np.savez('M.npz', M, Mp, k)
            self.assertLessEqual(count_nonzero(Mp[:,i]), k)
        self.assertAlmostEqual(norm(Mp), 1)

    def test_supp(self):
        from pyfaust.proj import supp
        from numpy.random import rand
        from numpy.random import randn, permutation as randperm
        from numpy.linalg import norm
        from random import randint
        min_n, min_m = 5, 5
        m = randint(min_m, 128)
        n = randint(min_n, 128)
        M = rand(m,n)
        k = randint(1, min(m,n))
        nnz_rinds = randperm(m)[:k]
        nnz_cinds = randperm(n)[:k]
        S = np.zeros((m,n))
        S[nnz_rinds, nnz_cinds] = 1
        p = supp(S, normalized=True)
        pM = p(M)
        # same nnz number
        self.assertEqual(np.count_nonzero(pM), np.count_nonzero(S))
        # same support
        self.assertTrue(np.allclose(pM != 0, S != 0))
        # pM normalized (according to fro-norm)
        self.assertAlmostEqual(norm(pM), 1)
        # same projection without normalization
        p = supp(S, normalized=False)
        pM = p(M)
        self.assertTrue(np.allclose(pM[pM != 0], M[S != 0]))

    def test_const(self):
        from pyfaust.proj import const
        from numpy.random import rand
        from random import randint
        min_n, min_m = 5, 5
        m = randint(min_m, 128)
        n = randint(min_n, 128)
        M = rand(m,n)
        C = rand(m,n)
        p = const(C, normalized=False)
        pM = p(M)
        self.assertTrue(np.allclose(C, pM))

    def test_normcol(self):
        from pyfaust.proj import normcol
        from numpy.random import rand
        from numpy.random import randn, permutation as randperm
        from numpy.linalg import norm
        from random import randint, random
        min_n, min_m = 5, 5
        m = randint(min_m, 128)
        n = randint(min_n, 128)
        M = rand(m,n)*random()*50
        k = random()*50
        p = normcol(M.shape,k)
        pM = p(M)
        for i in range(pM.shape[1]):
            self.assertAlmostEqual(norm(pM[:,i]), k)

    def test_normlin(self):
        from pyfaust.proj import normlin
        from numpy.random import rand
        from numpy.random import randn, permutation as randperm
        from numpy.linalg import norm
        from random import randint, random
        min_n, min_m = 5, 5
        m = randint(min_m, 128)
        n = randint(min_n, 128)
        M = rand(m,n)*random()*50
        k = random()*50
        p = normlin(M.shape,k)
        pM = p(M)
        for i in range(pM.shape[0]):
            self.assertAlmostEqual(norm(pM[i,:]), k)


    def test_toeplitz(self):
        print("test_toeplitz")
        M = np.array([[0.0912995 , 0.94030087, 0.83110585, 0.53118248, 0.93799911],
                   [0.80432231, 0.63536838, 0.85907982, 0.13178579, 0.33610065],
                   [0.65570143, 0.80996454, 0.58217456, 0.90521081, 0.69633737],
                   [0.67242761, 0.04777302, 0.61627327, 0.47974862, 0.06107523],
                   [0.97367686, 0.27061922, 0.28952798, 0.22093562, 0.56537747]])
        pM_ref = np.array([[0.47079371, 0.69141668, 0.55307634, 0.43364156, 0.93799911],
                        [0.61287393, 0.47079371, 0.69141668, 0.55307634, 0.43364156],
                        [0.33100081, 0.61287393, 0.47079371, 0.69141668, 0.55307634],
                        [0.47152341, 0.33100081, 0.61287393, 0.47079371, 0.69141668],
                        [0.97367686, 0.47152341, 0.33100081, 0.61287393,
                         0.47079371]])
        from pyfaust.proj import toeplitz
        from numpy.random import rand
        p = toeplitz(M.shape, normalized=False)
        self.assertTrue(np.allclose(p(M), pM_ref))

    def test_hankel(self):
        print("test_hankel")
        M = np.array([[0.0912995 , 0.94030087, 0.83110585, 0.53118248, 0.93799911],
                   [0.80432231, 0.63536838, 0.85907982, 0.13178579, 0.33610065],
                   [0.65570143, 0.80996454, 0.58217456, 0.90521081, 0.69633737],
                   [0.67242761, 0.04777302, 0.61627327, 0.47974862, 0.06107523],
                   [0.97367686, 0.27061922, 0.28952798, 0.22093562, 0.56537747]])
        pM_ref = np.array([[0.0912995 , 0.87231159, 0.70739189, 0.71816361, 0.53468187],
                        [0.87231159, 0.70739189, 0.71816361, 0.53468187, 0.53205099],
                        [0.70739189, 0.71816361, 0.53468187, 0.53205099, 0.48853799],
                        [0.71816361, 0.53468187, 0.53205099, 0.48853799, 0.14100543],
                        [0.53468187, 0.53205099, 0.48853799, 0.14100543, 0.56537747]])
        from pyfaust.proj import hankel
        p = hankel(M.shape, normalized=False)
        self.assertTrue(np.allclose(p(M), pM_ref))

    def test_circ(self):
        print("test_circ")
        M = np.array([[0.0912995 , 0.94030087, 0.83110585, 0.53118248, 0.93799911],
                   [0.80432231, 0.63536838, 0.85907982, 0.13178579, 0.33610065],
                   [0.65570143, 0.80996454, 0.58217456, 0.90521081, 0.69633737],
                   [0.67242761, 0.04777302, 0.61627327, 0.47974862, 0.06107523],
                   [0.97367686, 0.27061922, 0.28952798, 0.22093562, 0.56537747]])
        pM_ref = np.array([[0.47079371, 0.74786872, 0.52045517, 0.37205711, 0.67789897],
                        [0.67789897, 0.47079371, 0.74786872, 0.52045517, 0.37205711],
                        [0.37205711, 0.67789897, 0.47079371, 0.74786872, 0.52045517],
                        [0.52045517, 0.37205711, 0.67789897, 0.47079371, 0.74786872],
                        [0.74786872, 0.52045517, 0.37205711, 0.67789897, 0.47079371]])
        from pyfaust.proj import circ
        from numpy.random import rand
        p = circ(M.shape, normalized=False)
        self.assertTrue(np.allclose(p(M), pM_ref))

    def test_blockdiag(self):
        print("test_blockdiag")
        from pyfaust.proj import blockdiag
        from numpy.random import rand
        from numpy.linalg import norm
        import numpy as np
        M = rand(15,15)
        shapes = [(1,1), (2,2), (3,3), (4,4), (5,5)]
        p = blockdiag(M.shape, shapes, normalized=False)
        pM = p(M)
        boundaries = [(1,1), (3,3), (6,6), (10,10), (15,15)]
        M_blocks = [ M[0:1, 0:1], M[1:3,1:3], M[3:6, 3:6], M[6:10, 6:10],
                  M[10:15, 10:15] ]
        pM_blocks = [ pM[0:1, 0:1], pM[1:3,1:3], pM[3:6, 3:6], pM[6:10, 6:10],
                  pM[10:15, 10:15] ]
        for i in range(len(M_blocks)):
            self.assertTrue(np.allclose(M_blocks[i], pM_blocks[i]))
        self.assertTrue(not np.allclose(M,pM))
        # TODO: verify the fro norm of pM is equal to the norm of the vector
        # composed of all entries of the blocks gathered from M directly
        v = np.zeros((1,)) 
        for i in range(len(M_blocks)):
           v = np.concatenate((v, M_blocks[i].flatten()))
        self.assertAlmostEqual(norm(v), norm(pM))

    def test_skperm(self):
        print("test_skperm")
        from pyfaust.proj import skperm
        M = np.array([[-0.04440802, -0.17569296, -0.02557815, -0.15559154],
                   [-0.0083095,  -3.38725936, -0.78484126, -0.4883618 ],
                   [-1.48942563, -1.71787215, -0.84000212, -3.71752454],
                   [-0.88957883, -0.19107863, -5.92900636, -6.51064175]])
        k = 2
        p = skperm(M.shape, k, normalized=False)
        ref_pM = [[-0.04440802,0.,-0.02557815,0.,],
                  [-0.0083095,-3.38725936,0.,0.,],
                  [ 0.,-1.71787215,0.,-3.71752454],
                  [ 0.,0.,-5.92900636,-6.51064175]]
        self.assertTrue(np.allclose(p(M), ref_pM))

    def test_proj_id(self):
        from pyfaust.proj import proj_id
        from numpy.random import rand
        M = rand(32, 33)
        p = proj_id(M.shape)
        self.assertTrue(np.allclose(p(M), M))

    def test_default_proj_cons_attributes(self):
        print("Test default normalized/pos attributes of"
              " pyfaust.proj*/pyfaust.factparams.Constraint*")
        from pyfaust.proj import (proj_id, toeplitz, hankel, circ, blockdiag,
                                  supp, const, sp, splin, spcol, splincol,
                                  skperm, normcol, normlin)
        from pyfaust.factparams import ConstraintMat, ConstraintInt, ConstraintReal
        shape = (10,11)
        # Matrix Constraint/projectors
        c = ConstraintMat('id', shape=shape)
        p = proj_id(shape)
        self.assertTrue(c.normalized == False  == p.constraint.normalized and
                        c.pos == False == p.constraint.pos)
        c = ConstraintMat('toeplitz', shape=shape)
        p = toeplitz(shape)
        self.assertTrue(c.normalized == p.constraint.normalized == True and
                        c.pos == p.constraint.pos == False)
        c = ConstraintMat('circ', shape=shape)
        p = circ(shape)
        self.assertTrue(c.normalized == p.constraint.normalized == True and
                        c.pos == p.constraint.pos == False)
        c = ConstraintMat('hankel', shape=shape)
        p = hankel(shape)
        self.assertTrue(c.normalized == p.constraint.normalized == True and
                        c.pos == p.constraint.pos == False)
        S = np.zeros(shape)
        S[np.random.rand(*shape) > .5] = 1
        c = ConstraintMat('supp', S)
        p = supp(S)
        self.assertTrue(c.normalized == p.constraint.normalized == True and
                        c.pos == p.constraint.pos == False)
        C = np.random.rand(*shape)
        c = ConstraintMat('const', C)
        p = const(C)
        self.assertTrue(c.normalized == p.constraint.normalized == False and
                        c.pos == p.constraint.pos == False)
        p = blockdiag(C.shape, [(1,1), (2,2), (3,3), (4,5)])
        self.assertTrue(p.constraint.normalized == True and
                        p.constraint.pos == False)
        # Int Constraint/projectors
        c = ConstraintInt('sp', shape[0], shape[1], 15)
        p = sp(shape, 15)
        self.assertTrue(c.normalized == p.constraint.normalized == True and
                        c.pos == p.constraint.pos == False)
        c = ConstraintInt('splin', shape[0], shape[1], 2)
        p = splin(shape, 2)
        self.assertTrue(c.normalized == p.constraint.normalized == True and
                        c.pos == p.constraint.pos == False)
        c = ConstraintInt('spcol', shape[0], shape[1], 2)
        p = spcol(shape, 2)
        self.assertTrue(c.normalized == p.constraint.normalized == True and
                        c.pos == p.constraint.pos == False)
        c = ConstraintInt('splincol', shape[0], shape[1], 2)
        p = splincol(shape, 2)
        self.assertTrue(c.normalized == p.constraint.normalized == True and
                        c.pos == p.constraint.pos == False)
        c = ConstraintInt('skperm', shape[0], shape[1], 2)
        p = skperm(shape, 2)
        self.assertTrue(c.normalized == p.constraint.normalized == True and
                        c.pos == p.constraint.pos == False)
        # Real Constraint/projectors
        c = ConstraintReal('normlin', shape[0], shape[1], .6)
        p = normlin(shape, .6)
        self.assertTrue(c.normalized == p.constraint.normalized == False and
                        c.pos == p.constraint.pos == False)
        c = ConstraintReal('normcol', shape[0], shape[1], .6)
        p = normcol(shape, .6)
        self.assertTrue(c.normalized == p.constraint.normalized == False and
                        c.pos == p.constraint.pos == False)

    def test_butterfly(self):
        print("test pyfaust.fact.butterfly")
        from pyfaust import wht, dft
        from pyfaust.fact import butterfly
        H = wht(64).toarray()
        for dir in ['right', 'left', 'bbtree']:
            FH = butterfly(H, type=dir)
            self.assertAlmostEqual((FH-H).norm()/norm(H), 0)
        D = dft(64).toarray()
        for dir in ['right', 'left', 'bbtree']:
            FD = butterfly(D, type=dir, perm='bitrev')
            self.assertAlmostEqual((FD-D).norm()/norm(D), 0)

    def test_palm4msa_mhtp(self):
        print("Test palm4msa_mhtp")
        # call Palm4MSA specifying params
        from os import dup2, pipe # for
        from pyfaust.fact import palm4msa
        from pyfaust.factparams import (ParamsPalm4MSA, ConstraintList,
                                        StoppingCriterion,
                                        ConstraintInt,
                                        ConstraintReal, ParamsFact)
        import numpy as np
        from tempfile import gettempdir
        from os.path import join
        from pyfaust.proj import splin, normcol
        from pyfaust.factparams import MHTPParams
        from pyfaust.fact import palm4msa_mhtp
        M = np.random.rand(500, 32)
        stop_crit = StoppingCriterion(num_its=200)
        cons = [ splin((500,32), 5), normcol((32,32), 1.0)]
        param = ParamsPalm4MSA(cons, stop_crit)
        param.is_verbose = True
        # MHTP will run every 100 iterations of PALM4MSA (that is 2 times) for 50
        # iterations on each factor
        mhtp_param = MHTPParams(num_its=50, palm4msa_period=100)
        tmp_dir = gettempdir()
        tmp_file = join(tmp_dir, "verbose_output_of_palm4msa_mhtp_test")
        print("tmp_file:", tmp_file)
        f = open(tmp_file, 'w')
        dup2(1,2)
        dup2(f.fileno(), 1)
        F = palm4msa_mhtp(M, param, mhtp_param)
        print()
        f.close()
        dup2(2,1)
        factor0_line = factor1_line = False
        with open(tmp_file, 'r') as lines:
            for line in lines.readlines():
                print(line, end='')
                factor0_line = factor0_line or line.startswith('Starting a MHTP pass (50 iterations) for factor #0')
                factor1_line = factor1_line or line.startswith('Starting a MHTP pass (50 iterations) for factor #1')
            self.assertTrue(factor0_line)
            self.assertTrue(factor1_line)
            self.assertLessEqual(norm((F.toarray()-M), "fro")/norm(M,"fro"), 0.4)
        os.remove(tmp_file)


    def test_hierarchical_mhtp(self):
       print("Test hierarchical_mhtp")
       from os import dup2, pipe # for
       from tempfile import gettempdir
       from os.path import join 
       from pyfaust.fact import hierarchical_mhtp
       from pyfaust.factparams import ParamsHierarchical, StoppingCriterion
       from pyfaust.factparams import MHTPParams
       from pyfaust.proj import sp, normcol, splin
       import numpy as np
       M = np.random.rand(500, 32)
       fact_cons = [splin((500, 32), 5), sp((32,32), 96), sp((32,32), 96)]
       res_cons = [normcol((32,32), 1), sp((32,32), 666), sp((32,32), 333)]
       stop_crit1 = StoppingCriterion(num_its=200)
       stop_crit2 = StoppingCriterion(num_its=200)
       # 50 iterations of MHTP will run every 100 iterations of PALM4MSA (each time PALM4MSA is called by the hierarchical algorithm)
       mhtp_param = MHTPParams(num_its=50, palm4msa_period=100)
       param = ParamsHierarchical(fact_cons, res_cons, stop_crit1, stop_crit2)
       param.is_verbose = True
       F = hierarchical_mhtp(M, param, mhtp_param)
       self.assertLessEqual(norm((F.toarray()-M), "fro")/norm(M,"fro"), 0.5)

    def test_hierarchical_dft(self):
        print("Test hierarchical dft")
        from pyfaust import dft
        from pyfaust.fact import hierarchical
        DFT = dft(32).toarray()
        F = hierarchical(DFT, 'dft', backend=2020)
        err = norm(F.toarray()-DFT)/norm(DFT)
        self.assertLessEqual(err, 1e-6)
        F = hierarchical(DFT, 'dft', backend=2016)
        err = norm(F.toarray()-DFT)/norm(DFT)
        self.assertLessEqual(err, 1e-6)

    def test_single_bsr_Faust(self):
        from scipy.sparse import bsr_matrix
        from pyfaust import Faust
        from numpy import allclose
        from numpy.random import rand
        nzblocks = rand(3, 2, 3) # 3 blocks of size 2x3
        # create a scipy BSR matrix
        B = bsr_matrix((nzblocks, [0, 1, 2], [0, 1, 2, 3, 3, 3]), shape=(10,9))
        # create the single factor Faust
        F = Faust(B)
        self.assertTrue(allclose(F.toarray(), B.toarray()))

    def test_save_bsr_Faust(self):
        from scipy.sparse import bsr_matrix, random
        from pyfaust import Faust
        from numpy import allclose
        from numpy.random import rand
        nzblocks = rand(3, 2, 3) # 3 blocks of size 2x3
        # create a scipy BSR matrix
        B = bsr_matrix((nzblocks, [0, 1, 2], [0, 1, 2, 3, 3, 3]), shape=(10,9))
        # create the single factor Faust
        F = Faust([B, B.T, rand(10, 18), random(18, 18, .2)])
        F.save('test_bsr_faust.mat')
        G = Faust('test_bsr_faust.mat')
        self.assertTrue(allclose(F.toarray(), G.toarray()))
        G = Faust.load('test_bsr_faust.mat')
        self.assertTrue(allclose(F.toarray(), G.toarray()))

    def test_bsr_Faust(self):
        from random import randint
        from scipy.sparse import random
        from numpy.random import rand
        def rand_bsr(m, n, bm, bn, bnnz):
            nblocks_per_col = m//bm
            nblocks_per_row = n//bn
            nblocks = nblocks_per_col*nblocks_per_row
            # 1D possible nz block indices
            BI = list(range(nblocks))
            # choose bnnz ones among them
            nzBI = []
            for i in range(bnnz):
                ri = randint(0,len(BI)-1)
                nzBI.append(BI[ri])
                del BI[ri]
            nzBI.sort()
            indices = np.array(nzBI)%nblocks_per_row
            bcolinds_ind = 0
            bindptr_ind = 1
            indptr = np.zeros((int(nblocks_per_col+1)))
            for bi in nzBI:
                while bi//nblocks_per_row+1 > bindptr_ind:
                    bindptr_ind += 1
                indices[bcolinds_ind] = bi%nblocks_per_row
                bcolinds_ind += 1
                indptr[bindptr_ind] += 1
            for i in range(1, int(nblocks_per_col)+1):
                indptr[i] += indptr[i-1]
            data = rand(bnnz, bm, bn)
            return data, indices, indptr

        F_factors = []
        F_sp_factors = []
        for i in range(10):
            nrows = ncols = 1024
            bnnz = 20
            bm = bn = 64
            data, indices, indptr = rand_bsr(nrows,ncols, bm, bn, bnnz)
            F_factors += [bsr_matrix((data, indices, indptr), shape=(nrows, ncols))]
            F_sp_factors.append(F_factors[i].tocsr())
        bsrF = Faust(F_factors)
        spF = Faust(F_sp_factors)
        print("toarray err:", np.linalg.norm(bsrF.toarray()-spF.toarray())/np.linalg.norm(spF.toarray()))
        self.assertTrue(np.allclose(bsrF.toarray(), spF.toarray()))
        M = rand(bsrF.shape[1], bsrF.shape[0])
        print("err ds:", np.linalg.norm(bsrF@M-spF@M)/np.linalg.norm(spF@M))
        self.assertTrue(np.allclose(bsrF@M, spF@M))
        M = random(bsrF.shape[1], bsrF.shape[0], .02, format='csr')
        print("err sp:", np.linalg.norm(bsrF@M-spF@M)/np.linalg.norm(spF@M))
        self.assertTrue(np.allclose(bsrF@M, spF@M))


if __name__ == "__main__":
    if(len(sys.argv)> 1):
        # argv[1] is for adding a directory in PYTHONPATH
        # (to find pyfaust module)
        #sys.path.append(sys.argv[1])
        sys.path = [sys.argv[1]]+sys.path # gives priority to pyfaust path
        del sys.argv[1] # deleted to avoid interfering with unittest
    from pyfaust import Faust
    if(len(sys.argv) > 1):
        #ENOTE: test only a single test if name passed on command line
        try:
            singleton = unittest.TestSuite()
            singleton.addTest(TestFaustPy(sys.argv[1]))
            unittest.TextTestRunner().run(singleton)
        except:
            # if TestFaustPy failed to launch the unit test
            # maybe it exists in TestFaustFactory
            singleton = unittest.TestSuite()
            singleton.addTest(TestFaustFactory(sys.argv[1]))
            unittest.TextTestRunner().run(singleton)
    else:
        unittest.main()
