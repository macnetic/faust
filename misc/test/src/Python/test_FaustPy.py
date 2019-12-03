import unittest
from scipy import sparse
import random
import tempfile
import os
import sys
import numpy as np
from scipy.io import savemat,loadmat
from numpy.linalg import norm
from scipy.sparse import spdiags
import math

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
            ref_full_NF = F.todense()
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
                       msg="\nref_full_F=\n"+str(F.todense())+"\nvec_ref_full_NF=\n"+str(vec_ref_full_NF)+"\nvec_test_NF=\n"+ \
                       str(vec_test_NF.toarray())+"\nF=\n"+ \
                       str(str(F.save('/tmp/normalize_test.mat'))+", i="+str(i)))


    def testNormInf(self):
        print("testNormInf()")
        ref_norm = norm(self.mulFactors(), np.inf)
        test_norm = self.F.norm(np.inf)
        self.F.save("norminf_test.mat");
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
        test_prod = self.F[::,::].todense()
        self.assertProdEq(prod, test_prod)
        # reminder the Faust is of minimal size 3 for the both dims (MIN_DIM_SIZE)
        # test_prod = self.F[...,...].todense() # forbidden (only one index can
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
        row = self.F[rand_i,...].todense()
        for j in range(0,self.F.shape[1]):
            if(row[0,j] == 0):
                self.assertEqual(prod[rand_i,j], 0)
            else:
                self.assertLessEqual(abs(row[0,j]-(prod[rand_i,j]))/prod[rand_i,j],10**-6)
        # test one random col
        rand_j = self.r.randint(0,self.F.shape[1]-2)
        col = self.F[..., rand_j].todense()
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
        n1 = norm(sF.toarray()-self.F.todense()[rand_i:rand_ii,rand_j:rand_jj])
        n2 = norm(self.F.todense()[rand_i:rand_ii,rand_j:rand_jj])
        if(n2 == 0): # avoid nan error
            self.assertLessEqual(n1,0.0005)
        else:
            self.assertLessEqual(n1/n2,0.005)
        # test a second slice on sF only if the shape allows it # min shape (1,1)
        if(sF.shape[0] >= 2 and sF.shape[1] >= 2):
            rand_i, rand_j = \
            self.r.randint(0,sF.shape[0]-2),self.r.randint(0,sF.shape[1]-2)
            rand_ii, rand_jj = \
                    self.r.randint(rand_i+1,sF.shape[0]-1),self.r.randint(rand_j+1,sF.shape[1]-1)
            sF2 = sF[rand_i:rand_ii,rand_j:rand_jj]
            n1 = norm(sF2.todense()-sF.todense()[rand_i:rand_ii,rand_j:rand_jj])
            n2 = norm(sF.todense()[rand_i:rand_ii,rand_j:rand_jj])
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
        col_ids = [ self.r.randint(0,F.shape[1]-1) for i in
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
        [ frand(5,F.shape[0])*Faust(rand(F.shape[0],F.shape[1])),
         frand(5,F.shape[0], .5,
                           field='complex')*Faust(rand(F.shape[0],F.shape[1]))]
        for i in range(0,len(fausts)):
            F2 = fausts[i]
            self.assertAlmostEqual(norm((F+F2).toarray()-
                                        (F.toarray()+F2.toarray())), 0,
                                   places=2)

    def testMinus(self):
        from pyfaust import rand as frand
        from numpy.random import rand
        print("testMinus()")
        print("substraction of a Faust-scalar")
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
        print("substraction Faust-Faust")
        fausts = \
        [ frand(5,F.shape[0])*Faust(rand(F.shape[0],F.shape[1])),
         frand(5,F.shape[0], .5,
                           field='complex')*Faust(rand(F.shape[0],F.shape[1]))]
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
        r_fausts = [ frand(self.r.randint(1,10), F.shape[1]),
                    frand(self.r.randint(1,10), F.shape[1], field='complex')]
        for i in range(0,len(r_fausts)):
            rF = r_fausts[i]
            assert(isinstance(rF, Faust))
            test_prod = F*rF
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
        test_prod = F*D
        self.assertTrue(np.allclose(test_prod, F.toarray().dot(D.toarray())))
        print("test mul of a Faust by a complex dia_matrix")
        D = \
        dia_matrix((np.random.rand(1,F.shape[1])+np.random.rand(1,F.shape[1])*np.complex(0,1),np.array([0])),
                       shape=(F.shape[1],F.shape[1]))
        test_prod = F*D
        self.assertTrue(np.allclose(test_prod, F.toarray().dot(D.toarray())))
        Mr = csr_matrix(rand(F.shape[1],10))
        Mc = csr_matrix(rand(F.shape[1],10)+np.complex(0,1)*rand(F.shape[1],10))
        print("test mul Faust-csr_matrix")
        self.assertTrue(np.allclose(F*Mr, F.toarray().dot(Mr.toarray())))
        print("test mul Faust-complex csr_matrix")
        self.assertTrue(np.allclose(F*Mc, F.toarray().dot(Mc.toarray())))

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
                G = rand(self.r.randint(1,TestFaustPy.MAX_NUM_FACTORS),
                                  F.shape[(cat_axis+1)%2],
                                  field=field)
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
                loadmat(sys.path[0]+"/../../../misc/data/mat/config_compared_palm2.mat")['data']
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
        E = F.todense()-M
        #print("err.:",norm(F.todense(), "fro"),  norm(E,"fro"), norm (M,"fro"))
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
        loadmat(sys.path[0]+"/../../../misc/data/mat/matrix_hierarchical_fact.mat")['matrix']
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
        E = F.todense()-M
        #print("err.:",norm(F.todense(), "fro"),  norm(E,"fro"), norm (M,"fro"),
        print("err: ", norm(E,"fro")/norm(M,"fro"))
        # matrix to factorize and reference relative error come from
        # misc/test/src/C++/hierarchicalFactorization.cpp
        self.assertAlmostEqual(norm(E,"fro")/norm(M,"fro"), 0.268513, places=5)

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
        loadmat(sys.path[0]+"/../../../misc/data/mat/matrix_hierarchical_fact.mat")['matrix']
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
        #print(F.todense())
        E = F.todense()-M
        #print("err.:",norm(F.todense(), "fro"),  norm(E,"fro"), norm (M,"fro"),
        print("err: ", norm(E,"fro")/norm(M,"fro"))
        # matrix to factorize and reference relative error come from
        # misc/test/src/C++/hierarchicalFactorization.cpp
        self.assertAlmostEqual(norm(E,"fro")/norm(M,"fro"), 1.0063, places=4)

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
        loadmat(sys.path[0]+"/../../../misc/data/mat/config_compared_palm2.mat")['data']
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
        #print(F.todense())
        E = F.todense()-M
        #print("err.:",norm(F.todense(), "fro"),  norm(E,"fro"), norm (M,"fro"))
        print("err:", norm(E,"fro")/norm(M,"fro"))
        print("lambda:", _lambda)
        # matrix to factorize and reference relative error come from
        # misc/test/src/C++/test_palm4MSA.cpp
        self.assertAlmostEqual(norm(E,"fro")/norm(M,"fro"), 0.29177, places=4)

    def testHadamard(self):
        print("Test pyfaust.wht()")
        from pyfaust import wht
        pow2_exp = random.Random().randint(1,10)
        n = 2**pow2_exp
        H = wht(n, False)
        fH = H.todense()
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
        fF = F.todense()
        ref_fft = fft(np.eye(n))
        self.assertAlmostEqual(norm(ref_fft-fF)/norm(ref_fft),0)
        assert(np.allclose(dft(n).toarray(),
                           dft(n, False).normalize().toarray()))

    def testFGFTGivens(self):
        print("Test fact.fgft_givens()")
        print(sys.path)
        import pyfaust.fact
        L = loadmat(sys.path[0]+"/../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat")['Lap']
        L = L.astype(np.float64)
        J = \
        int(loadmat(sys.path[0]+"/../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat")['J'])
        D, F = pyfaust.fact.fgft_givens(L, J, nGivens_per_fac=1, verbosity=0)
        D = spdiags(D, [0], L.shape[0], L.shape[0])
        print("Lap norm:", norm(L, 'fro'))
        err = norm((F*D.todense())*F.T.todense()-L,"fro")/norm(L,"fro")
        print("err: ", err)
        # the error reference is from the C++ test,
        # misc/test/src/C++/GivensFGFT.cpp.in
        self.assertAlmostEqual(err, 0.0326529, places=7)

    def testFGFTGivensParallel(self):
        print("Test pyfaust.fact.fgft_givens() -- parallel")
        from pyfaust.fact import eigtj, fgft_givens
        L = loadmat(sys.path[0]+"/../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat")['Lap']
        L = L.astype(np.float64)
        J = \
                int(loadmat(sys.path[0]+"/../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat")['J'])
        t = int(L.shape[0]/2)
        D, F = fgft_givens(L, J, nGivens_per_fac=t, verbosity=0)
        D = spdiags(D, [0], L.shape[0], L.shape[0])
        print("Lap norm:", norm(L, 'fro'))
        err = norm((F*D.todense())*F.T.todense()-L,"fro")/norm(L,"fro")
        print("err: ", err)
        # the error reference is from the C++ test,
        # misc/test/src/C++/GivensFGFTParallel.cpp.in (_double version)
        self.assertAlmostEqual(err, 0.0398154, places=7)
        D2, F2 = eigtj(L, J, nGivens_per_fac=t, verbosity=0)
        D2 = spdiags(D, [0], L.shape[0], L.shape[0])
        print("Lap norm:", norm(L, 'fro'))
        err2 = norm((F2*D.todense())*F2.T.todense()-L,"fro")/norm(L,"fro")
        print("err2: ", err2)
        # the error reference is from the C++ test,
        # misc/test/src/C++/GivensFGFTParallel.cpp.in
        self.assertEqual(err, err2)

    def testFGFTGivensSparse(self):
        print("Test fact.fgft_givens_sparse()")
        print(sys.path)
        from scipy.sparse import csr_matrix
        import pyfaust.fact
        L = loadmat(sys.path[0]+"/../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat")['Lap']
        L = L.astype(np.float64)
        J = \
        int(loadmat(sys.path[0]+"/../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat")['J'])
        D, F = pyfaust.fact.fgft_givens(csr_matrix(L), J, nGivens_per_fac=0, verbosity=0)
        D = spdiags(D, [0], L.shape[0], L.shape[0])
        print("Lap norm:", norm(L, 'fro'))
        err = norm((F*D.todense())*F.T.todense()-L,"fro")/norm(L,"fro")
        print("err: ", err)
        # the error reference is from the C++ test,
        # misc/test/src/C++/GivensFGFT.cpp.in
        self.assertAlmostEqual(err, 0.0326529, places=7)

    def testFGFTGivensParallelSparse(self):
        print("Test pyfaust.fact.fgft_givens_sparse() -- parallel")
        from pyfaust.fact import eigtj, fgft_givens
        from scipy.sparse import csr_matrix
        L = loadmat(sys.path[0]+"/../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat")['Lap']
        L = L.astype(np.float64)
        J = \
                int(loadmat(sys.path[0]+"/../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat")['J'])
        t = int(L.shape[0]/2)
        D, F = fgft_givens(csr_matrix(L), J, nGivens_per_fac=t, verbosity=0)
        D = spdiags(D, [0], L.shape[0], L.shape[0])
        print("Lap norm:", norm(L, 'fro'))
        err = norm((F*D.todense())*F.T.todense()-L,"fro")/norm(L,"fro")
        print("err: ", err)
        # the error reference is from the C++ test,
        # misc/test/src/C++/GivensFGFTParallel.cpp.in (_double version)
        self.assertAlmostEqual(err, 0.0398154, places=7)
        D2, F2 = eigtj(L, J, nGivens_per_fac=t, verbosity=0)
        D2 = spdiags(D, [0], L.shape[0], L.shape[0])
        print("Lap norm:", norm(L, 'fro'))
        err2 = norm((F2*D.todense())*F2.T.todense()-L,"fro")/norm(L,"fro")
        print("err2: ", err2)
        # the error reference is from the C++ test,
        # misc/test/src/C++/GivensFGFTParallel.cpp.in
        self.assertEqual(err, err2)

    def testeigtj_abserr(self):
        from numpy.random import rand
        from pyfaust.fact import eigtj
        from scipy.io import loadmat
        import sys
        err = .01
        M = rand(128,128)
        M = M.dot(M.T)
        D,U = eigtj(M, tol=err, relerr=False)
        self.assertAlmostEqual(norm(M-U*np.diag(D)*U.T), err, places=3 )
        L = loadmat(sys.path[0]+"/../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat")['Lap']
        L = L.astype(np.float64)
        M_cplx = L*np.complex(1,0) + L*np.complex(0,1)
        M_cplx = M_cplx.dot(np.matrix(M_cplx).H)
        err=.1
        D,U = eigtj(M_cplx, tol=err, relerr=False, verbosity=0)
        self.assertLessEqual(norm(M_cplx-U*np.diag(D)*U.H), err)

    def testeigtj_relerr(self):
        from numpy.random import rand
        from pyfaust.fact import eigtj
        err = .01
        M = rand(128,128)
        M = M.dot(M.T)
        D,U = eigtj(M, tol=err)
        self.assertLessEqual(norm(M-U*np.diag(D)*U.T)/norm(M), err)
        M_cplx = M + rand(128,128)*np.complex(0,1)
        M_cplx = M_cplx.dot(np.matrix(M_cplx).H)
        D,U = eigtj(M_cplx, tol=err)
        self.assertLessEqual(norm(M_cplx-U*np.diag(D)*U.H)/norm(M_cplx), err)

    def testFactPalm4MSA_fgft(self):
        print("Test pyfaust.fact._palm4msa_fgft()")
        from pyfaust.factparams import ConstraintReal,\
                ConstraintInt, ConstraintName, StoppingCriterion
        #from pyfaust.factparams import ParamsPalm4MSAFGFT, StoppingCriterion
        from numpy import diag,copy
        from pyfaust.fact import _palm4msa_fgft
        from pyfaust.factparams import ParamsPalm4MSAFGFT

        L = \
        loadmat(sys.path[0]+"/../../../misc/data/mat/ref_test_PALM4SMA_FFT2")['data']
        init_D = \
        loadmat(sys.path[0]+"/../../../misc/data/mat/ref_test_PALM4SMA_FFT2")['p_init_D']
        init_D = copy(diag(init_D))
        init_facts1 = \
                loadmat(sys.path[0]+"/../../../misc/data/mat/ref_test_PALM4SMA_FFT2")['p_init_facts1']
        init_facts2 = \
                loadmat(sys.path[0]+"/../../../misc/data/mat/ref_test_PALM4SMA_FFT2")['p_init_facts2']
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
        err = norm((F.todense()*D)*F.T.todense()-L,"fro")/norm(L,"fro")
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
        loadmat(sys.path[0]+"/../../../misc/data/mat/HierarchicalFactFFT_test_U_L_params.mat")['U']
        Lap = \
        loadmat(sys.path[0]+"/../../../misc/data/mat/HierarchicalFactFFT_test_U_L_params.mat")['Lap'].astype(np.float)
        init_D = \
        loadmat(sys.path[0]+"/../../../misc/data/mat/HierarchicalFactFFT_test_U_L_params.mat")['init_D']
        params_struct = \
        loadmat(sys.path[0]+'/../../../misc/data/mat/HierarchicalFactFFT_test_U_L_params.mat')['params']
        nfacts = params_struct['nfacts'][0,0][0,0] #useless
        niter1 = params_struct['niter1'][0,0][0,0]
        niter2 = params_struct['niter2'][0,0][0,0]
        verbose = params_struct['verbose'][0,0][0,0]==1
        is_update_way_R2L = params_struct['update_way'][0,0][0,0]==1
        init_lambda = params_struct['init_lambda'][0,0][0,0]
        stepsize = params_struct['stepsize'][0,0][0,0]
        factside = params_struct['fact_side'][0,0][0,0] == 1 # ignored
        # for convenience I set the constraints manually and don't take them
        # from mat file, but they are the same
        # default step_size
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
        err = norm((F.todense()*D)*F.T.todense()-Lap,"fro")/norm(Lap,"fro")
        # matrix to factorize and reference relative error come from
        # misc/test/src/C++/hierarchicalFactorizationFFT.cpp
        self.assertAlmostEqual(err, 0.084417, places=5)


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
        singleton = unittest.TestSuite()
        singleton.addTest(TestFaustPy(sys.argv[1]))
        unittest.TextTestRunner().run(singleton)
    else:
        unittest.main()
