import unittest
from scipy import sparse
import random
import tempfile
import os
import sys
import numpy as np
from scipy.io import savemat  # , loadmat


class TestFaustPy(unittest.TestCase):

    MAX_NUM_FACTORS = 64  # for the tested Faust

    def setUp(self):
        """ Initializes the tests objects """
        r = random.Random()  # initialized from time or system
        num_factors = r.randint(1, TestFaustPy.MAX_NUM_FACTORS)
        factors = []
        d2 = r.randint(1, 1000)
        for i in range(0, num_factors):
            d1, d2 = d2, r.randint(1, 1000)
            factors += [sparse.random(d1, d2, density=0.1, format='csr',
                        dtype=np.float64).todense()]
        self.F = Faust(factors)
        print("Tests on random Faust with dims=", self.F.get_nb_rows(),
              self.F.get_nb_cols())

    def testSave(self):
        tmp_dir = tempfile.gettempdir()+os.sep
        # save the Faust through Faust core API
        test_file = tmp_dir+"A.mat"
        self.F.save(test_file, format="Matlab")
        # save the Faust relying on numpy API
        ref_file = tmp_dir+"A_ref.mat"
        mdict = {'faust_factors':
                 np.ndarray(shape=(1, self.F.get_nb_factors()), dtype=object)}
        # self.F.display()
        for i in range(0, self.F.get_nb_factors()):
            mdict['faust_factors'][0, i] = self.F.get_factor(i)
        savemat(ref_file, mdict)
        # open the two saved files and compare the fausts
        F_test = Faust(test_file)
        F_ref = Faust(ref_file)
        # print("self.F.get_nb_factors()=", self.F.get_nb_factors())
        self.assertEqual(F_test.get_nb_factors(), F_ref.get_nb_factors())
        for i in range(0, F_ref.get_nb_factors()):
            fact_ref = F_ref.get_factor(i)
            fact_test = F_test.get_factor(i)
            self.assertEqual(fact_ref.shape, fact_test.shape)
            self.assertTrue((fact_ref == fact_test).all())


if __name__ == "__main__":
    if(len(sys.argv)> 1):
        # argv[1] is for adding a directory in PYTHONPATH
        # (to find FaustPy module)
        sys.path.append(sys.argv[1])
        del sys.argv[1] # deleted to avoid interfering with unittest
    from FaustPy import Faust
    unittest.main()
