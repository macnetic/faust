from __future__ import print_function
from pylab import *
import os, sys

rcgs = [2, 4, 8]
dims = [128, 256, 512]
nb_facts = [2, 4, 8]
nb_mults = 500
constraint = True # 'per_row' # per_col
dims_len = len(dims)
rcgs_len = len(rcgs)
nb_facts_len = len(nb_facts)

def runtime_comparison():
    """
    TODO
    """
    from pyfaust import Faust, FaustFactory

    matrix_or_vector = 'vector'


    fausts = ndarray(shape=(dims_len, rcgs_len, nb_facts_len), dtype=Faust)
    dense_mats =  list()


    # loading all the different Fausts and dense matrices
    for j in range(0,dims_len):
        dim = dims[j]
        A = rand(dim, dim)
        dense_mats.append(A)

        for k in range(0,rcgs_len):
            rcg = rcgs[k]

            for l in range(0,nb_facts_len):
                nf = nb_facts[l]
                F = FaustFactory.rand(nf , dim, 1./(nf*rcg),
                                      per_row=constraint, fac_type='sparse')
                assert(F.rcg() == rcg)
                fausts[j,k,l] = F
                assert(F.shape == (dim,dim))

    # time comparison
    from time import time, clock
    if sys.platform == 'win32':
        timer = clock
    else:
        timer = time

    tdense = ndarray(shape=(nb_mults, dims_len, 2))
    tfaust = ndarray(shape=(nb_mults, dims_len, rcgs_len, nb_facts_len, 2))

    for i in range(0, nb_mults):
        print("\r\r #muls =",i+1,'/',nb_mults, end='')
        for j in range(0,dims_len):
            dim = dims[j]

            if(matrix_or_vector == 'matrix'):
                dim2 = dim # mul. by a square matrix
            elif(matrix_or_vector == 'vector'):
                dim2 = 1
            else:
                raise("matrix_or_vector string must be equal to matrix or"
                      " vector")

            for k in range(0,rcgs_len):
                rcg = rcgs[k]

                for l in range(0,nb_facts_len):
                    nf = nb_facts[l]
                    x = rand(dim,dim2)

                    if(k == 0 and l == 0):
                        A = dense_mats[j]
                        t = timer()
                        y = A.dot(x)
                        tdense[i,j,0] = timer()-t
                        t = timer()
                        y_trans = A.T.dot(x)
                        tdense[i,j,1] = timer()-t

                    F = fausts[j,k,l]
                    t = timer()
                    yfaust = F*x
                    tfaust[i,j,k,l,0] = timer()-t
                    t = timer()
                    yfaust_trans = F.T*x
                    tfaust[i,j,k,l,1] = timer()-t

    print()
    output_dir = 'pyfaust_demo_output'
    import os.path
    if(not os.path.exists(output_dir)):
       os.mkdir(output_dir)
    path_tfaust = output_dir+os.sep+'runtime_cmp_tfaust-'+matrix_or_vector+'.txt'
    path_tdense = output_dir+os.sep+'runtime_cmp_tdense-'+matrix_or_vector+'.txt'
    savetxt(path_tdense, tdense.reshape(tdense.size))
    savetxt(path_tfaust, tfaust.reshape(tfaust.size))

    # test
    #tfaust_r = loadtxt(path_tfaust).reshape(tfaust.shape)
    #assert(all(tfaust_r == tfaust))
    #tdense_r = loadtxt(path_tdense).reshape(tdense.shape)
    #assert(all(tdense_r == tdense))



