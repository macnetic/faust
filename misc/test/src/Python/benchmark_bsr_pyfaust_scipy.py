import numpy as np
from numpy.random import rand
from random import randint
from scipy.sparse import bsr_matrix, random
from pyfaust import Faust
from time import process_time

"""
This scripts purpose is to benchmark BSR Fausts by comparing its times of
calculation to those obtained with CSR Fausts and scipy BSR Faust (as a Python
list). The operations used in the benchmark are all products: Faust-array,
Faust-CSR matrix and Faust-vector products.
"""

names = ['BSR Faust', 'CSR Faust', 'Array Faust', 'scipy BSR']


def mul_scipy_faust(SF, M):
    """
    Args:
        SF: a list of scipy matrices (aka a scipy Faust).
        M: the matrix to multiply by.
    """
    R = SF[-1]@M
    for i in range(len(SF)-2, -1, -1):
        R = SF[i]@R
    return R


def rand_bsr(m, n, bm, bn, bnnz):
    """
    Generates a random BSR matrix.

    Note: The following assertion must be True: m%bm == n%bn == 0.

    Args:
        m: the matrix number of rows.
        n: the matrix number of columns.
        bm: the matrix nonzero blocks number of rows.
        bn: the matrix nonzero blocks number of columns.
        bnnz: the number of nonzero blocks.
    """
    if not (m % bm == 0 and n % bn == 0):
        raise ValueError("bm and bn must evenly divide m and n respectively.")
    nblocks_per_col = m//bm
    nblocks_per_row = n//bn
    nblocks = nblocks_per_col*nblocks_per_row
    # 1D possible nz block indices
    BI = list(range(nblocks))
    # choose bnnz ones among them
    nzBI = []
    for i in range(bnnz):
        ri = randint(0, len(BI)-1)
        nzBI.append(BI[ri])
        del BI[ri]
    nzBI.sort()
    indices = np.array(nzBI) % nblocks_per_row
    bcolinds_ind = 0
    bindptr_ind = 1
    indptr = np.zeros((int(nblocks_per_col+1)))
    for bi in nzBI:
        while bi//nblocks_per_row+1 > bindptr_ind:
            bindptr_ind += 1
        indices[bcolinds_ind] = bi % nblocks_per_row
        bcolinds_ind += 1
        indptr[bindptr_ind] += 1
    for i in range(1, int(nblocks_per_col)+1):
        indptr[i] += indptr[i-1]
    data = rand(bnnz, bm, bn)
    return data, indices, indptr


def benchmark(bsrF, csrF, dsF, scipyF, M, label, nruns=10):
    """
    Measures time of Faust by M multiplication using the different Fausts in
    the first four arguments.
    """
    print('Benchmarking', label)
    times = np.zeros((4, nruns))
    Fausts = [bsrF, csrF, dsF, scipyF]
    for r in range(nruns):
        for i, F in enumerate(Fausts):
            print("run=", r, "i=", i, "F=", '\n', F)
            if isinstance(F, list):
                s = process_time()
                mul_scipy_faust(F, M)
                e = process_time()
            else:
                s = process_time()
                F@M
                e = process_time()
            times[i, r] = e-s
    np.savetxt(label+'.txt', times)


if __name__ == '__main__':
    # The Faust list of factors for the BSR, CSR and dense formats
    F_bsr_factors = []
    F_csr_factors = []
    F_ds_factors = []
    m = n = 1024
    bm = bn = 128
    bnnz = 20
    nfactors = 10
    for i in range(nfactors):
        data, indices, indtpr = rand_bsr(m, n, bm, bn, bnnz)
        F_bsr_factors += [bsr_matrix((data, indices, indtpr), shape=(m, n))]
        F_csr_factors.append(F_bsr_factors[i].tocsr())
        F_ds_factors.append(F_bsr_factors[i].toarray())
    bsrF = Faust(F_bsr_factors)
    csrF = Faust(F_csr_factors)
    dsF = Faust(F_ds_factors)
    scipyF = F_bsr_factors
    benchmark(bsrF, csrF, dsF, scipyF, rand(m, n), label='mul_dense')
    benchmark(bsrF, csrF, dsF, scipyF, random(m, n, .02, format='csr'),
              label='mul_csr')
    benchmark(bsrF, csrF, dsF, scipyF, rand(m), label='mul_vec')
