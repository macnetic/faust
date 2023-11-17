from __future__ import print_function
import sys
import os
if(len(sys.argv) > 1):
    sys.path.append(sys.argv[1])

from pyfaust import *
from numpy import empty, allclose, zeros, tril
from scipy.io import loadmat
from pyfaust.tools import (UpdateCholeskyFull, UpdateCholesky,
                           UpdateCholeskySparse)
from pyfaust.demo import get_data_dirpath
from scipy.sparse import issparse

if __name__ == '__main__':
    datap = os.path.join(get_data_dirpath(), 'faust_MEG_rcg_8.mat')

    d = loadmat(datap)

    facts = d['facts']
    facts = [facts[0,i] for i in range(facts.shape[1]) ]
    FD = Faust(facts)
    D = FD.toarray()

    I = [125, 132, 1000, 155]
    for P, Pt in [(lambda x: D@x, lambda x: D.T.conj()@x),
                  (lambda x: np.matrix(FD@x),
                   lambda x: np.matrix(FD.H@x))]:
        for s, R in enumerate([empty((0, 0)), csr_matrix((0, 0))]):
            for i in range(1, len(I)+1):
                R = UpdateCholesky(R[:i,:i], P, Pt, I[:i], 8193)
                print("R:", R.shape, issparse(R))
                print(R, end='\n\n')
                if s and issparse(R):
                    R = R.toarray()
                assert allclose(D[:,I[:min(i,len(I))]].T.conj()@D[:,I[:min(i,len(I))]],
                                R.T.conj()@R)

                assert allclose((R.T.conj()@R).T.conj(), R.T.conj()@R)
                assert allclose(tril(R, -1), zeros(R.shape))
                print(D[:,I].T.conj()@D[:,I])
                print(R.T.conj()*R)
                if s:
                    R = csr_matrix(R)
