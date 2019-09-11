from __future__ import print_function
import sys
import os
if(len(sys.argv) > 1):
    sys.path.append(sys.argv[1])

from pyfaust import *
from numpy import empty, allclose, zeros, tril
from scipy.io import loadmat
from pyfaust.tools import UpdateCholeskyFull
from pyfaust.demo import get_data_dirpath

datap = os.path.join(get_data_dirpath(), 'faust_MEG_rcg_8.mat')

d = loadmat(datap)

facts = d['facts']
facts = [facts[0,i] for i in range(facts.shape[1]) ]
FD = Faust(facts)
D = FD.todense()

I = [125, 132, 1000, 155]
for P, Pt in [(lambda x: D*x, lambda x: D.H*x),
              (lambda x: np.matrix(FD*x),
               lambda x: np.matrix(FD.H*x))]:
    R = empty((0,0))
    for i in range(1, len(I)+1):
        R = UpdateCholeskyFull(R[0:i,0:i], P, Pt, I[:i], 8193)
        print(R, end='\n\n')
        assert(allclose(D[:,I[:min(i,len(I))]].H*D[:,I[:min(i,len(I))]], R.H*R))
        assert(allclose((R.H*R).H, R.H*R))
        assert(allclose(tril(R, -1), zeros(R.shape)))
        print(D[:,I].H*D[:,I])
        print(R.H*R)

