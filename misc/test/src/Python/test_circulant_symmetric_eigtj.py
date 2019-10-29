from pyfaust.fact import eigtj
from numpy.fft import fft
import numpy as np
from numpy import eye, sqrt, diag, size
from scipy.linalg import toeplitz
from numpy.linalg import norm
import sys
from lapsolver import solve_dense

def compute_cost_matrix(U, V):
    d = U.shape[0]
    C = np.empty((d,2*d))
    for i in range(0, d):
        for j in range(0, d):
            diff1 = U[:,i] + V[:,j]
            diff2 = U[:,i] - V[:,j]
            diff1 *= diff1
            diff2 *= diff2
            C[i,[j,d+j]] = np.sum(abs(diff1)), np.sum(abs(diff2))
    return C

def symmetrized_norm(U, V):
    C = compute_cost_matrix(U, V)
    row_idx, col_idx = solve_dense(C)
    #print("linear sum prob sol:", row_idx, col_idx)
    best_frobenius_norm = np.sqrt(abs(C[row_idx, col_idx].sum()))
    return best_frobenius_norm

def best_permutation(U, V):
    C = compute_cost_matrix(U, V)
    row_idx, col_idx = solve_dense(C)
    pV = np.zeros_like(U)
    pV = pV.astype(V.dtype)
    for i,j in enumerate(col_idx):
        if(j >= U.shape[0]):
            j = j%U.shape[0]
            pV[:,i] = V[:,j]
        else:
            pV[:,i] = -V[:,j]
    return pV

if(len(sys.argv) > 1):
    n = int(sys.argv[1])
else:
    n = 128

F = fft(eye(n))/sqrt(n)

C = [ np.random.rand() for i in range(0,n//2-1) ]
C = [ np.random.rand() ] + C[0:n//2-1] + [np.random.rand()] + C[n//2-2::-1]


T = toeplitz(C,C)


# verify T is circulant
for i in range(1,n-1):
    d1 = diag(T, i)
    for j in range(0,size(d1)-1):
        assert(d1[j]==d1[j+1])
    d2 = diag(T, n-i)
    for j in range(0,size(d2)-1):
        assert(d2[j]==d2[j+1])
    e1 = d1[0]
    e2 = d2[0]
    assert(e1==e2)

# verify T is symmetric
assert((T == T.T).all())

Dhat, Uhat = eigtj(T.astype(np.complex), int(np.floor(n*np.log2(n))), nGivens_per_fac=n//2)

print("Dhat=", Dhat)
print("Uhat=", Uhat)

print("Uhat.toarray()=", Uhat.toarray())
print("F=", F)
print("err:", norm(Uhat.toarray()-F)/norm(F))
print("err best permu.:", symmetrized_norm(F,Uhat.toarray())/norm(F))
print("err best permu. (verif):", norm(best_permutation(F,Uhat.toarray())-F)/norm(F))
print("norm(F):", norm(F, 2))
print("Uhat.norm():", Uhat.norm(2))


#from itertools import permutations
#from numpy.random import permutation
#min_err=1000
#min_p = []
#for p in permutations(permutation(n)):
#    p = list(p)
#    err = norm(F-Uhat[:,p].toarray())/norm(F)
#    if(err < min_err):
#        min_err = err
#        min_p = p
#print(p)
#print(min_err)
