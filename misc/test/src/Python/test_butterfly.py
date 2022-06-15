from pyfaust.fact import butterfly
from pyfaust import rand_butterfly, Faust, wht
from pyfaust.tools import bitrev_perm
import numpy as np
#import utils
from scipy.sparse import csr_matrix
import re
try :
    from hierarchical_fact import project_BP_model_P_fixed
    using_tung_leons_code = True
except:
    using_tung_leons_code = False

def err(a, F):
    if isinstance(F, Faust):
        err = (F-a).norm()/Faust(a).norm()
    else:
        err = np.linalg.norm(F-a)/ np.linalg.norm(a)
    return ("small error" if err < 1e-6 else "large error", err)

def is_as_expected(err, expected):
    return err < 1e-12 and re.match('.*small.*', expected) \
       or err > 1e-12 and re.match('.*large.*', expected)

def rand_perm(n):
    J = np.random.permutation(n)
    P = np.zeros((n, n))
    P[J, np.arange(n)] = 1
    P = csr_matrix(P)
    return P

def perm2indices(P):
    return P.T.nonzero()[1]

def indices2perm(I):
    n = len(I)
    P = np.zeros((n, n))
    P[I, np.arange(n)] = 1
    P = csr_matrix(P)
    return P

N = 32
em = ['expectation: should be small', 'expectation: should be large']
expectations = [em[0], em[1], em[1], em[1], em[0], em[1], em[1], em[1], em[0]]

for P, ptype in [(bitrev_perm(N), 'bitrev_perm'), (rand_perm(N),
                                                    'rand_perm')]:
    for H, mat_type in [(rand_butterfly(N), 'rand_butterly_real'),
                     (rand_butterfly(N, dtype='complex'), 'rand_butterly_complex'), (wht(N), 'hadamard')]:
        print("="*20, "matrix:", mat_type, '/', "P=", ptype)
        assert(np.allclose(indices2perm(perm2indices(P)).toarray(),
                           P.toarray()))
        Ha = H.toarray()
        P = P.astype(H.dtype)
        #P = csr_matrix(utils.bit_reversal_permutation_matrix(int(np.log2(N))))
        G1 = H@P
        G2 = H@P.T


        ref_F = np.empty((9), dtype=object)
        F = np.empty((9), dtype=object)
        for i in range(len(ref_F)):
            if i < 3:
                print("mat=H")
                mat = Ha
            elif i < 6:
                print("mat=H@P")
                mat = G1
            else:
                print("mat=H@P.T")
                mat = G2
            if i % 3 == 0:
                print("perm=None")
                p = None
            elif i % 3 == 1:
                print("perm=P")
                p = P
            else:
                print("perm=P.T")
                p = P.T
#            print("mat.flags['F_CONTIGUOUS']=", mat.flags['F_CONTIGUOUS'])
            if using_tung_leons_code:
                F_prod, F_factors = project_BP_model_P_fixed(mat, 'balanced',
                                                             p=(p.toarray() if p is not None else p), return_factors=True)
                ref_F[i] = Faust(F_factors)
            F[i] = butterfly(mat, type='bbtree', perm=(perm2indices(p) if p is not None else p))
            if using_tung_leons_code:
                if p is not None:
                    assert(np.allclose(ref_F[i]@p, F_prod, rtol=1e-6))
                else:
                    assert(np.allclose(ref_F[i].toarray(), F_prod, rtol=1e-6))
                print("ref_F"+str(i+1)+" err", err(mat, F_prod), expectations[i])
            print("F"+str(i+1)+" err", err(mat, F[i]), expectations[i])
            if not is_as_expected(err(mat, F[i])[1], expectations[i]):
                np.savez(mat_type+'mat-'+ptype+'-'+str(i+1)+'.npz', mat=mat,
                         perm=p, P=P)
                F[i].save('F'+str(i+1)+'.mat')
                print("(NOT OK)")
            print("---")
