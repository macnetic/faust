import pygsp.graphs
from pygsp.graphs import community, erdosrenyi, path, sensor, randomring, ring
import matplotlib.pyplot as plt
import random
from random import randint
import numpy as np
from numpy import empty,log2, size, count_nonzero, diag, savez, load
from pyfaust import *
import pyfaust.fact
from pyfaust.factparams import *
from numpy.linalg import eigh, norm, eig
from time import clock
from os.path import exists
from scipy.sparse import spdiags
from scipy.sparse import csr_matrix
from scipy.io import savemat
#sys.path.append(sys.environ['HOME']+'/givens-factorization') # TODO: output arg
#from util import symmetrized_norm
from lapsolver import solve_dense
import math

graph_ctors = [ erdosrenyi.ErdosRenyi, community.Community, path.Path,
               sensor.Sensor, randomring.RandomRing, ring.Ring ]

graph_names = [ 'erdosrenyi', 'community', 'path', 'sensor', 'randring', 'ring'
               ]

dim_size = 128

plt.rcParams['lines.markersize'] = .7


nruns = 20
plotting = False

# types of data for the benchmark
FOURIER_ERR_ID = 0
LAP_ERR_ID = 1
D_ERR_ID = 2
TIME_ID = 3

data_type_str = [ 'Fourier Error', 'Laplacian Error', 'D error', 'Execution Time (sec)' ]

# FGFT algos/confs tested
GIVENS_REAL = 0
GIVENS_REAL_SPARSE = 1
GIVENS_REAL_AS_CPLX = 2
GIVENS_REAL_AS_CPLX_SPARSE = 3

PARGIVENS_REAL = 4
PARGIVENS_REAL_SPARSE = 5
PARGIVENS_REAL_AS_CPLX = 6
PARGIVENS_REAL_AS_CPLX_SPARSE = 7

GIVENS_HERMIT = 8
GIVENS_HERMIT_SPARSE = 9
PARGIVENS_HERMIT = 10
PARGIVENS_HERMIT_SPARSE = 11

algo_str = [ 'Givens (real dense)', 'Givens (real sparse)'
            , 'Givens (cplx dense)',  'Givens (cplx sparse)']
algo_str += [ 'Givens //(real dense)', 'Givens //(real sparse)'
            , 'Givens //(cplx dense)',  'Givens //(cplx sparse)' ]

algo_str += [ 'Givens (hermit. dense)', 'Givens (hermit. sparse)',
             'Givens // (hermit. dense)', 'Givens // (hermit. sparse)' ]

all_data = empty((len(graph_ctors), 4, len(algo_str), nruns), dtype=np.float) # arg 4 is
                                                                  # for type of
                                                                  # data:
                                                                      # times,
                                                                      # fourier_errs,
                                                                      # lap_errs,
                                                                      # D err


#times = all_data[:, : ,TIME_ID, :]
#fourier_errs = all_data[:, :, FOURIER_ERR_ID, :]
#lap_errs = all_data[:, :, LAP_ERR_ID, :]

data_file = 'benchmark_eigtj_pyfaust.npz'


# get already computed data from file
old_nruns = 0
if(exists(data_file)):
    od = load(data_file)['arr_0']
    all_data[:od.shape[0],:od.shape[1],:od.shape[2],:od.shape[3]] = od
    old_nruns = od.shape[3]


J = 0

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

def best_permutation(U, V, D):
    C = compute_cost_matrix(U, V)
    row_idx, col_idx = solve_dense(C)
    pV = np.zeros_like(U)
    pV = pV.astype(V.dtype)
    pD = np.zeros_like(D)
    for i,j in enumerate(col_idx):
        if(j >= U.shape[0]):
            j = j%U.shape[0]
            pV[:,i] = V[:,j]
            pD[i] = D[j]
        else:
            pV[:,i] = -V[:,j]
            pD[i] = pD[j]
    return pD, pV


def permu_near_id(U):
    I = [ i for i in range(0,U.shape[1]) ]
    J = []
    for j in range(0,U.shape[1]):
        J += [I[np.argmax(U[j,I])]]
        I.remove(J[-1])
    return U[:,J]


def hermitian_from_lap(Lap):
    """
    builds an hermitian matrix from the real Laplacian
    by multiplying nonzero coeffs by a random phase
    """
    H = Lap.astype(np.complex)
    for r in range(0,Lap.shape[0]):
        for c in range(0,r+1):
            if(Lap[r,c] != 0):
                theta = np.random.rand()*math.pi*2
                H[r,c] = Lap[r,c]*np.exp(theta*np.complex(0,1))
    H = np.tril(H) + np.tril(H).conj().T
    assert((H.conj().T == H).all())
    assert((H.imag != 0).any())
    return H

for j in range(old_nruns,nruns):
    print("================= #run:", j)
    i = 0
    while(i < len(graph_ctors)):
        print("=========================== (#run,family graph)", j,
              graph_names[i])
        kwargs = dict()
        if(not graph_ctors[i] in [community.Community, path.Path, ring.Ring ]):
            kwargs['seed'] = randint(1,10000);
        G = graph_ctors[i](N=dim_size, **kwargs)
        G.set_coordinates(**kwargs) # kind = 'spring'
        Lap = G.L.toarray().astype(np.float)
        #print(L)
        if((Lap != Lap.T).any()): exit("Lap. not valid")
        if(plotting):
            fig, axes = plt.subplots(1, 2)
            _ = axes[0].spy(G.W,markersize=.7)
            _ = G.plot(ax=axes[1], vertex_size=.9)
        D, U = eigh(Lap)
        dim = Lap.shape[0]

        print("running Trunc. Jacobi")
        #print("J=", J)
        #J = 64*64
        #J = Lap.shape[0]//2*(Lap.shape[0]-1)
        for a in range(0, 4):
            print(algo_str[a])
            J = round(dim_size*(dim_size-1)/2)
            Lapa = Lap
            if(a in [GIVENS_REAL_AS_CPLX, GIVENS_REAL_AS_CPLX_SPARSE]):
                Lapa = Lapa.astype(np.complex)
            if(a in [GIVENS_REAL_SPARSE, GIVENS_REAL_AS_CPLX_SPARSE]):
                Lapa = csr_matrix(Lapa)
            t = clock()
            Dhat, F = pyfaust.fact.eigtj(Lapa, J, nGivens_per_fac=0)
            t = clock()-t
            Dhata, Fa = best_permutation(U, F.toarray(), Dhat)
            all_data[i,D_ERR_ID,a, j] = norm(D-Dhat)/norm(D)
            Dhat = spdiags(Dhat, [0], Lapa.shape[0], Lapa.shape[0])
            all_data[i,TIME_ID,a,j] = t

            if(a in [GIVENS_REAL_SPARSE, GIVENS_REAL_AS_CPLX_SPARSE]):
                Lapa = Lapa.todense()
            givens_err = \
            norm(Fa@diag(Dhata)@Fa.T.conj()-Lap,'fro')/norm(Lap,'fro')
            givens_err2 = norm(F@Dhat.todense()@F.toarray().T.conj()-Lap,'fro')/norm(Lap,'fro')
            all_data[i,LAP_ERR_ID,a,j] = min(givens_err, givens_err2)
            print("lap err:", givens_err, givens_err2)
            all_data[i,FOURIER_ERR_ID,a,j] = symmetrized_norm(U, F.toarray())/norm(U,'fro')
            print("fourier err:",  all_data[i,FOURIER_ERR_ID,a,j])
            print("fourier err2:",
                  norm(permu_near_id(Fa.T.conj()@U)-np.eye(*(U.shape)))/norm(np.eye(*(U.shape))))

#            print("errs Uhat permuted, Uhat:", symmetrized_norm(U,
#                                                               Fa)/norm(U,'fro'),
#                 norm(Fa-U, 'fro')/norm(U,'fro'))
            #all_data[i,FOURIER_ERR_ID,a, j] = norm(U-Fa,'fro')/norm(U,'fro')

        print("running Parallel Trunc. Jacobi")
        for a in range(4, 8):
            print(algo_str[a])
            J = round(dim_size*(dim_size-1)/2)
            Lapa = Lap
            if(a in [PARGIVENS_REAL_AS_CPLX, PARGIVENS_REAL_AS_CPLX_SPARSE]):
                Lapa = Lapa.astype(np.complex)
            if(a in [PARGIVENS_REAL_SPARSE, PARGIVENS_REAL_AS_CPLX_SPARSE]):
                Lapa = csr_matrix(Lapa)

            t = clock()
            Dhat, F = pyfaust.fact.eigtj(Lap, J, nGivens_per_fac=int(dim/2))
            t = clock()-t
            Dhata, Fa = best_permutation(U, F.toarray(), Dhat)
            all_data[i,D_ERR_ID,a, j] = norm(D-Dhat)/norm(D)
            Dhat = spdiags(Dhat, [0], Lap.shape[0], Lap.shape[0])
            all_data[i,TIME_ID,a,j] = t
            if(a in [PARGIVENS_REAL_SPARSE, PARGIVENS_REAL_AS_CPLX_SPARSE]):
                Lapa = Lapa.todense()
            print("J=", J)
            pargivens_err1 = \
                    norm(Fa@diag(Dhata)@Fa.T.conj()-Lap,'fro')/norm(Lap,'fro')
            pargivens_err2 = norm(F@Dhat.todense()@F.toarray().T.conj()-Lap,'fro')/norm(Lap,'fro')
            all_data[i,LAP_ERR_ID,a,j] = min(pargivens_err1, pargivens_err2)
            print("lap err:", pargivens_err1, givens_err2)
            all_data[i,FOURIER_ERR_ID,a,j] = symmetrized_norm(U, F.toarray())/norm(U,'fro')
            print("fourier err:",  all_data[i,FOURIER_ERR_ID,a,j])
            print("fourier err2:",
                  norm(permu_near_id(Fa.T.conj()@U)-np.eye(*(U.shape)))/norm(np.eye(*(U.shape))))
            #print("errs Uhat permuted, Uhat:", symmetrized_norm(U, Fa)/norm(U,'fro'),
            #     norm(U-Fa, 'fro')/norm(U,'fro'))
            #all_data[i,FOURIER_ERR_ID,a, j] = norm(U-Fa,'fro')/norm(U,'fro')
            print("err:", norm(Lap @ F.toarray() - F.toarray()@Dhat)/norm(Lap))

        H = hermitian_from_lap(Lap)
        D, U = eigh(H)
        # benchmark truncated jacobi cplx algo on it
        for a in range(8, 12):
            if(a < 10):
                print("Running truncated jacobi on a hermitian matrix"
                      " (deduced from Laplacian)")
                nGivens_per_fac=0
            else:
                print("Running parallel truncated jacobi on a hermitian"
                      " matrix (deduced from Laplacian)")
                nGivens_per_fac=H.shape[0]//2
            print(algo_str[a])
            if(a in [GIVENS_HERMIT_SPARSE, PARGIVENS_HERMIT_SPARSE]):
                H = csr_matrix(H)
            J = round(dim_size*(dim_size-1)/2)
            t = clock()
            Dhat, F = pyfaust.fact.eigtj(H, J,
                                         nGivens_per_fac=nGivens_per_fac)
            t = clock()-t
            Dhata, Fa = best_permutation(U, F.toarray(), Dhat)
            D_err1 = norm(D-Dhat)/norm(D)
            D_err2 = norm(D-Dhata)/norm(D)
            print('D_err1, D_err2:', D_err1, D_err2)
            all_data[i,D_ERR_ID,a, j] = D_err1
            print("Greatest/Last eigenvalues of Dhat, Dhata, and D:", Dhat[-1],
                  Dhata[-1], D[-1])
            Dhat = spdiags(Dhat, [0], H.shape[0], H.shape[0])
            all_data[i,TIME_ID,a,j] = t

            if(a in [GIVENS_HERMIT_SPARSE, PARGIVENS_HERMIT_SPARSE]):
                H = H.todense()
            givens_err = \
            norm(Fa@diag(Dhata)@Fa.T.conj()-H,'fro')/norm(H,'fro')
            givens_err2 = \
            norm(F@Dhat.todense()@F.toarray().T.conj()-H,'fro')/norm(H,'fro')
            all_data[i,LAP_ERR_ID,a,j] = min(givens_err, givens_err2)
            print("lap err:", givens_err, givens_err2)
            all_data[i,FOURIER_ERR_ID,a,j] = symmetrized_norm(U, F.toarray())/norm(U,'fro')
            print("fourier err:",  all_data[i,FOURIER_ERR_ID,a,j])
            print("fourier err2:",
                  norm(permu_near_id(Fa.T.conj()@U)-np.eye(*(U.shape)))/norm(np.eye(*(U.shape))))
            #print("errs Uhat permuted, Uhat:", symmetrized_norm(U, Fa)/norm(U,'fro'),
            #     norm(U-Fa, 'fro')/norm(U,'fro'))
            #all_data[i,FOURIER_ERR_ID,a, j] = norm(U-Fa,'fro')/norm(U,'fro')
            print("err:", norm(H @ F.toarray() - F.toarray()@Dhat)/norm(H))


        savez(data_file, all_data[:,:,:,:j+1])
        i += 1

plt.rcParams['figure.figsize'] = [12.0, 8]
for j in range(0,4): # j : index for type of data
    for a in range(0,2): # a == 0 (Givens), a == 1 (Givens //)
        nsubplots = 4
        f, axes = plt.subplots(1, nsubplots, sharey=True)
        plt.suptitle("pyfaust eigtj benchmark on "+data_type_str[j]+"\n"+str(nruns)+" runs "
                     "(Laplacian size: "+str(dim_size)+'^2)'
                     ' J='+str(round(dim_size*(dim_size-1)/2)), size=14)
        for i in range(0,nsubplots): # i : index for type of algo
            # k is for type of graph
            # plt.figure()
            axes[i].set_title("Algo: "+ algo_str[i+a*nsubplots])
            if(j == TIME_ID):
                axes[i].semilogy()
            axes[i].boxplot([all_data[k,j,i+a*nsubplots, :] for k in
                             range(0,len(graph_ctors))],showfliers=True)
            axes[i].set_xticklabels(list(graph_names[:]), rotation='vertical')

        axes[0].set_ylabel(data_type_str[j])
        plt.savefig('pyfaust_eigtj_benchmark-'+data_type_str[j]+'-'+str(a)+'.png')

    # plots results for GIVENS_HERMIT, GIVENS_HERMIT_SPARSE, PARGIVENS_HERMIT,
    # PARGIVENS_HERMIT_SPARSE
    nsubplots = 4
    f, axes = plt.subplots(1, nsubplots, sharey=True)
    plt.suptitle("pyfaust eigtj benchmark on "+data_type_str[j]+"\n"+str(nruns)+" runs "
                 "(Laplacian size: "+str(dim_size)+'^2)'
                 ' J='+str(round(dim_size*(dim_size-1)/2)), size=14)
    for i in range(0,nsubplots): # i : offset for type of algo
        # k is for type of graph
        # plt.figure()
        axes[i].set_title("Algo: "+ algo_str[8+i])
        if(j == TIME_ID):
            axes[i].semilogy()
        axes[i].boxplot([all_data[k,j,8+i, :] for k in
                         range(0,len(graph_ctors))],showfliers=True)
        axes[i].set_xticklabels(list(graph_names[:]), rotation='vertical')

    axes[0].set_ylabel(data_type_str[j])
    plt.savefig('pyfaust_eigtj_benchmark-'+data_type_str[j]+'-'+str(2)+'.png')

plt.show(block=True)


#G = erdosrenyi.ErdosRenyi(N=64, seed=42)
#G.set_coordinates(kind='spring', seed=42)
#fig, axes = plt.subplots(1, 2)
#_ = axes[0].spy(G.W, markersize=2)
#_ = G.plot(ax=axes[1])

