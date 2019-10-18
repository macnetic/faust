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
from numpy.linalg import eig, norm
from time import clock
from os.path import exists
from scipy.sparse import spdiags
from scipy.sparse import csr_matrix
from scipy.io import savemat
#sys.path.append(sys.environ['HOME']+'/givens-factorization') # TODO: output arg
#from util import symmetrized_norm
from lapsolver import solve_dense

graph_ctors = [ erdosrenyi.ErdosRenyi, community.Community, path.Path,
               sensor.Sensor, randomring.RandomRing, ring.Ring ]

graph_names = [ 'erdosrenyi', 'community', 'path', 'sensor', 'randring', 'ring'
               ]

dim_size = 128

plt.rcParams['lines.markersize'] = .7


nruns = 17
plotting = False

# types of data for the benchmark
FOURIER_ERR_ID = 0
LAP_ERR_ID = 1
TIME_ID = 2

data_type_str = [ 'Fourier Error', 'Laplacian Error', 'Execution Time (sec)' ]

# FGFT algos tested
GIVENS = 0
PARGIVENS = 1
HIERPALM = 2
COORDDESCENT_TFRERIX = 3

algo_str = [ 'Givens', 'Parallel Givens', 'Hierarchical PALM', 'L1 Coordinate'
            ' Descent' ]

all_data = empty((len(graph_ctors), 3, 4, nruns), dtype=np.float) # 1st 3 is
                                                                  # for type of data: times, fourier_errs, lap_errs
                                                                  # 3 is for algos used


#times = all_data[:, : ,TIME_ID, :]
#fourier_errs = all_data[:, :, FOURIER_ERR_ID, :]
#lap_errs = all_data[:, :, LAP_ERR_ID, :]

data_file = 'benchmark_FGFT_pyfaust.npz'


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
            C[i,[j,j*2]] = np.sum(diff1), np.sum(diff2)
    return C

def symmetrized_norm(U, V): 
    C = compute_cost_matrix(U, V)
    row_idx, col_idx = solve_dense(C)
    print("linear sum prob sol:", row_idx, col_idx)
    best_frobenius_norm = np.sqrt(C[row_idx, col_idx].sum())
    return best_frobenius_norm

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
        D, U = eig(Lap)
        dim = Lap.shape[0]
        nfacts = int(round(log2(dim))-3)
        over_sp = 1.5 # sparsity overhead
        dec_fact = .5 # decrease of the residum sparsity
        fact_cons, res_cons = [], []
        for k in range(1, nfacts):
            fact_cons += [ ConstraintInt('sp', dim, dim,
                                         min(int(round(dec_fact**k*dim**2*over_sp)),Lap.shape[0]))
                         ]
            res_cons += [ ConstraintInt('sp', dim, dim,
                                        min(int(round(2*dim*over_sp)),Lap.shape[0])) ]

        params = ParamsHierarchical(fact_cons, res_cons,
                                    StoppingCriterion(num_its=50),
                                    StoppingCriterion(num_its=100),
                                    step_size=1.0000e-06,
                                    constant_step_size=False,
                                    init_lambda=1.0,
                                    is_fact_side_left=True,
                                    is_verbose=False)

        #t = clock()
        print("running Hierarchical PALM4MSA on U")
        try:
            F = pyfaust.fact.hierarchical(U,params)
            pass
        except:
            # PALM failing on this Laplacian, repeat this iteration on another
            continue
        #times[i,TIME_ID,HIERPALM,j] = clock()-t
        complexity_global = F.nnz_sum()

        #diag_init_D = diag(D)
        print("running Hierarchical PALM4MSA for FGFT")
        diag_init_D = np.copy(D)
        t = clock()
        try:
            F, Dhat = pyfaust.fact.fgft_palm(U, Lap,
                                         params,
                                         init_D=diag_init_D)
        except:
            # PALM failing on this Laplacian, repeat this iteration on another
            # one
            continue
        t = clock()-t
        all_data[i,TIME_ID,HIERPALM,j] = t
        hier_fgft_err = norm((F.todense()*diag(Dhat))*F.T.todense()-Lap,"fro")/norm(Lap,"fro")
        all_data[i,LAP_ERR_ID,HIERPALM,j]= hier_fgft_err
        all_data[i,FOURIER_ERR_ID,HIERPALM,j] = symmetrized_norm(U,
                                                                 F.toarray())/norm(U,'fro') #(F-U).norm("fro")/norm(U,"fro")





        #nfacts = complexity_global/(2*dim)
        print("running Trunc. Jacobi")#J = round(nfacts*(dim/2))
        J = round(complexity_global/4)
        #print("J=", J)
        #J = 64*64
        #J = Lap.shape[0]//2*(Lap.shape[0]-1)
        t = clock()
        Dhat, F = pyfaust.fact.fgft_givens(Lap, J, nGivens_per_fac=0)
        t = clock()-t
        Dhat = spdiags(Dhat, [0], Lap.shape[0], Lap.shape[0])
        all_data[i,TIME_ID,GIVENS,j] = t

        givens_err = norm((F*Dhat.todense())*F.T.todense()-Lap,'fro')/norm(Lap,'fro')
        all_data[i,LAP_ERR_ID,GIVENS,j] = givens_err
        all_data[i,FOURIER_ERR_ID,GIVENS,j] = symmetrized_norm(U,
                                                               F.toarray())/norm(U,'fro')
        print("errs Uhat permuted, Uhat:", symmetrized_norm(U,
                                                           F.toarray())/norm(U,'fro'),
             (F-U).norm('fro')/norm(U,'fro'))

# (F-U).norm("fro")/norm(U,"fro")

        print("running Parallel Trunc. Jacobi")
        t = clock()
        Dhat, F = pyfaust.fact.fgft_givens(Lap, J, nGivens_per_fac=int(dim/2))
        t = clock()-t
        Dhat = spdiags(Dhat, [0], Lap.shape[0], Lap.shape[0])
        all_data[i,TIME_ID,PARGIVENS,j] = t
        print("J=", J)
        print("nnz_sum FGFT givens parallel=", F.nnz_sum(), "num of facts:",
              F.numfactors())
        all_data[i,LAP_ERR_ID,PARGIVENS,j] = norm((F*Dhat.todense())*F.T.todense()-Lap,'fro')/norm(Lap,'fro')
        all_data[i,FOURIER_ERR_ID,PARGIVENS,j] = symmetrized_norm(U,
                                                                  F.toarray())/norm(U,'fro')
        print("errs Uhat permuted, Uhat:", symmetrized_norm(U,
                                                           F.toarray())/norm(U,'fro'),
             (F-U).norm('fro')/norm(U,'fro'))

 # (F-U).norm("fro")/norm(U,"fro")
        print("err:", norm(Lap @ F.toarray() - F.toarray()@Dhat)/norm(Lap))

#        # Thomas Frerix L1 Coordinate Descent algorithm
#        f = lambda x: norm(x,1)
#        from coordinate_descent import coordinate_descent
#        import scipy.sparse
#        print("running L1 coordinate_descent")
#        t = clock()
#        res_tuple = coordinate_descent(U,f, J)
#        t = clock()-t
#        # construct Givens factors and the Faust from them
#        Gs = []
#        for alpha,p,q in list(reversed(res_tuple)):
#            G = scipy.sparse.eye(*U.shape).tocsr()
#            G[p,p] = np.cos(alpha)
#            G[p,q] = -np.sin(alpha)
#            G[q,p]= np.sin(alpha)
#            G[q,q] = np.cos(alpha)
#            Gs += [G]
#        print("Gs=")
#        Gs = Faust(Gs)
#        print(Gs)
#        Gs.save('Gs.mat')
#        dU = { 'U': U }
#        savemat('U.mat', dU)
#        all_data[i,TIME_ID,COORDDESCENT_TFRERIX,j] = t
#        all_data[i,LAP_ERR_ID,COORDDESCENT_TFRERIX,j] = 1 # not computed
#        abs_err = symmetrized_norm(U, Gs.toarray())
#        print("sym_norm:", abs_err, "rel err", abs_err/norm(U,'fro'), "err"
#              " without permutation:", (Gs-U).norm('fro')/norm(U, 'fro'))
#        all_data[i,FOURIER_ERR_ID,COORDDESCENT_TFRERIX,j] = abs_err/norm(U,"fro")

        savez(data_file, all_data[:,:,:,:j+1])
        i += 1

plt.rcParams['figure.figsize'] = [12.0, 8]
for j in range(0,3): # j : index for type of data
    f, axes = plt.subplots(1, 3, sharey=True)
    plt.suptitle("pyfaust Benchmark on "+data_type_str[j]+"\n"+str(nruns)+" runs "
                 "(Laplacian size: "+str(dim_size)+'^2)', size=14)
    for i in range(0,3): # i : index for type of algo
        # k is for type of graph
        # plt.figure()
        axes[i].set_title("Algo: "+ algo_str[i])
        if(j == TIME_ID):
            axes[i].semilogy()
        axes[i].boxplot([all_data[k,j,i, :] for k in
                     range(0,len(graph_ctors))],showfliers=True)
        axes[i].set_xticklabels(list(graph_names[:]), rotation='vertical')

    axes[0].set_ylabel(data_type_str[j])
    plt.savefig('pyfaust_FGFT_benchmark-'+data_type_str[j]+'.png')
plt.show(block=True)


#G = erdosrenyi.ErdosRenyi(N=64, seed=42)
#G.set_coordinates(kind='spring', seed=42)
#fig, axes = plt.subplots(1, 2)
#_ = axes[0].spy(G.W, markersize=2)
#_ = G.plot(ax=axes[1])

