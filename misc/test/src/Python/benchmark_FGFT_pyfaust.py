import pygsp.graphs
from pygsp.graphs import community, erdosrenyi, path, sensor, randomring, ring
import matplotlib.pyplot as plt
import random
from random import randint
import numpy as np
from numpy import empty,log2, size, count_nonzero, diag, savez, load
from pyfaust import *
from pyfaust.factparams import *
from numpy.linalg import eig, norm
from time import clock
from os.path import exists


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

algo_str = [ 'Givens', 'Parallel Givens', 'Hierarchical PALM' ]

all_data = empty((len(graph_ctors), 3, 3, nruns), dtype=np.float) # 1st 3 is
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

        params = ParamsHierarchicalFact(fact_cons, res_cons,
                                        StoppingCriterion(num_its=50),
                                        StoppingCriterion(num_its=100),
                                        step_size=1.0000e-06,
                                        constant_step_size=False,
                                        init_lambda=1.0,
                                        is_fact_side_left=True, is_verbose=True)

        #t = clock()
        try:
            F = FaustFactory.fact_hierarchical(U,params)
        except:
            # PALM failing on this Laplacian, repeat this iteration on another
            continue
        #times[i,TIME_ID,HIERPALM,j] = clock()-t
        complexity_global = F.nnz_sum()

        #diag_init_D = diag(D)
        diag_init_D = np.copy(D)
        t = clock()
        try:
            F, Dhat = FaustFactory.fgft_palm(U, Lap,
                                         params,
                                         init_D=diag_init_D)
        except:
            # PALM failing on this Laplacian, repeat this iteration on another
            continue
        t = clock()-t
        all_data[i,TIME_ID,HIERPALM,j] = t
        hier_fgft_err = norm((F.todense()*diag(Dhat))*F.T.todense()-Lap,"fro")/norm(Lap,"fro")
        all_data[i,LAP_ERR_ID,HIERPALM,j]= hier_fgft_err
        all_data[i,FOURIER_ERR_ID,HIERPALM,j] =  (F-U).norm("fro")/norm(U,"fro")

        #nfacts = complexity_global/(2*dim)
        #J = round(nfacts*(dim/2))
        J = round(complexity_global/4)
        t = clock()
        F, Dhat = FaustFactory.fgft_givens(Lap, J, 0)
        t = clock()-t
        all_data[i,TIME_ID,GIVENS,j] = t

        givens_err = norm((F*Dhat.todense())*F.T.todense()-Lap,'fro')/norm(Lap,'fro')
        all_data[i,LAP_ERR_ID,GIVENS,j] = givens_err
        all_data[i,FOURIER_ERR_ID,GIVENS,j] = (F-U).norm("fro")/norm(U,"fro")

        t = clock()
        F, Dhat = FaustFactory.fgft_givens(Lap, J, int(dim/2))
        t = clock()-t
        all_data[i,TIME_ID,PARGIVENS,j] = t
        print("J=", J)
        print("nnz_sum FGFT givens parallel=", F.nnz_sum(), "num of facts:",
              F.get_num_factors())
        all_data[i,LAP_ERR_ID,PARGIVENS,j] = norm((F*Dhat.todense())*F.T.todense()-Lap,'fro')/norm(Lap,'fro')
        all_data[i,FOURIER_ERR_ID,PARGIVENS,j] = (F-U).norm("fro")/norm(U,"fro")
        print("err:", norm(Lap * F.toarray() - F.toarray()*Dhat)/norm(Lap))
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

