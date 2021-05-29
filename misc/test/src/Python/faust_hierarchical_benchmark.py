from pyfaust.fact import hierarchical_py
from pyfaust.fact import hierarchical, hierarchical2020
from pyfaust.factparams import ParamsHierarchicalSquareMat
from pyfaust import wht, Faust
import numpy as np
from numpy.linalg import norm
from pyfaust.proj import splincol
from timeit import timeit
import os
import sys
import math
import matplotlib.pyplot as plt
from qkmeans.core.utils import build_constraint_set_smart
from qkmeans.data_structures import SparseFactors
import time
import logging
import daiquiri
from numpy import savez, load

from qkmeans.palm.palm_fast import hierarchical_palm4msa
#from qkmeans.palm.palm import hierarchical_palm4msa
from qkmeans.core.utils import build_constraint_set_smart


# loop on all pyfaust implementations of hierarchical Palm4MSA algorithm
# use default values for norm2 parameters if environement variables
# NORM2_THRESHOLD or NORM2_MAX_ITER are not available
# store times in array and errors in another

# one function for each implementation that returns tuple of time, error and
# name

def benchmark_hierarchical2020(H):
    p = ParamsHierarchicalSquareMat.createParams(H, 'squaremat')
    err = [math.inf]
    packing_RL = False
    if('PACKING_RL' in os.environ):
        packing_RL=bool(int(os.environ['PACKING_RL']))
    def exec_algo():
        p.norm2_max_iter = NORM2_MAX_ITER
        p.norm2_threshold = NORM2_THRESHOLD
        p.packing_RL = False
        print("packing_RL:", p.packing_RL)
        print("factor_format:", p.factor_format)
#        F, _lambda = hierarchical2020(H, p.stop_crits[0].num_its, p.constraints,
#                                      p.is_update_way_R2L, p.is_fact_side_left, True,
#                                      p.norm2_threshold, p.norm2_max_iter,
#                                      packing_RL)
        F_faust, _lambda = hierarchical(H, p, ret_lambda=True, backend=2020)
        #F_faust = F*_lambda
        err[0] = error(H,F_faust)
    time = timeit(exec_algo, number=1)
    name = "hierarchical2020"
    if(packing_RL):
        name += "p"
    else:
        name += "np"
    return time, err[0], name

def benchmark_hierarchical2016(H):
    p = ParamsHierarchicalSquareMat.createParams(H, 'squaremat')
    err = [math.inf]
    def exec_algo():
        p = ParamsHierarchicalSquareMat.createParams(H, 'squaremat')
        #p.is_verbose = True
        #F_faust = hierarchical(H, "squaremat")
        #p.norm2_max_iter = 10000
        #p.norm2_threshold = 1e-16
        p.norm2_max_iter = NORM2_MAX_ITER
        p.norm2_threshold = NORM2_THRESHOLD
        F_faust = hierarchical(H, p)
        err[0] = error(H,F_faust)
    time = timeit(exec_algo, number=1)
    return time, err[0], "hierarchical2016"

def benchmark_hierarchical2020_py(H):
    if('NORM2_ON_ARRAY' in os.environ.keys()):
        compute_2norm_on_array=True
    else:
        compute_2norm_on_array=False
    print("hierarchical_py")
    print("compute_2norm_on_array=", compute_2norm_on_array)
    err = [math.inf]
    def exec_algo():
        if(H.shape[0] == 1024): # fix
            compute_2norm_on_arrays=True
        F_faust = hierarchical_py(H, J=int(np.log2(dim)),
                              N=30,fac_proxs=[splincol((dim,dim), 2,
                                                       normalized=True) for i in
                                              range(l2dim-1)],
                              res_proxs=[splincol((dim,dim), int(dim/2**(i+1)),
                                                 normalized=True) for i in
                                         range(l2dim-1)], is_update_way_R2L=True,
                              is_fact_side_left=False,
                              compute_2norm_on_arrays=compute_2norm_on_array,
                                  norm2_max_iter=NORM2_MAX_ITER,
                                  norm2_threshold=NORM2_THRESHOLD)
        err[0] = error(H,F_faust)
    time = timeit(exec_algo, number=1)
    return time, err[0], "hierarchical2020_py"

def benchmark_pyqalm_hierarchical(H):
    err = [math.inf]
    def exec_algo():
        dim = H.shape[0]
        nb_factors = int(np.log2(dim))
        sparsity_factor = 2
        lst_factors = [np.eye(dim) for _ in range(nb_factors)]
        lst_factors[-1] = np.zeros((dim, dim))
        _lambda = 1
        nb_iter = 30

        lst_proj_op_by_fac_step, lst_proj_op_by_fac_step_desc = build_constraint_set_smart(left_dim=dim,
                                                                                               right_dim=dim,
                                                                                               nb_factors=nb_factors,
                                                                                               sparsity_factor=sparsity_factor,
                                                                                               residual_on_right=True,
                                                                                               fast_unstable_proj=False, constant_first=False)



        final_lambda, final_factors, final_X, _, _ = hierarchical_palm4msa(
            arr_X_target=H,
            lst_S_init=lst_factors,
            lst_dct_projection_function=lst_proj_op_by_fac_step,
            f_lambda_init=_lambda,
            nb_iter=nb_iter,
            update_right_to_left=True,
            residual_on_right=True, track_objective_palm=False)

        if(isinstance(final_factors, SparseFactors)):
            F = Faust([final_factors._lst_factors[i] for i in range(len(final_factors))])
        else:
            F = Faust([final_factors[i] for i in range(len(final_factors))])
        err[0] = norm(final_X-H)/norm(H)
    time = timeit(exec_algo, number=1)
    return time, err[0], "hierarchical_pyqalm"


def error(H, FH):
    assert(isinstance(FH, Faust))
    return (FH-H).norm()/norm(H)

# at most 4 algos tested
markers= [ '+', '^', '*', 'x' ]

if(len(sys.argv) <= 1):
    print("""USAGE:
          """+sys.argv[0]+""" benchmark_res.npz|<log2dim_min>[<log2dim_max>] """)
    exit(1)

try:
    l2dim0 = int(sys.argv[1])
    l2dim1 = l2dim0
    if(len(sys.argv) > 2):
        l2dim1 = int(sys.argv[2])


    benchmarks = [ benchmark_pyqalm_hierarchical,benchmark_hierarchical2020,
                  benchmark_hierarchical2016, ]

    if(not "NO_HIERAR_PY" in os.environ or os.environ['NO_HIERAR_PY'] != '1'):
        benchmarks += [ benchmark_hierarchical2020_py ]

    benchmarks = np.array(benchmarks)

    names = np.empty(len(benchmarks), dtype=list)
    times = np.empty((benchmarks.shape[0], l2dim1-l2dim0+1))
    errs = np.empty((benchmarks.shape[0], l2dim1-l2dim0+1))

    if('NORM2_THRESHOLD' in os.environ):
        NORM2_THRESHOLD = float(os.environ['NORM2_THRESHOLD'])
    else:
        NORM2_THRESHOLD = 1e-16
    if('NORM2_MAX_ITER' in os.environ):
        NORM2_MAX_ITER = int(os.environ['NORM2_MAX_ITER'])
    else:
        NORM2_MAX_ITER = 100



    for dim_index, l2dim in enumerate(range(l2dim0, l2dim1+1)):
        dim = 2**l2dim
        H = wht(dim).toarray()
        for i,bench in enumerate(benchmarks):
            time, err, name = bench(H)
            times[i, dim_index] = time
            errs[i, dim_index] = err
            names[i] = name
            print(name, "time:", time, "err:", err)

    savez("faust_hierarchical_benchmark_results.npz", **{ 'times': times, 'errs':
                                                       errs, 'names': names,
                                                         'l2dims':
                                                         np.array([l2dim0,
                                                                   l2dim1]),
                                                         'norm_params':
                                                         np.array([NORM2_MAX_ITER,
                                                                  NORM2_THRESHOLD])})
except:
    # expected error: sys.argv[1] is not an int but a filename
    # this is a pre-computed set of data times, errs and names to render
    d = load(sys.argv[1], allow_pickle=True)
    times = d['times']
    errs = d['errs']
    names = d['names']
    l2dim0 = int(d['l2dims'][0])
    l2dim1 = int(d['l2dims'][1])
    NORM2_MAX_ITER, NORM2_THRESHOLD = d['norm_params'][:]

plt.rcParams['figure.figsize'] = [12.0, 8]

assert(len(markers) >= len(names))

for bi in range(len(names)):
    if(names[bi] == "hierarchical2020_py" and "NO_HIERAR_PY" in os.environ and os.environ['NO_HIERAR_PY'] == '1'):
        print("skipping hierarchical_py")
        continue
    plt.plot(2**np.arange(l2dim0, l2dim1+1), times[bi,:], label=names[bi],
             marker=markers[bi])

plt.grid(True)
plt.legend()
plt.title("hierarchical times vs hadamard matrix dimension (norm2 "
          "max_iter="+str(NORM2_MAX_ITER)+
          "threshold="+str(NORM2_THRESHOLD)+")")
plt.ylabel("secs")
plt.tight_layout()
plt.xticks(2**np.arange(l2dim0, l2dim1+1))
plt.savefig(sys.argv[0]+"_times.png")
plt.show(block=False)

plt.figure()
for bi in range(len(names)):
    if(names[bi] == "hierarchical2020_py" and "NO_HIERAR_PY" in os.environ and os.environ['NO_HIERAR_PY'] == '1'):
        print("skipping hierarchical_py")
        continue
    plt.plot(2**np.arange(l2dim0, l2dim1+1), errs[bi,:], label=names[bi],
             marker=markers[bi])

plt.grid(True)
plt.legend()
plt.title("hierarchical errors vs hadamard matrix dimension (norm2 "
          "max_iter="+str(NORM2_MAX_ITER)+
          " threshold="+str(NORM2_THRESHOLD)+")")
plt.tight_layout()
plt.xticks(2**np.arange(l2dim0, l2dim1+1))
plt.savefig(sys.argv[0]+"_errors.png")
plt.show(block=True)

