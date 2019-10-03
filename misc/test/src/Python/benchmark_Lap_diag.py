from scipy.io import loadmat
from sys import argv
from os.path import dirname, sep, join, exists
from numpy import diag, copy, log2, count_nonzero, size, loadtxt, savetxt
import numpy, re
from numpy.linalg import norm
from pyfaust.factparams import *
from pyfaust import *
from pylab import *
from time import clock




if __name__ == '__main__':

    data_file = join(dirname(argv[0]),
                     '../../../data/mat/Laps_U_Ds-6_graph_types-dims_128_to_1024.mat')
    if(not exists(data_file)):
        raise Exception(data_file+" not found. You need to run"
                        " misc/test/src/Matlab/get_graph_test_UDLap_matrices.m"
                        " for generation.")
    matfile = loadmat(data_file)
    Laps = matfile['Laps'][0]
    Ds = matfile['Ds'][0]
    Us = matfile['Us'][0]
    start_i=0
    num_laps = Laps.shape[0]
    num_laps = 45
    hier_fgft_errs = empty(num_laps)
    givens_errs = empty(num_laps)
    par_givens_errs = empty(num_laps)

    U_hier_fgft_errs = empty(num_laps)
    U_givens_errs = empty(num_laps)
    U_par_givens_errs = empty(num_laps)

    hier_palm_times = empty(num_laps)
    hier_fgft_times = empty(num_laps)
    givens_times = empty(num_laps)
    par_givens_times = empty(num_laps)

    if exists('benchmark_pyfaust_fgft.txt'):
        saved_lap_errs = loadtxt('benchmark_pyfaust_fgft.txt')
        len_saved = saved_lap_errs[0].shape[0]
        givens_errs[:len_saved], par_givens_errs[:len_saved],\
        hier_fgft_errs[:len_saved] = saved_lap_errs
        saved_U_errs = loadtxt('benchmark_pyfaust_fgft_U.txt')
        U_givens_errs[:len_saved], U_par_givens_errs[:len_saved], \
        U_hier_fgft_errs[:len_saved] = \
        saved_U_errs
        saved_fgft_times = loadtxt('benchmark_pyfaust_fgft_times.txt')
        givens_times[:len_saved], par_givens_times[:len_saved],\
                hier_fgft_times[:len_saved], hier_palm_times[:len_saved] = saved_fgft_times
        start_i = len_saved
        if(saved_U_errs[0].shape[0] != saved_lap_errs[0].shape[0]):
            raise Exception("Error files not consistent, delete and recompute")
    for i in range(start_i, num_laps):
        print('Lap ',i,'/',num_laps)
        U,D,Lap = Us[i],Ds[i],Laps[i].astype(numpy.float)
        dim = Lap.shape[0]
        nfacts = int(round(log2(dim))-3)
        over_sp = 1.5 # sparsity overhead
        dec_fact = .5 # decrease of the residum sparsity
        fact_cons, res_cons = [], []
        for j in range(1, nfacts):
            fact_cons += [ ConstraintInt('sp', dim, dim,
                                         min(int(round(dec_fact**j*dim**2*over_sp)),size(Lap)))
                         ]
            res_cons += [ ConstraintInt('sp', dim, dim,
                                        min(int(round(2*dim*over_sp)),size(Lap))) ]

        params = ParamsHierarchicalFact(fact_cons, res_cons,
                                        StoppingCriterion(num_its=50),
                                        StoppingCriterion(num_its=100),
                                        step_size=1.0000e-06,
                                        constant_step_size=False,
                                        init_lambda=1.0,
                                        is_fact_side_left=False)

        t = clock()
        F = FaustFactory.fact_hierarchical(U,params)
        hier_palm_times[i] = clock()-t
        complexity_global = F.nnz_sum()
        rc_palm = complexity_global / count_nonzero(U)

        diag_init_D = diag(D)
        diag_init_D = numpy.copy(diag_init_D)
        t = clock()
        F, Dhat = FaustFactory.fgft_palm(U, Lap,
                                                      params,
                                                      diag_init_D)
        t = clock()-t
        hier_fgft_times[i] = t
        hier_fgft_err = norm((F.todense()*diag(Dhat))*F.T.todense()-Lap,"fro")/norm(Lap,"fro")
        hier_fgft_errs[i] = hier_fgft_err
        U_hier_fgft_errs[i] =  (F-U).norm("fro")/norm(U,"fro")



        #nfacts = complexity_global/(2*dim)
        #J = round(nfacts*(dim/2))
        J = round(complexity_global/4)
        t = clock()
        F, Dhat = FaustFactory.fgft_givens(Lap, J, 0)
        t = clock()-t
        givens_times[i] = t

        givens_err = norm((F*Dhat.todense())*F.T.todense()-Lap,'fro')/norm(Lap,'fro')
        givens_errs[i] = givens_err
        U_givens_errs[i] = (F-U).norm("fro")/norm(U,"fro")

        t = clock()
        F, Dhat = FaustFactory.fgft_givens(Lap, J, int(dim/2))
        t = clock()-t
        par_givens_times[i] = t
        print("J=", J)
        print("nnz_sum FGFT givens parallel=", F.nnz_sum(), "num of facts:",
              F.numfactors())
        par_givens_err = norm((F*Dhat.todense())*F.T.todense()-Lap,'fro')/norm(Lap,'fro')
        par_givens_errs[i] = par_givens_err
        U_par_givens_errs[i] = (F-U).norm("fro")/norm(U,"fro")
        print("err:", norm(Lap @ F.toarray() - F.toarray()@Dhat)/norm(Lap))



        savetxt('benchmark_pyfaust_fgft.txt', np.array([givens_errs[:i+1],
                                                        par_givens_errs[:i+1],
                                                        hier_fgft_errs[:i+1]]))
        savetxt('benchmark_pyfaust_fgft_U.txt', np.array([U_givens_errs[:i+1],
                                                          U_par_givens_errs[:i+1],
                                                          U_hier_fgft_errs[:i+1]]))
        savetxt('benchmark_pyfaust_fgft_times.txt', np.array([givens_times[:i+1],
                                                          par_givens_times[:i+1],
                                                          hier_fgft_times[:i+1],
                                                          hier_palm_times[:i+1]]))


matlab_matfile = None
if(len(argv) > 1 and exists(argv[1]) and re.match(".*.mat", argv[1])):
    matlab_matfile = loadmat(argv[1])

    matlab_hier_fgft_errs = [matlab_matfile['hier_errs'][0,i][0,0] for i in
                             range(matlab_matfile['hier_errs'].shape[1]) ]

    matlab_givens_errs = [matlab_matfile['givens_errs'][0,i][0,0] for i in
                          range(matlab_matfile['givens_errs'].shape[1]) ]


    matlab_par_givens_errs = [matlab_matfile['par_givens_errs'][0,i][0,0] for i in
                              range(matlab_matfile['par_givens_errs'].shape[1]) ]


    matlab_hier_palm_times = [ matlab_matfile['hier_palm_times'][0,i][0,0] for i in
                             range(matlab_matfile['hier_palm_times'].shape[1]) ]

    matlab_hier_fgft_times = [ matlab_matfile['hier_fgft_times'][0,i][0,0] for i in
                             range(matlab_matfile['hier_fgft_times'].shape[1]) ]
    matlab_givens_times = [ matlab_matfile['givens_fgft_times'][0,i][0,0] for i in
                             range(matlab_matfile['givens_fgft_times'].shape[1]) ]
    matlab_par_givens_times = [ matlab_matfile['par_givens_fgft_times'][0,i][0,0] for i in
                             range(matlab_matfile['par_givens_fgft_times'].shape[1]) ]




print('hier_fgft_errs = ', hier_fgft_errs)
print('givens_errs = ', givens_errs)
print('par_givens_errs = ', par_givens_errs)
print('U_hier_fgft_errs = ', U_hier_fgft_errs)
print('U_givens_errs = ', U_givens_errs)
print('U_par_givens_errs = ', U_par_givens_errs)
print('hier_palm_times = ', hier_palm_times)
print('hier_fgft_times = ', hier_fgft_times)
print('givens_times = ', givens_times)
print('par_givens_times = ', par_givens_times)

plt.rcParams['figure.figsize'] = [18.0, 12]
lap_indices = arange(num_laps)

grid(True, which='both', axis='both')
xticks(lap_indices)
scatter(lap_indices, givens_errs, c='b', marker='o', label='pyfaust_givens_errs')
scatter(lap_indices, par_givens_errs, c='r', marker='o', label='pyfaust_par_givens_errs')
scatter(lap_indices, hier_fgft_errs, c='y', marker='o', label='pyfaust_hier_fgft_errs')

if matlab_matfile:
    scatter(lap_indices, matlab_hier_fgft_errs[:len(lap_indices)],
            label='matlab_hier_fgft_errs', c='y', marker='+', s=315)
    scatter(lap_indices, matlab_givens_errs[:len(lap_indices)],
            label='matlab_givens_errs', c='b', marker='+', s=315)
    scatter(lap_indices, matlab_par_givens_errs[:len(lap_indices)],
            label='matlab_par_givens_errs', c='r', marker='+', s=315)

ylabel("Relative error on Laplacian (reconstructed with FGFT)")
xlabel("Laplacians\n i-th Lap. is for"
      " erdos-renyi if i==0 (mod 6)\ncommunity if i==1 (mod 6)\nsensor if i=="
      "2 (mod 6)\npath if i==3 (mod 6)\nrandom ring if i==4 (mod6), ring if i=="
      "5 (mod6)")

title("Error benchmark on "+str(num_laps)+" Laplacians (128x128)")

legend()
savefig('benchmark_Lap_diag_pyfaust_Laplacian_figure.png')

figure(2)
grid(True, which='both', axis='both')

if matlab_matfile:
    matlab_U_hier_fgft_errs = [matlab_matfile['U_hier_errs'][0,i][0,0] for i in
                             range(matlab_matfile['U_hier_errs'].shape[1]) ]

    matlab_U_givens_errs = [matlab_matfile['U_givens_errs'][0,i][0,0] for i in
                          range(matlab_matfile['U_givens_errs'].shape[1]) ]


    matlab_U_par_givens_errs = [matlab_matfile['U_par_givens_errs'][0,i][0,0] for i in
                          range(matlab_matfile['U_par_givens_errs'].shape[1]) ]
    scatter(lap_indices, matlab_U_hier_fgft_errs[:len(lap_indices)],
            label='matlab_hier_fgft_errs', c='y', marker='+', s=315)
    scatter(lap_indices, matlab_U_givens_errs[:len(lap_indices)],
            label='matlab_givens_errs', c='b', marker='+', s=315)
    scatter(lap_indices, matlab_U_par_givens_errs[:len(lap_indices)],
            label='matlab_par_givens_errs', c='r', marker='+', s=315)

xticks(lap_indices)

scatter(lap_indices, U_givens_errs, c='b', marker='o', label='pyfaust_givens_errs')
scatter(lap_indices, U_par_givens_errs, c='r', marker='o', label='pyfaust_par_givens_errs')
scatter(lap_indices, U_hier_fgft_errs, c='y', marker='o', label='pyfaust_hier_fgft_errs')

ylabel("Relative error of FGFT versus U (Fourier matrix)")
xlabel("Fourier matrices")
title("Error benchmark on "+str(num_laps)+" Graph Fourier matrices (128x128)")
legend()

savefig('benchmark_Lap_diag_pyfaust_Fourier_figure.png')
if matlab_matfile:
    print('matlab_hier_fgft_errs = ', matlab_hier_fgft_errs)
    print('matlab_givens_errs = ', matlab_givens_errs)
    print('matlab_par_givens_errs = ', matlab_par_givens_errs)
    print('matlab_U_hier_fgft_errs = ', matlab_U_hier_fgft_errs)
    print('matlab_U_givens_errs = ', matlab_U_givens_errs)
    print('matlab_U_par_givens_errs = ', matlab_U_par_givens_errs)

figure(3)
semilogy(lap_indices, givens_times, c='b', label='pyfaust Givens FGFT Time',
     marker='o')
semilogy(lap_indices, par_givens_times, c='r', label='pyfaust // Givens FGFT Time',
    marker='o')
semilogy(lap_indices, hier_fgft_times, c='y', label='pyfaust PALM FGFT Time',
     marker='o')
semilogy(lap_indices, hier_palm_times, c='black', label='pyfaust PALM U Fac. Time',
     marker='o')

if matlab_matfile:
    plot(lap_indices, matlab_givens_times[:len(lap_indices)], c='b', label='Matlab Givens FGFT Time',
         marker='+', markersize=20)
    plot(lap_indices, matlab_par_givens_times[:len(lap_indices)], c='r', label='Matlab // Givens FGFT Time',
         marker='+', markersize=20)
    plot(lap_indices, matlab_hier_fgft_times[:len(lap_indices)], c='y', label='Matlab PALM FGFT Time',
         marker='+', markersize=20)
    plot(lap_indices, matlab_hier_palm_times[:len(lap_indices)], c='black', label='Matlab PALM U Fac. Time',
         marker='+', markersize=20)

ylabel("Factorization Compute Time (sec)")
xlabel("Laplacians\n i-th Lap. is for"
      " erdos-renyi if i==0 (mod 6)\ncommunity if i==1 (mod 6)\nsensor if i=="
      "2 (mod 6)\npath if i==3 (mod 6)\nrandom ring if i==4 (mod6), ring if i=="
      "5 (mod6)")

xticks(lap_indices)
grid(True, which='both', axis='both')

legend()
savefig('benchmark_pyfaust_fac_time_figure.png')


show()


