
from __future__ import print_function
from pylab import *
import os,sys


DEFT_DATA_DIR = 'pyfaust_demo_output'
DEFT_FIG_DIR = 'pyfaust_demo_figures'

class quickstart:

    @staticmethod
    def run():
        #import the pyfaust package
        import pyfaust;

        # import module to generate data
        import scipy
        from scipy import sparse as sp
        import numpy as np


        # generate the factors of the Faust
        dim1 = 1000;
        dim2 = 2000;
        nb_factor = 2;
        list_factor_sparse=[0]*nb_factor
        list_factor=[0]*nb_factor
        int_max = 120
        density_per_fact=0.1;
        list_factor_sparse[0]=int_max*sp.random(dim1,dim1,density=density_per_fact,format='csr',dtype=np.float64);
        list_factor[0]=list_factor_sparse[0].toarray();
        list_factor_sparse[1]=int_max*sp.random(dim1,dim2,density=density_per_fact,format='csr',dtype=np.float64);
        list_factor[1]=list_factor_sparse[1].toarray();


        #print(list_factor[0])
        #print(list_factor[1])


        # create a Faust named F from its factors
        A = pyfaust.Faust(list_factor)

        # get the size of the Faust
        print("dimension of the Faust : ", A.shape)

        # transpose a Faust
        A_trans = A.transpose()

        # multiplication a Numpy matrix by a Faust
        x = np.random.randint(int_max, size=(dim2,1))
        y = A * x

        # convert a faust to numpy matrix
        A_numpy = A.toarray()

        # slicing
        coeff = A[0,0]
        col_2nd = A[:,1];
        submatrix_A = A[3:5,2:3]



        # speed-up multiplication
        import time
        nb_mult = 100
        t_dense = 0.0;
        t_faust = 0.0;

        for i in range(nb_mult):
            x=np.random.randint(int_max, size=(dim2,1))

            t_begin = time.time()
            y_dense = A_numpy.dot(x)
            t_elapsed = time.time() - t_begin
            t_dense += t_elapsed

            t_begin = time.time()
            y_faust=A*x;
            t_elapsed = time.time()-t_begin
            t_faust += t_elapsed

        print("multiplication SPEED-UP using Faust")
        print("Faust is "+str(t_dense/t_faust)+" faster than a full matrix")
        print("Faust nnz: "+str(A.nnz_sum()))
        print("Faust density: "+str(A.density()))
        print("Faust RCG: "+str(A.rcg()))
        print("Faust norm: "+str(A.norm()))
        print("Faust nb of factors: "+str(A.get_num_factors()))
        for i in range(0,A.get_num_factors()):
            #print("Faust size of factor ",i,"=",A.get_factor(i).shape)
            # test Faust gets back the same sparse factors given at init
            assert((A.get_factor(i).toarray() == list_factor_sparse[i]).all())
            #print(A.get_factor(i))

        # test Faust saving
        A.save("A.mat")
        As = pyfaust.Faust(filepath="A.mat")
        assert((A.get_factor(0).toarray() == As.get_factor(0).toarray()).all())
        assert((A.get_factor(1).toarray() == As.get_factor(1).toarray()).all())

        # test Faust transpose
        #print(A.get_factor(0))
        tA = A.transpose()
        tf1 = tA.get_factor(1).toarray()
        #print(tf1)
        f1 = np.transpose(tf1)
        assert(not (tf1 == A.get_factor(0).toarray()).all() or (tf1 == f1).all())
        assert((f1 == A.get_factor(0).toarray()).all())

        print("end quickstart.py")

class fft:

    _nb_mults = 500
    _log2_dims = arange(6,13)
    _dims = 2**_log2_dims

    _FFT_FAUST_FULL=0
    _FFT_FAUST=1
    _FFT_NATIVE=2

    _NUM_FFT_TYPES=3


    @staticmethod
    def fig_speed_up_fourier():
        """
        TODO
        """
        import os.path

        output_dir = 'pyfaust_demo_figures'

        input_dir = 'pyfaust_demo_output'
        input_path = input_dir+os.sep+'speed_up_fourier.txt'
        mult_times = loadtxt(input_path).reshape(fft._nb_mults, len(fft._dims), fft._NUM_FFT_TYPES)
        # fft._NUM_FFT_TYPES == 3 is for fft._FFT_FAUST, fft._FFT_FAUST_FULL, fft._FFT_NATIVE

        mean_mult_times = squeeze(mult_times.mean(axis=0))

        curve_thickness = 2
        legend_loc = 'upper left'

        plt.rcParams['figure.figsize'] = [12.0, 8]
        subplot(211)

        hold(True)
        line_marker_types = [ 'ro-', 'bo-', 'go-', 'r+-', 'b+-', 'g+-' ]

        title('Runtime Fourier A*x')
        for i in range(0,fft._NUM_FFT_TYPES):
            semilogy(fft._log2_dims, mean_mult_times[:,i],line_marker_types[i],
                     lw=curve_thickness)
            ylabel('Computed Time (sec)')
        # legend(['dense', 'Faust', 'fft'], loc=legend_loc) # in first subplot
        grid(True)
        axes([fft._log2_dims[0], fft._log2_dims[-1], mean_mult_times.min(),
              mean_mult_times.max()])
        xticks(fft._log2_dims)

        subplot(212)
        grid(True)
        title('Speedup Fourier A*x')
        for i in range(1,fft._NUM_FFT_TYPES):
            semilogy(fft._log2_dims, mean_mult_times[:,0]/mean_mult_times[:,i],
                     line_marker_types[i], lw=curve_thickness)
            ylabel('Speedup Fourier A*x')
            xlabel ('log(dim)')

        figlegend(figure(1).get_axes()[0].get_lines(),['dense', 'faust',
                                                       'fft'],loc='best')


        if(not os.path.exists(output_dir)):
            os.mkdir(output_dir)
        savefig(output_dir+os.sep+'Fourier-RuntimeComp-Speedup.png',
                dpi=200)
        show()



    @staticmethod
    def speed_up_fourier():
        from pyfaust import Faust, FaustFactory
        treshold = 10**-10

        print('Speedup Fourier')
        print("===============")
        print('Generating data...')

        fft_mats = []
        fft_fausts = []

        for k in range(0,len(fft._dims)):
            print('\rFFT fft._dims processed: ', fft._dims[0:k+1], end='')
            F = FaustFactory.fourier(fft._log2_dims[k])
            fft_mat = F.toarray()#fft(eye(fft._dims[k])) # or TODO: F.toarray() ?
            fft_fausts += [ F ]
            fft_mats += [ fft_mat ]
        print()

        print("Gathering computation times...")

        mult_times = ndarray(shape=(fft._nb_mults, len(fft._dims), 3))
        # 3 for: fft._FFT_FAUST_FULL, fft._FFT_FAUST, fft._FFT_NATIVE


        for i in range(0,fft._nb_mults):
            print('\r#muls:',i+1,'/', fft._nb_mults, end='')
            for k in range(0,len(fft._dims)):
               dim = fft._dims[k]
               fft_mat = fft_mats[k]
               fft_faust = fft_fausts[k]
               ydense = empty(dim)
               yfaust = empty(dim)
               yfft = empty(dim)
               ydense_trans = empty(dim)
               yfaust_trans = empty(dim)
               x = rand(dim,1)
               t = timer()
               ydense = fft_mat.dot(x)
               mult_times[i,k,fft._FFT_FAUST_FULL] = timer()-t
               t = timer()
               yfaust = fft_faust*x
               mult_times[i,k,fft._FFT_FAUST] = timer()-t
               t = timer()
               yfft = fft2(x)
               mult_times[i,k,fft._FFT_NATIVE] = timer()-t
               if(norm(ydense-yfaust)>treshold):
                   raise Exception('Multiplication error: larger error than '
                                   'treshold for ydense.')

               from numpy.fft import fft as npfft
               n1 = norm(npfft(np.eye(fft._dims[k]))-fft_faust.toarray())
               assert(n1 < treshold)
               if(norm(yfft-yfaust)/norm(yfft)>treshold):
                   print('\nerror:', norm(yfft-yfaust))
                   raise Exception('Multiplication error: larger error than '
                                   'treshold for yfaust')

        print()


        output_dir = 'pyfaust_demo_output'
        import os.path
        if(not os.path.exists(output_dir)):
           os.mkdir(output_dir)
        output_path = output_dir+os.sep+'speed_up_fourier.txt'
        savetxt(output_path, mult_times.reshape(mult_times.size))

        # test
        # mult_times_r = loadtxt(output_path).reshape(mult_times.shape)
        # assert(all(mult_times_r == mult_times))



class runtimecmp:

    _rcgs = [2, 4, 8]
    _dims = [128, 256, 512]
    _nb_facts = [2, 4, 8]
    _nb_mults = 500
    _constraint = True # 'per_row' # per_col
    _dims_len = len(_dims)
    _rcgs_len = len(_rcgs)
    _nb_facts_len = len(_nb_facts)

    @staticmethod
    def runtime_comparison():
        """
        TODO
        """
        from pyfaust import Faust, FaustFactory

        matrix_or_vector = 'vector'


        fausts = ndarray(shape=(runtimecmp._dims_len, runtimecmp._rcgs_len, runtimecmp._nb_facts_len), dtype=Faust)
        dense_mats =  list()


        # loading all the different Fausts and dense matrices
        for j in range(0,runtimecmp._dims_len):
            dim = runtimecmp._dims[j]
            A = rand(dim, dim)
            dense_mats.append(A)

            for k in range(0,runtimecmp._rcgs_len):
                rcg = runtimecmp._rcgs[k]

                for l in range(0,runtimecmp._nb_facts_len):
                    nf = runtimecmp._nb_facts[l]
                    F = FaustFactory.rand(nf , dim, 1./(nf*rcg),
                                          per_row=runtimecmp._constraint, fac_type='sparse')
                    assert(F.rcg() == rcg)
                    fausts[j,k,l] = F
                    assert(F.shape == (dim,dim))



        tdense = ndarray(shape=(runtimecmp._nb_mults, runtimecmp._dims_len, 2))
        tfaust = ndarray(shape=(runtimecmp._nb_mults, runtimecmp._dims_len,
                                runtimecmp._rcgs_len, runtimecmp._nb_facts_len, 2))

        for i in range(0, runtimecmp._nb_mults):
            print("\r\r #muls =",i+1,'/',runtimecmp._nb_mults, end='')
            for j in range(0,runtimecmp._dims_len):
                dim = runtimecmp._dims[j]

                if(matrix_or_vector == 'matrix'):
                    dim2 = dim # mul. by a square matrix
                elif(matrix_or_vector == 'vector'):
                    dim2 = 1
                else:
                    raise("matrix_or_vector string must be equal to matrix or"
                          " vector")

                for k in range(0,runtimecmp._rcgs_len):
                    rcg = runtimecmp._rcgs[k]

                    for l in range(0,runtimecmp._nb_facts_len):
                        nf = runtimecmp._nb_facts[l]
                        x = rand(dim,dim2)

                        if(k == 0 and l == 0):
                            A = dense_mats[j]
                            t = timer()
                            y = A.dot(x)
                            tdense[i,j,0] = timer()-t
                            t = timer()
                            y_trans = A.T.dot(x)
                            tdense[i,j,1] = timer()-t

                        F = fausts[j,k,l]
                        t = timer()
                        yfaust = F*x
                        tfaust[i,j,k,l,0] = timer()-t
                        t = timer()
                        yfaust_trans = F.T*x
                        tfaust[i,j,k,l,1] = timer()-t

        print()
        output_dir = 'pyfaust_demo_output'
        import os.path
        if(not os.path.exists(output_dir)):
           os.mkdir(output_dir)
        path_tfaust = output_dir+os.sep+'runtime_cmp_tfaust-'+matrix_or_vector+'.txt'
        path_tdense = output_dir+os.sep+'runtime_cmp_tdense-'+matrix_or_vector+'.txt'
        savetxt(path_tdense, tdense.reshape(tdense.size))
        savetxt(path_tfaust, tfaust.reshape(tfaust.size))

        # test
        #tfaust_r = loadtxt(path_tfaust).reshape(tfaust.shape)
        #assert(all(tfaust_r == tfaust))
        #tdense_r = loadtxt(path_tdense).reshape(tdense.shape)
        #assert(all(tdense_r == tdense))

    @staticmethod
    def fig_runtime_comparison():
        """
        TODO
        """
        import os.path, os
        output_dir = 'pyfaust_demo_figures'
        input_dir = 'pyfaust_demo_output'
        input_file_existing = False

        for matrix_or_vector in  ('vector', 'matrix'):
            path_tfaust = input_dir+os.sep+'runtime_cmp_tfaust-'+matrix_or_vector+'.txt'
            path_tdense = input_dir+os.sep+'runtime_cmp_tdense-'+matrix_or_vector+'.txt'
            if(os.path.exists(path_tfaust) and os.path.exists(path_tdense)):
                input_file_existing = True
                break

        if(not input_file_existing):
            raise Exception("Input files don't exist, please run"
                            " runtime_comparison.py first.")

        if(not os.path.exists(output_dir)):
            os.mkdir(output_dir)

        tfaust = loadtxt(path_tfaust).reshape(runtimecmp._nb_mults,
                                              runtimecmp._dims_len, runtimecmp._rcgs_len, runtimecmp._nb_facts_len, 2)
        tdense = loadtxt(path_tdense).reshape(runtimecmp._nb_mults, runtimecmp._dims_len, 2)


        # average the mul. times according the number of runs
        mean_tdense = squeeze(tdense.mean(axis=0))
        mean_tfaust = squeeze(tfaust.mean(axis=0))

        # avoid legend overlapping axis
        plt.rcParams['figure.figsize'] = [12.0, 8]
        #print("mean_tdense=", mean_tdense, mean_tdense.shape)

        # plot the time computed in logscale with a fixed number of factors in a given figure
        # in each figure we have in x-axis the log2(dimension) of the square matrix
        #                         in y-axis the time
        # all times for faust multiplication with different RCG (density)
        # and the times for dense matrix multiplication are plotted

        curve_thickness = 2

        # hold differents figures in the same box
        ymin = min(mean_tdense.min(), mean_tfaust.min())
        ymax = max(mean_tdense.max(), mean_tfaust.max())

        legendy = [ 'Time (A*x)', 'Time (A.T*x)']

        fig, ax = subplots(2,runtimecmp._nb_facts_len, sharex=True, sharey=True)#, h*runtimecmp._nb_facts_len+nf+1)
        for h in arange(0,2):
            for nf in range(0,runtimecmp._nb_facts_len):

                legend_curve = []
                lines = []
                for k in range(0, runtimecmp._rcgs_len):
                    lines.append(*ax[h,nf].semilogy(log2(runtimecmp._dims), mean_tfaust[:, k, nf, h],
                             '-+', lw=curve_thickness))
                    legend_curve.append('Faust RCG '+str(runtimecmp._rcgs[k]))
                    hold(True)

                lines.append(*ax[h,nf].semilogy(log2(runtimecmp._dims), squeeze(mean_tdense[:,h]), '-+', c=(0, .8, .8),
                         lw=curve_thickness))

                legend_curve.append('Dense ')

                ax[h,nf].grid(True)
                axes([log2(runtimecmp._dims[0]), log2(runtimecmp._dims[-1]), ymin, ymax])
                if(h == 0):
                    ax[h,nf].set_title('#factors: '+str(runtimecmp._nb_facts[nf]))
                # legend for it axis (a bit heavy to read)
                #ax[h,nf].legend(legend_curve)
                #legend('best')
                if(nf == 0):
                    if(h == 1):
                        ax[h,nf].set_xlabel("log2(Dimension)")
                    ax[h,nf].set_ylabel(legendy[h])

        figlegend(lines,legend_curve,loc='upper left')
        # TODO: figure.Name matlab ?
        constraint_str = 'sp'
        # fig_name = 'Faust-'+matrix_or_vector+' multiplication '
        #        '(_constraint: 'constraint_str+')')
        savefig(output_dir+os.sep+'RuntimeComp-'+matrix_or_vector+'_multiplication_constraint_'+constraint_str+'.png',
               dpi=200)
        #tight_layout()
        #show()

class hadamardfact:

    _nb_mults = 500
    _n = 5
    _nfacts = _n
    _tdense_fname = 'hadamardfact_mul_runtime_dense.txt'
    _tfaust_fname = 'hadamardfact_mul_runtime_faust.txt'
    _had_faust_fname = 'hadamardfact_faust.mat'
    _tspeedup_times_fname = 'hadamardfact_speedup_times.txt'
    _fig1_fname = 'Hadamard-factorization.png'
    _fig2_fname = 'Hadamard-factorization_nnz_coeff.png'
    _fig_speedup = 'Hadamard-speedup.png'
    _HAD_DENSE, _HAD_TRANS_DENSE, _HAD_FAUST, _HAD_TRANS_FAUST = range(0,4)
    _NUM_TIME_TYPES = 4
    _log2_dims = arange(6,13)
    _dims = 2**_log2_dims

    @staticmethod
    def run_fact(output_dir=DEFT_DATA_DIR):
        from pyfaust import FaustFactory
        from pyfaust.factparams import ParamsHierarchicalFact, ConstraintInt, \
        ConstraintName, StoppingCriterion

        # generate a Hadamard transform and factorize its full matrix
        n = hadamardfact._n
        d = 2**n
        H = FaustFactory.hadamard(n)
        full_H = H.toarray()

        params = ParamsHierarchicalFact(hadamardfact._nfacts, is_update_way_R2L=True, init_lambda=1.0,
                                        fact_constraints=[ConstraintInt(ConstraintName(ConstraintName.SPLINCOL),d,d,2)
                                        for i in range(0,n-1)],
                                        res_constraints=[ConstraintInt(ConstraintName(ConstraintName.SPLINCOL),d,d,int(d/2.**(i+1)))
                                         for i in range(0,n-1)],
                                        data_num_rows=d, data_num_cols=d,
                                        stop_crits=[StoppingCriterion(num_its=30),StoppingCriterion(num_its=30)],
                                        is_verbose=False)
        had_faust = FaustFactory.fact_hierarchical(H.toarray(), params)
        full_had_faust = had_faust.toarray()
        rel_err = norm(full_had_faust-full_H)/norm(full_H)
        print("\n\nRelative error between hadamard matrix and its transform: ",
              rel_err)

        # gather computation times
        _nb_mults = hadamardfact._nb_mults
        dense_times = empty(_nb_mults)
        faust_times = empty(_nb_mults)

        for i in range(0,_nb_mults):
            print("\r\r #muls =",i+1,'/', _nb_mults, end='')
            x = rand(d,1)

            t = timer()
            y_X = full_H.dot(x)
            dense_times[i] = timer()-t

            t = timer()
            y_faust = had_faust*x
            faust_times[i] = timer()-t



        print()


        path_tdense = _write_array_in_file(output_dir,
                                           hadamardfact._tdense_fname,
                                           dense_times)
        path_tfaust = _write_array_in_file(output_dir,
                                           hadamardfact._tfaust_fname,
                                           faust_times)
        had_faust.save(_prefix_fname_with_dir(output_dir,
                                              hadamardfact._had_faust_fname))

        # test
#        tfaust_r = loadtxt(path_tfaust)
#        assert(all(tfaust_r == faust_times))
#        tdense_r = loadtxt(path_tdense)
#        assert(all(tdense_r == dense_times))


    @staticmethod
    def fig_fact(input_dir=DEFT_DATA_DIR, output_dir=DEFT_FIG_DIR):
        from pyfaust import Faust
        fig_dir = DEFT_FIG_DIR
        if(not os.path.exists(fig_dir)):
            os.mkdir(fig_dir)
        hold(True)

        had_faust = Faust(filepath=_prefix_fname_with_dir(input_dir,
                                                          hadamardfact._had_faust_fname))
        fig1 = figure(1)
        subplot("1"+str(had_faust.get_num_factors()+1)+'1')
        imshow(had_faust.toarray())
        xticks([])
        yticks([])
        facts = [];
        for i in range(0,had_faust.get_num_factors()):
            subplot("1"+str(hadamardfact._nfacts+1)+str(i+2))
            # all factors are normally sparse
            fac = had_faust.get_factor(i)
            facts.append(fac)
            if(not isinstance(fac,ndarray)):
                fac = fac.toarray()
            imshow(fac)
            xticks([])
            yticks([])

        fig2 = figure(2)
        subplot("1"+str(had_faust.get_num_factors()+1)+'1')
        imshow(had_faust.toarray())
        xticks([])
        yticks([])
        for i in range(0,had_faust.get_num_factors()):
            subplot("1"+str(hadamardfact._nfacts+1)+str(i+2))
            title("nz = "+str(count_nonzero(fac)))
            spy(facts[i], markersize=1)
            xticks([])
            yticks([])

        _write_fig_in_file(output_dir, hadamardfact._fig1_fname, fig1)
        _write_fig_in_file(output_dir, hadamardfact._fig2_fname, fig2)
        show()

    @staticmethod
    def  run_speedup_hadamard(output_dir=DEFT_DATA_DIR):
       from pyfaust import FaustFactory
       treshold = 10.**-10
       print("Speedup Hadamard")
       print("================")
       print("Generating data...")

       had_mats = []
       had_fausts = []
       _dims, _log2_dims = hadamardfact._dims, hadamardfact._log2_dims
       for k in range(0, len(_dims)):
           print("\rHadamard dims processed: ", _dims[0:k+1], end='')
           F = FaustFactory.hadamard(_log2_dims[k])
           had_fausts += [ F ]
           had_mats += [ F.toarray() ]
       print()

       print("Gathering multiplication computation times...")

       _nb_mults = hadamardfact._nb_mults
       _NUM_TIME_TYPES = hadamardfact._NUM_TIME_TYPES
       _HAD_DENSE, _HAD_TRANS_DENSE, _HAD_FAUST, _HAD_TRANS_FAUST = \
       hadamardfact._HAD_DENSE, hadamardfact._HAD_TRANS_DENSE, \
       hadamardfact._HAD_FAUST, hadamardfact._HAD_TRANS_FAUST
       mult_times = ndarray(shape=(_nb_mults, len(_dims), _NUM_TIME_TYPES))

       for i in range(0,_nb_mults):
           print('\r#muls:',i+1,'/', _nb_mults, end='')
           for k in range(0,len(_dims)):
               dim = _dims[k]
               had_mat = had_mats[k]
               had_faust = had_fausts[k]

               x = rand(dim,1)

               t = timer()
               ydense = had_mat.dot(x)
               mult_times[i,k,_HAD_DENSE] = timer()-t

               t = timer()
               yfaust = had_faust*x
               mult_times[i,k,_HAD_FAUST] = timer()-t

               if(norm(ydense-yfaust) > treshold):
                   raise Exception("speedup hadamard: mul. error greater than "
                                   "treshold")

               t = timer()
               ydense_trans = had_mat.T.dot(x)
               mult_times[i,k,_HAD_TRANS_DENSE] = timer()-t

               t = timer()
               yfaust_trans = had_faust.T*x
               mult_times[i,k,_HAD_TRANS_FAUST] = timer()-t

               if(norm(yfaust_trans-ydense_trans) > treshold):
                   raise Exception("speedup_hadamard: mul. error on transpose "
                                   "faust/mat.")


       print()

       path_mult_times = _write_array_in_file(output_dir,
                                              hadamardfact._tspeedup_times_fname,
                                              mult_times.reshape(mult_times.size))
       mult_times_r = loadtxt(path_mult_times)
       assert(all(mult_times_r.reshape(mult_times.shape) == mult_times))

    @staticmethod
    def fig_speedup_hadamard(output_dir=DEFT_FIG_DIR, input_dir=DEFT_DATA_DIR):
        times_txt_fpath = _prefix_fname_with_dir(input_dir, hadamardfact._tspeedup_times_fname)
        if(not os.path.exists(times_txt_fpath)):
            raise Exception("Input file doesn't exist, please call "
                            "run_speedup_hadamard() before calling "
                            "fig_speedup_hadamard()")
        _nb_mults, _dims, _NUM_TIME_TYPES = hadamardfact._nb_mults, \
                hadamardfact._dims, hadamardfact._NUM_TIME_TYPES
        mult_times = loadtxt(times_txt_fpath).reshape(_nb_mults, len(_dims),
                                                      _NUM_TIME_TYPES)

        _HAD_DENSE, _HAD_TRANS_DENSE, _HAD_FAUST, _HAD_TRANS_FAUST = \
               hadamardfact._HAD_DENSE, hadamardfact._HAD_TRANS_DENSE, \
               hadamardfact._HAD_FAUST, hadamardfact._HAD_TRANS_FAUST
        _log2_dims = hadamardfact._log2_dims

        # mean times of each category and speedups
        mean_mult_times = mult_times.mean(axis=0)
        print(mult_times)
        speedup_trans = mean_mult_times[:,_HAD_TRANS_DENSE] / \
        mean_mult_times[:,_HAD_TRANS_FAUST]
        speedup = mean_mult_times[:,_HAD_DENSE] / mean_mult_times[:,_HAD_FAUST]
        plt.rcParams['figure.figsize'] = [12.0, 8]

        # plot results
        line_width = 2.
        for t in [("1","A*x", _HAD_FAUST, _HAD_DENSE),("3", "A.T*x",
                                                       _HAD_TRANS_FAUST,
                                                       _HAD_TRANS_DENSE)]:
            ymin = min(min(mean_mult_times[:,t[2]]),
                       min(mean_mult_times[:,t[3]]))
            ymax = max(max(mean_mult_times[:,t[2]]),
                       max(mean_mult_times[:,t[3]]))

            subplot("22"+t[0])
            title('Runtime Hadamard '+t[1])
            grid(True)
            hold(True)
            semilogy(_log2_dims, mean_mult_times[:,t[2]], lw=line_width)
            semilogy(_log2_dims, mean_mult_times[:, t[3]], lw=line_width)
            ylabel('Computed Time (sec)')
            if(t[0] == "3"): xlabel('log(dim)')
            legend(['Faust', 'Dense'])
            axes([_log2_dims[0], _log2_dims[-1], ymin, ymax])

        for t in [("2","A*x", speedup),("4", "A.T*x", speedup_trans)]:
            ymin = min(t[2])
            ymax = max(t[2])

            subplot("22"+t[0])
            title('Speedup Hadamard '+t[1])
            grid(True)
            hold(True)
            semilogy(_log2_dims, t[2], lw=line_width)
            semilogy(_log2_dims, t[2]/t[2], lw=line_width, c='black')
            ylabel('Speedup')
            if(t[0] == "4"): xlabel('log(dim)')
            legend(['Faust', 'Neutral'])
            axes([_log2_dims[0], _log2_dims[-1], ymin, ymax])

        _write_fig_in_file(output_dir, hadamardfact._fig_speedup, figure(1))
        show()

def _write_array_in_file(output_dir, fname, array):
    """
    output_dir: directory folder created if doesn't exist.
    fname: output file in output_dir.
    array: the numpy ndarray is flattened before saving. (N-D to 1D). You need
    to keep the dimensions to read-reshape to original array.
    """
    _create_dir_if_doesnt_exist(output_dir)
    fpath = _prefix_fname_with_dir(output_dir, fname)
    savetxt(fpath, array.reshape(array.size))
    return fpath

def _write_fig_in_file(output_dir, fname, fig, dpi=200):
    """
    output_path: directory folder created if doesn't exist.
    fname: output file in output_dir.
    fig: the figure to save.
    """
    _create_dir_if_doesnt_exist(output_dir)
    fpath = _prefix_fname_with_dir(output_dir, fname)
    if(not isinstance(fig, Figure)):
        raise Exception("fig must be a Figure object")
    fig.savefig(fpath,dpi=dpi)

def _prefix_fname_with_dir(dir, fname):
    import os.path
    return dir+os.sep+fname

def _create_dir_if_doesnt_exist(output_dir):
    import os.path
    if(not os.path.exists(output_dir)):
        os.mkdir(output_dir)

# time comparison function to use
from time import time, clock
if sys.platform == 'win32':
    timer = clock
else:
    timer = time



