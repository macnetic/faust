
from __future__ import print_function
from pylab import *
import os,sys



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
        from time import time, clock
        if sys.platform == 'win32':
            timer = clock
        else:
            timer = time

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

        # time comparison
        from time import time, clock
        if sys.platform == 'win32':
            timer = clock
        else:
            timer = time

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


