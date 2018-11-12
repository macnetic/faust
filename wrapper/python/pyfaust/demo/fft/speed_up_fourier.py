from __future__ import print_function
from pylab import *
import os,sys

nb_mults = 500
log2_dims = arange(6,12)
dims = 2**log2_dims

FFT_FAUST_FULL=0
FFT_FAUST=1
FFT_NATIVE=2

NUM_FFT_TYPES=3

def speed_up_fourier():
    from pyfaust import Faust, FaustFactory
    treshold = 10**-10

    print('Speedup Fourier')
    print("===============")
    print('Generating data...')

    fft_mats = []
    fft_fausts = []

    for k in range(0,len(dims)):
        print('\rFFT dims processed: ', dims[0:k+1], end='')
        F = FaustFactory.fourier(log2_dims[k])
        fft_mat = F.toarray()#fft(eye(dims[k])) # or TODO: F.toarray() ?
        fft_fausts += [ F ]
        fft_mats += [ fft_mat ]
    print()

    print("Gathering computation times...")
    from time import time, clock
    if sys.platform == 'win32':
        timer = clock
    else:
        timer = time

    mult_times = ndarray(shape=(nb_mults, len(dims), 3))
    # 3 for: FFT_FAUST_FULL, FFT_FAUST, FFT_NATIVE


    for i in range(0,nb_mults):
        print('\r#muls:',i+1,'/', nb_mults, end='')
        for k in range(0,len(dims)):
           dim = dims[k]
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
           mult_times[i,k,FFT_FAUST_FULL] = timer()-t
           t = timer()
           yfaust = fft_faust*x
           mult_times[i,k,FFT_FAUST] = timer()-t
           t = timer()
           yfft = fft2(x)
           mult_times[i,k,FFT_NATIVE] = timer()-t
           if(norm(ydense-yfaust)>treshold):
               raise Exception('Multiplication error: larger error than '
                               'treshold for ydense.')

           n1 = norm(fft(np.eye(dims[k]))-fft_faust.toarray())
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



