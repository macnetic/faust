
from pylab import *


def fig_speed_up_fourier():
    from pyfaust.demo.fft.speed_up_fourier import nb_mults, log2_dims, dims, FFT_FAUST, \
    FFT_FAUST_FULL, FFT_NATIVE, NUM_FFT_TYPES
    import os.path, os, sys

    output_dir = 'pyfaust_demo_figures'

    input_dir = 'pyfaust_demo_output'
    input_path = input_dir+os.sep+'speed_up_fourier.txt'
    mult_times = loadtxt(input_path).reshape(nb_mults, len(dims), NUM_FFT_TYPES)
    # NUM_FFT_TYPES == 3 is for FFT_FAUST, FFT_FAUST_FULL, FFT_NATIVE

    mean_mult_times = squeeze(mult_times.mean(axis=0))

    curve_thickness = 2
    legend_loc = 'upper left'

    plt.rcParams['figure.figsize'] = [12.0, 8]
    subplot(211)

    hold(True)
    line_marker_types = [ 'ro-', 'bo-', 'go-', 'r+-', 'b+-', 'g+-' ]

    title('Runtime Fourier A*x')
    for i in range(0,NUM_FFT_TYPES):
        semilogy(log2_dims, mean_mult_times[:,i],line_marker_types[i],
                 lw=curve_thickness)
        ylabel('Computed Time (sec)')
    # legend(['dense', 'Faust', 'fft'], loc=legend_loc) # in first subplot
    grid(True)
    axes([log2_dims[0], log2_dims[-1], mean_mult_times.min(),
          mean_mult_times.max()])
    xticks(log2_dims)

    subplot(212)
    grid(True)
    title('Speedup Fourier A*x')
    for i in range(1,NUM_FFT_TYPES):
        semilogy(log2_dims, mean_mult_times[:,0]/mean_mult_times[:,i],
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





