from pylab import *
from pyfaust.demo.runtimecmp.runtime_comparison import rcgs, dims, nb_facts, nb_mults,\
        dims_len, rcgs_len, nb_facts_len

def fig_runtime_comparison():
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

    tfaust = loadtxt(path_tfaust).reshape(nb_mults, dims_len, rcgs_len, nb_facts_len, 2)
    tdense = loadtxt(path_tdense).reshape(nb_mults, dims_len, 2)


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

    fig, ax = subplots(2,nb_facts_len, sharex=True, sharey=True)#, h*nb_facts_len+nf+1)
    for h in arange(0,2):
        for nf in range(0,nb_facts_len):

            legend_curve = []
            lines = []
            for k in range(0, rcgs_len):
                lines.append(*ax[h,nf].semilogy(log2(dims), mean_tfaust[:, k, nf, h],
                         '-+', lw=curve_thickness))
                legend_curve.append('Faust RCG '+str(rcgs[k]))
                hold(True)

            lines.append(*ax[h,nf].semilogy(log2(dims), squeeze(mean_tdense[:,h]), '-+', c=(0, .8, .8),
                     lw=curve_thickness))

            legend_curve.append('Dense ')

            ax[h,nf].grid(True)
            axes([log2(dims[0]), log2(dims[-1]), ymin, ymax])
            if(h == 0):
                ax[h,nf].set_title('#factors: '+str(nb_facts[nf]))
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
    #        '(constraint: 'constraint_str+')')
    savefig(output_dir+os.sep+'RuntimeComp-'+matrix_or_vector+'_multiplication_constraint_'+constraint_str+'.png',
           dpi=200)
    #tight_layout()
    #show()
