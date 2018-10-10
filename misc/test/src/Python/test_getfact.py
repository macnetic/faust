#!/usr/bin/env python3

from pyfaust import *
from time import clock
import numpy
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt

NUM_OF_ITERS = 100


mat_sizes = numpy.linspace(10,1000, 30)
densities = numpy.linspace(0.1, 1, 4)


for is_square in [True]:#,False]:
    for is_transpose in [True,False]:
        for mat_format in ["Sparse", "Dense"]:
            fig, axes = plt.subplots(2,2, sharey=True, sharex=True)
            square_str = "Square"
            if(not is_square): square_str = "Non-"+square_str
            if(is_transpose): square_str += " transpose"
            plt.suptitle("Execution time VS Size of "+mat_format+" "+square_str+\
                         " Matrix", color='blue')

            for i in range(0,len(densities)):
                density = densities[i]
                opt_times = []
                nonopt_times = []
                for j in range(0,len(mat_sizes)):
                    if(is_square):
                        size = int(mat_sizes[j])
                    else:
                        if(j >= 1):
                            size = [mat_sizes[j-1],mat_sizes[j]]
                        else:
                            size = [1, mat_sizes[j]]
                    F = FaustFactory.rand(5, size, fac_type=mat_format.lower(), density=density)
                    if(is_transpose):
                        F = F.T
                    m1 = clock()
                    for k in range(0, NUM_OF_ITERS):
                        F.get_factor(1)
                    m2 = clock()
                    opt_times += [m2-m1]
                    if(mat_format == 'Sparse'):
                        m1 = clock()
                        for k in range(0, NUM_OF_ITERS):
                            csr_matrix(F.get_factor_nonopt(1))
                    else:
                        m1 = clock()
                        for k in range(0, NUM_OF_ITERS):
                            F.get_factor_nonopt(1)
                    m2 = clock()
                    nonopt_times += [m2-m1]


                ax_i, ax_j = int(i/2), int(i%2)
                axes[ax_i,ax_j].set_title(" Matrix"
                                          " density: "+str(density),fontdict={'fontweight':'bold'})
                axes[ax_i, ax_j].plot(mat_sizes, opt_times, c='blue',
                                      label='get_factor_opt' )
                axes[ax_i, ax_j].plot(mat_sizes, nonopt_times, c='red',
                                      label='get_factor_nonopt')
                if(ax_i == 1):
                    axes[ax_i,ax_j].set_xlabel("Matrix Size")
                if(ax_j == 0):
                    axes[ax_i,ax_j].set_ylabel("Time (s) (for "
                                               +str(NUM_OF_ITERS)+" calls)")
                axes[ax_i, ax_j].legend()
            #plt.show()
            figname = \
                   ("get_factor_exec_time-"+mat_format+"_"+square_str+".png").replace(" ","_")
            fig.savefig(figname,dpi=200)
            print("Figure saved: "+figname)

plt.show()
