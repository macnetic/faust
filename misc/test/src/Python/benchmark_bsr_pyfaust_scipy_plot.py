from pylab import *

"""
This script goal is to plot the data produced by
benchmark_bsr_pyfaust_scipy_plot.py script.
"""

names = ['pyfaust BSR Faust', 'pyfaust CSR Faust', 'Array Faust', 'scipy BSR Faust']

fig, ax = subplots(1,3)#, sharey=True)

for a in ax:
    #a.set_yscale('log')
    a.set_xlabel('Faust type')
    a.set_xticklabels(names[0:2]+names[3:], rotation=90)
mul_times = loadtxt('mul_dense.txt')
bsr_csr_speedup = median(mul_times[1,:]/mul_times[0,:])
bsr_scipy_speedup = median(mul_times[3,:]/mul_times[0,:])
ax[0].set_xlabel("\nBSR Faust median speedups:\n- vs CSR:%.2f  " % bsr_csr_speedup + '\n'
+ "- vs scipy: %.2f" % bsr_scipy_speedup, fontsize=9)
mul_times = vstack((mul_times[0:2], mul_times[3:]))
ax[0].boxplot(mul_times.T)
ax[0].set_title('Faust-dense matrix\nproduct', fontsize=10)
ax[0].set_ylabel('time (s)')
mul_times = loadtxt('mul_csr.txt')
bsr_csr_speedup = median(mul_times[1,:]/mul_times[0,:])
bsr_scipy_speedup = median(mul_times[3,:]/mul_times[0,:])
ax[1].set_xlabel("\nBSR Faust median speedups:\n- vs CSR:%.2f  " % bsr_csr_speedup + '\n'
+ "- vs scipy: %.2f" % bsr_scipy_speedup, fontsize=9)
mul_times = vstack((mul_times[0:2], mul_times[3:]))
ax[1].boxplot(mul_times.T)
ax[1].set_title('Faust-CSR sparse matrix\nproduct', fontsize=10)
mul_times = loadtxt('mul_vec.txt')
bsr_csr_speedup = median(mul_times[1,:]/mul_times[0,:])
bsr_scipy_speedup = median(mul_times[3,:]/mul_times[0,:])
ax[2].set_xlabel("\nBSR Faust median speedups:\n- vs CSR:%.2f  " % bsr_csr_speedup + '\n'
+ "- vs scipy: %.2f" % bsr_scipy_speedup, fontsize=9)
mul_times = vstack((mul_times[0:2], mul_times[3:]))
ax[2].boxplot(mul_times.T)
ax[2].set_title('Faust-vector product', fontsize=10)
suptitle("BSR Faust Benchmark ("+str(mul_times[0].size)+" runs)", fontweight="bold")
tight_layout()
savefig('mul_times_box_plot.png')
