from FaustPy import RandFaustType, Faust
import numpy as np
import os
import matplotlib.pyplot as plt

if __name__ == '__main__':
    nfactors = 3
    startd = 0.01
    endd = 1
    min_dim_sz = 1000
    max_dim_sz = 1000
    sizes = []
    rcgs = []
    ntests = 100
    for i, d in zip(
                    list(range(0, ntests)),
        np.linspace(startd, endd, ntests)):
        F = Faust.randFaust(RandFaustType.SPARSE, RandFaustType.REAL, nfactors,
                            nfactors,
                            min_dim_sz, max_dim_sz, d)
        filepath = 'test_faust_size'+str(i)+'.mat'
        F.save(filepath)
        stat = os.stat(filepath)
        sizes.append(stat.st_size)
        rcgs.append(F.RCG())
        os.remove(filepath)
        plt.title('File Size vs RCG for Pure-Sparse Fausts \n('+str(nfactors)+' '
                  'random factors '+str(min_dim_sz)+'x'+str(max_dim_sz)+')')
    plt.xlabel('Relative Complexity Gain (RCG)')
    plt.ylabel('File Size (bytes)')
    #plt.scatter(rcgs, sizes, s=1)
    plt.plot(rcgs, sizes)
    plt.show()
