
# this script aims to show how SVD could be useful to get a better RCG than
# the dense matrix while not totally losing the accuracy

# The SVD case here can be seen as a particular Faust

from sys import argv
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import norm
from numpy.random import rand
from numpy import diag,copy,savetxt,empty
from threading import Thread
from scipy.linalg import svd


class TruncSvd(Thread):

    def __init__(self, offset, M, U, S, V, r, errs, rcs):
        Thread.__init__(self)
        self.offset = offset
        self.M = M
        self.U = U
        self.V = V
        self.S = S
        self.r = r
        self.errs = errs
        self.rcs = rcs

    def run(s):
        global norm_ord
        S = diag(s.S)
        while(s.r < min(s.M.shape[0], s.M.shape[1])):
            Sr = copy(S)
            Sr[s.r:, s.r:] = 0
            if(norm_ord == 2):
                err = S[s.r, s.r]/S[0,0]
                s.errs[0, s.r-1] = err
                #print(norm(M-U@Sr@V,2)/norm(M,2),err)
            elif(norm_ord == "fro"):
                err = norm(s.S[s.r:s.S.shape[0]])/norm(s.S)
                s.errs[0, s.r-1] = err
                #print(norm(M-U@Sr@V,'fro')/norm(M,'fro'),err)
            s.rcs[0, s.r-1] = (s.r*(s.M.shape[0]+s.M.shape[1]))
            s.rcs[0, s.r-1] /= s.M.size
            s.r += s.offset

if __name__ == '__main__':
    if(len(argv) < 4):
        print("USAGE: ", argv[0], "<num_lines> <num_cols> <norm_ord>")
        print("norm_ord is fro or 2")
        exit(1)
    m = int(argv[1])
    n = int(argv[2])
    global norm_ord
    if(argv[3] == "2"):
        norm_ord=2
    elif(argv[3] == "fro"):
        norm_ord="fro"
    else:
        raise("Error: norm must be 2 or fro.")
    nthreads = 8
    errs = np.zeros([1, min(m, n)])
    rcs = np.zeros([1, min(m, n)])
    M = rand(m, n)
    U, S, V = svd(M)
    ths = []
    for i in range(1, nthreads+1):
        th = TruncSvd(nthreads, M, U, S, V, i, errs, rcs)
        ths.append(th)
        th.start()
    for th in ths:
        th.join()
    plt.scatter(rcs[0, :], errs[0, :], s=1)
    plt.xlabel('Density')
    plt.ylabel('Relative Error')
    plt.title('Relative Error vs Density of Truncated SVDs for a dense matrix M ('
              + str(m) + 'x' + str(n)+')')
    xy = empty((errs.shape[1],2))
    xy[:,0] = errs[0,:]
    xy[:,1] = rcs[0,:]
    savetxt("svd_err_vs_rc_output_"+str(nthreads)+"_err_"+str(norm_ord), xy)
    plt.tight_layout()
    plt.show()
