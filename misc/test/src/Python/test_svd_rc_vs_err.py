
# this script aims to show how SVD could be useful to get a better RCG than
# the full matrix while not totally losing the accuracy

# The SVD case here can be seen as a particular situation for a FAµST

from sys import argv
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import norm
from numpy.random import rand
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
        S = np.zeros([s.U.shape[1], s.V.shape[0]])
        if(norm_ord == "fro"): nfroM = 0
        for i in range(0, s.S.shape[0]):
            S[i, i] = s.S[i]
            if(norm_ord == "fro"): nfroM += S[i,i]**2
        if(norm_ord == "fro"): nfroM = np.sqrt(nfroM)
        # assert((abs(M - s.U.dot(S).dot(s.V))/abs(s.M) < .01).all())
        while(s.r < min(s.M.shape[0], s.M.shape[1])):
            #s.Mr = (s.U[:, 0:s.r].dot(S[0:s.r, 0:s.r])).dot(s.V[0:s.r, :])
            #s.errs[0, s.r-1] = norm(s.M-s.Mr, norm_ord)/norm(s.M, norm_ord)
            if(norm_ord == 2):
                err = S[s.r, s.r]/S[0,0]
                #assert(err == s.errs[0, s.r-1])
                #print(s.errs[0, s.r-1], err)
                s.errs[0, s.r-1] = err
            elif(norm_ord == "fro"):
                err = 0
                for j in range(s.r, s.S.shape[0]):
                    err += s.S[j]**2
                err = np.sqrt(err)/nfroM
                #print(s.errs[0, s.r-1], err)
                #assert(err == s.errs[0, s.r-1])
                s.errs[0, s.r-1] = err
            s.rcs[0, s.r-1] = (s.r*(s.M.shape[0]+s.M.shape[1]))
            s.rcs[0, s.r-1] /= s.M.shape[0] * s.M.shape[1]
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
    ths = list()
    for i in range(1, nthreads+1):
        th = TruncSvd(nthreads, M, U, S, V, i, errs, rcs)
        ths.append(th)
        th.start()
    for th in ths:
        th.join()
    plt.scatter(rcs[0, :], errs[0, :], s=1)
    # plt.gca().set_xscale('log', basex=.5)
    #plt.gca().invert_xaxis()
    plt.xlabel('Relative Complexity/Density')
    plt.ylabel('Relative Error')
    plt.title('Relative Error over RC of Truncated SVDs for a dense matrix M ('
              + str(m) + 'x' + str(n)+')')
    f = open("svd_err_vs_rc_output_"+str(nthreads), "w")
    for i in range(0,errs.shape[1]):
        f.write(str(errs[0,i])+" "+str(rcs[0,i])+"\n")
    f.close()
    #plt.plot(rcs, errs)
    plt.show()
