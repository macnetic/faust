from pyfaust.fact import palm4msa
from pyfaust.factparams import ParamsPalm4MSA, StoppingCriterion
from pyfaust.proj import splin, spcol
import pyfaust
from numpy.linalg import norm
import numpy as np

num_its = 100 # the total number of iterations to run PALM4MSA

# Generate a matrix to factorize
d = 64
S = pyfaust.rand(d, d, num_factors=2, density=0.1, per_row=True)
M = S.toarray()

# Set the set of parameters for PAL4MSA
projs = [splin((d,d), 10), spcol((d,d), 5)]

# set just one iteration of the step-by-step mode
stop_crit = StoppingCriterion(num_its=1)

# pack all these parameters into a ParamsPalm4MSA
p = ParamsPalm4MSA(projs, stop_crit)

rel_errs = [] # relative errors for all iterations

# Runs one iteration at a time and compute the relative error
for i in range(num_its):
    F, scale = palm4msa(M, p, ret_lambda=True)
    # backup the error
    rel_errs += [(F-M).norm()/np.linalg.norm(M)]
    # retrieve the factors from the Faust F obtained after the one-iteration execution
    # it's needed to convert them explicitely as Fortran (that is column major
    # order) numpy array, otherwise it would fail the next iteration
    facts = [np.asfortranarray(F.factors(0).toarray()),
             np.asfortranarray(F.factors(1).toarray())]
    # F is multiplied by lambda, we need to undo that before continuing to the
    # next iteration, the factor which was multiplied by scale can be any of
    # them but because all the others are normalized, it's easy to deduce which
    # one it is
    # NOTE: remember we want to set the factors exactly as they were at the end of
    # iteration 0, that's why all must be normalized (by dividind by scale)
    if norm(facts[0])  > norm(facts[1]):
        facts[0] /= scale
    else:
        facts[1] /= scale
        # all is ready to pack the parameters and initial state (the factors and
        # the scale) in a new ParamsPalm4MSA instance
        p = ParamsPalm4MSA(projs, stop_crit,
                           init_facts=facts,
                           init_lambda=scale)


print("relative errors along the iterations (step-by-step mode):", rel_errs)


stop_crit = StoppingCriterion(num_its=num_its)
p = ParamsPalm4MSA(projs, stop_crit)
G = palm4msa(M, p)
print("Relative error when running all PALM4MSA iterations at once: ", (G-M).norm() / norm(M))
print("Last relative error obtained when running PALM4MSA "
      "iteration-by-iteration: ", rel_errs[-1])
print("Relative error comparing the final Fausts obtained either by in"
      " step-by-step PALM4MSA versus all-iterations-at-once PALM4MSA: ", (F-G).norm() / G.norm())
