from pyfaust import rand
from pyfaust import *
import numpy as np

F1 = rand(5,50)
F2 = rand(5,50)


concatenate((F1, F2), axis=0)
concatenate((F1, F2), axis=1)

assert(np.allclose(hstack((F1, F2)).toarray(), concatenate((F1, F2),
                                                           axis=1).toarray()))
assert(np.allclose(vstack((F1, F2)).toarray(), concatenate((F1, F2),
                                                           axis=0).toarray()))

A1 = F1.toarray()
A2 = F2.toarray()
assert(np.allclose(np.concatenate((A1, A2), axis=0), concatenate((A1,A2),
                                                              axis=0)))
assert(np.allclose(np.vstack((A1, A2)), vstack((A1,A2))))
assert(np.allclose(np.hstack((A1, A2)), hstack((A1,A2))))

vstack((F1,F2,A2))
