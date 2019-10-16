# -*- coding: utf-8 -*-
# @PYFAUST_LICENSE_HEADER@

## @package pyfaust.fact @brief The pyfaust factorization module
##
##    This module gives access to the main factorization algorithms of
##    FAuST. Those algorithms can factorize a dense matrix to a sparse product
##    (i.e. a Faust object).
##
##    There are several factorization algorithms.
##
##    - The first one is Palm4MSA :
##    which stands for Proximal Alternating Linearized Minimization for
##    Multi-layer Sparse Approximation. Note that Palm4MSA is not
##    intended to be used directly. You should rather rely on the second algorithm.
##
##    - The second one is the Hierarchical Factorization algorithm:
##    this is the central algorithm to factorize a dense matrix to a Faust.
##    It makes iterative use of Palm4MSA to proceed with the factorization of a given
##    dense matrix.
##
##    - The third group of algorithms is for FGFT computing: fgft_palm, fgft_givens, eigtj
##
##
##



import numpy as np, scipy
from scipy.io import loadmat
from scipy.sparse import csr_matrix, csc_matrix
import _FaustCorePy
import pyfaust
import pyfaust.factparams
from pyfaust import Faust
import _FaustCorePy

def svdtj(M, maxiter, tol=0, relerr=True,  nGivens_per_fac=None, verbosity=0):
    """
        Performs a singular value decomposition and returns the left and right singular vectors as Faust transforms.

        NOTE: this function is based on fact.eigtj.

        Args:
            M: a real matrix.
            maxiter: see fact.eigtj
            tol: see fact.eigtj
            relerr: see fact.eigtj
            nGivens_per_fac: see fact.eigtj
            verbosity: see fact.eigtj

        Returns:
            The tuple U,S,V: U*S.todense()*V' being the approximation of M.
                - (sparse diagonal matrix) S the singular values in
                descendant order.
                - (Faust object) U the left-singular transform.
                - (Faust object) V the right-singular transform.

        Example:
            >>> from pyfaust.fact import svdtj
            >>> from numpy.random import rand
            >>> M = rand(128,128)
            >>> U,S,V = svdtj(M, 1024, 64)

     See also:
        eigtj
    """
    from scipy.sparse import spdiags
    from scipy import diag
    from pyfaust import Faust
    from numpy import argsort,sign,eye
    D1, W1 = eigtj(M.dot(M.T), maxiter, tol, relerr,
                   nGivens_per_fac, verbosity)
    D2, W2 = eigtj(M.T.dot(M), maxiter, tol, relerr,  nGivens_per_fac, verbosity)
    S = diag(W1.T*M*W2.todense())
    I = argsort(abs(S))[::-1]
    sign_S = spdiags(sign(S[I]), [0], S.shape[0], S.shape[0])
    S = spdiags(S[I], [0], S.shape[0], S.shape[0])
    S *= sign_S
    Id = eye(S.shape[0])
    U = W1[:,0:S.shape[0]]*Faust([Id[:,I],sign_S.toarray()])
    V = W2[:,0:S.shape[0]]*Faust(Id[:,I])
    return U,S,V

def eigtj(M, maxiter, tol=0, relerr=True,  nGivens_per_fac=None, verbosity=0,
          order='ascend'):
    """
    Runs the truncated Jacobi algorithm to compute the eigenvalues of M (returned in W) and the corresponding transform of eigenvectors (in Faust V columns).

    The output is such that dot(dot(V,diag(W)), V.T) approximates M.

    The trade-off between accuracy and sparsity can be set through the
    parameters maxiter and nGivens_per_fac.

    Args:
        M: (numpy.ndarray) the matrix to diagonalize. Must be real and symmetric.
        maxiter: (int) defines the number of Givens rotations that are computed in
        eigenvector transform V. The number of rotations per factor of V is
        defined by nGivens_per_fac.
        tol: (float) the tolerance error under what the algorithm stops. By default,
        it's zero for not stopping on error criterion.
        relErr: (bool) the type of error stopping criterion. Default to True to use
        relative error, otherwise (False) the absolute error is used.
        nGivens_per_fac: (int) the number of Givens rotations per factor of V, must be
        an integer between 1 to M.shape[0]/2 which is the default value (when
        nGivens_per_fac == None).
        verbosity: (int) the level of verbosity, the greater the value the more info.
        is displayed.
        order: (int) order of eigenvalues, default to 'ascend', other values
        are 'descend' or 'undef'.


    Returns:
        The tuple (W,V):
           - W (numpy.ndarray) the vector of the eigenvalues of M (in ascendant order).
           - V the Faust object representing the approximate eigenvector
            transform. The column V[:, i] is the eigenvector
            corresponding to the eigenvalue W[i].<br/>
            The last factor of V is a permutation matrix.
            The goal of this factor is to apply to the columns of V the same order as
            eigenvalues set in W.

    References:
    [1]   Le Magoarou L., Gribonval R. and Tremblay N., "Approximate fast
    graph Fourier transforms via multi-layer sparse approximations",
    IEEE Transactions on Signal and Information Processing
    over Networks 2018, 4(2), pp 407-420
    <https://hal.inria.fr/hal-01416110>

    Example:
        import numpy as np
        from pyfaust.fact import eigtj
        from scipy.io import loadmat
        from os.path import sep
        from pyfaust.demo import get_data_dirpath
        # get a graph Laplacian to diagonalize
        demo_path = sep.join((get_data_dirpath(),'Laplacian_256_community.mat'))
        data_dict = loadmat(demo_path)
        Lap = data_dict['Lap'].astype(np.float)
        Dhat, Uhat = eigtj(Lap, maxiter=Lap.shape[0]*100)
        # Uhat is the Fourier matrix/eigenvectors approximation as a Faust
        # (200 factors + permutation mat.)
        # Dhat the eigenvalues diagonal matrix approx.
        print(Uhat)
        print(Dhat)

    See also:
        fgft_givens, fgft_palm
    """
    return fgft_givens(M, maxiter, tol, relerr, nGivens_per_fac, verbosity,
                       order)

def fgft_givens(Lap, maxiter, tol=0.0, relerr=True, nGivens_per_fac=None,
                verbosity=0, order='ascend'):
    """
    Computes the FGFT of the Laplacian matrix Lap (using fact.eigtj).

    Args:
        Lap: the Laplacian matrix as a numpy array. Must be real and symmetric.
        maxiter: see fact.eigtj
        tol: see fact.eigtj
        relerr: see fact.eigtj
        nGivens_per_fac: see fact.eigtj
        verbosity: see fact.eigtj
        order: see fact.eigtj

    Returns:
        The tuple (D, FGFT):
        - with FGFT being the Faust object representing
        the Fourier transform and,
        - (numpy.ndarray) D a vector of the eigenvalues of the Laplacian (in
        ascendant order).


    References:
    [1]   Le Magoarou L., Gribonval R. and Tremblay N., "Approximate fast
    graph Fourier transforms via multi-layer sparse approximations",
    IEEE Transactions on Signal and Information Processing
    over Networks 2018, 4(2), pp 407-420.
    <https://hal.inria.fr/hal-01416110>

    See also:
        eigtj, fgft_palm
    """
    if(nGivens_per_fac == None): nGivens_per_fac = int(Lap.shape[0]/2)
    if((isinstance(Lap, np.ndarray) and (Lap.T != Lap).any()) or Lap.shape[0] != Lap.shape[1]):
        raise ValueError(" the matrix/array must be square and should be symmetric.")
    if(not isinstance(maxiter, int)): raise TypeError("maxiter must be a int")
    if(maxiter < 1): raise ValueError("maxiter must be >= 1")
    if(not isinstance(nGivens_per_fac, int)): raise TypeError("nGivens_per_fac must be a int")
    nGivens_per_fac = max(nGivens_per_fac, 1)
    nGivens_per_fac = min(nGivens_per_fac, maxiter)
    if(isinstance(Lap, np.ndarray)):
        core_obj,D = _FaustCorePy.FaustFact.fact_givens_fgft(Lap, maxiter, nGivens_per_fac,
                                                             verbosity, tol,
                                                             relerr, order)
    elif(isinstance(Lap, csr_matrix)):
        core_obj,D = _FaustCorePy.FaustFact.fact_givens_fgft_sparse(Lap, maxiter,
                                                                    nGivens_per_fac,
                                                                    verbosity,
                                                                    tol,
                                                                    relerr,
                                                                    order)
    else:
        raise TypeError("The matrix to diagonalize must be a"
                        " scipy.sparse.csr_matrix or a numpy array.")
    return D, Faust(core_obj=core_obj)

def _check_fact_mat(funcname, M):
    if(not isinstance(M, np.ndarray)):
        raise Exception(funcname+" 1st argument must be a numpy ndarray.")
    if(not isinstance(M[0,0], np.complex) and not isinstance(M[0,0],
                                                            np.float)):
        raise Exception(funcname+" 1st argument must be a float or complex "
                        "ndarray.")
    #if(isinstance(M[0,0], np.complex)):
    #   raise Exception(funcname+" doesn't yet support complex matrix "
    #                   "factorization.")

def palm4msa(M, p, ret_lambda=False):
    """
    Factorizes the matrix M with Palm4MSA algorithm using the parameters set in p.

    Args:
        M: the numpy matrix to factorize.
        p: the ParamsPalm4MSA instance to define the algorithm parameters.
        ret_lambda: set to True to ask the function to return the scale factor (False by default).

    Returns:
        The Faust object resulting of the factorization.
        if ret_lambda == True then the function returns a tuple (Faust, lambda).

    Examples:
    >>> from pyfaust import FaustFactory
    >>> from pyfaust.factparams import ParamsPalm4MSA, ConstraintList, StoppingCriterion
    >>> import numpy as np
    >>> M = np.random.rand(500, 32)
    >>> cons = ConstraintList('splin', 5, 500, 32, 'normcol', 1.0, 32, 32)
    >>> stop_crit = StoppingCriterion(num_its=200)
    >>> param = ParamsPalm4MSA(cons, stop_crit)
    >>> F = palm4msa(M, param)
    >>> F
    Faust size 500x32, density 0.22025, nnz_sum 3524, 2 factor(s):<br/>
    FACTOR 0 (real) SPARSE, size 500x32, density 0.15625, nnz 2500<br/>
    FACTOR 1 (real) SPARSE, size 32x32, density 1, nnz 1024<br/>
    """
    if(not isinstance(p, pyfaust.factparams.ParamsPalm4MSA)):
        raise TypeError("p must be a ParamsPalm4MSA object.")
    _check_fact_mat('palm4msa()', M)
    if(not p.is_mat_consistent(M)):
        raise ValueError("M's number of columns must be consistent with "
                         "the last residuum constraint defined in p. "
                         "Likewise its number of rows must be consistent "
                         "with the first factor constraint defined in p.")
    core_obj, _lambda = _FaustCorePy.FaustFact.fact_palm4msa(M, p)
    F = Faust(core_obj=core_obj)
    if(ret_lambda):
        return F, _lambda
    else:
        return F

def _palm4msa_fgft(Lap, p, ret_lambda=False):
    """
    """
    if(not isinstance(p, pyfaust.factparams.ParamsPalm4MSAFGFT)):
        raise TypeError("p must be a ParamsPalm4MSAFGFT object.")
    _check_fact_mat('_palm4msa_fgft()', Lap)
    if((Lap.T != Lap).any() or Lap.shape[0] != Lap.shape[1]):
        raise ValueError("Laplacian matrix must be square and symmetric.")
    if(not p.is_mat_consistent(Lap)):
        raise ValueError("M's number of columns must be consistent with "
                         "the last residuum constraint defined in p. "
                         "Likewise its number of rows must be consistent "
                         "with the first factor constraint defined in p.")
    core_obj, _lambda, D = _FaustCorePy.FaustFact.fact_palm4msa_fft(Lap, p)
    F = Faust(core_obj=core_obj)
    if(ret_lambda):
        return F, D, _lambda
    else:
        return F, D


def hierarchical(M, p, ret_lambda=False, ret_params=False):
    """
    Factorizes the matrix M with Hierarchical Factorization using the parameters set in p.
    @note This function has its shorthand pyfaust.faust_fact(). For
    convenience you might use it like this:<br/>
            <code>
            from pyfaust import *;
            F = faust_fact(M, p) # equiv. to hierarchical(M, p)
            </code>

    Args:
        M: the numpy matrix to factorize.
        p: is a set of factorization parameters. It might be a fully defined instance of parameters (pyfaust.factparams.ParamsHierarchical) or a simplified expression which designates a pre-defined parametrization:
            - 'squaremat' to use pre-defined parameters typically used to factorize a Hadamard square matrix of order a power of two (see pyfaust.demo.hadamard).
            - ['rectmat', j, k, s] to use pre-defined parameters used for instance in factorization of the MEG matrix which is a rectangular matrix of size m*n such that m < n (see pyfaust.demo.bsl); j is the number of factors, k the sparsity of the main factor's columns, and s the sparsity of rows for all other factors except the residuum (that is the first factor here because the factorization is made toward the left -- is_side_fact_left == true, cf. pyfaust.factparams.ParamsHierarchical).
            <br/>The residuum has a sparsity of P*rho^(num_facts-1). <br/> By default, rho == .8 and P = 1.4. It's possible to set custom values with for example p == ( ['rectmat', j, k, s], {'rho':.4, 'P':.7 }). <br/>The sparsity is here the number of non-zero elements.

        ret_lambda: set to True to ask the function to return the scale factor (False by default).
        ret_params: set to True to ask the function to return the
        ParamsHierarchical instance used (False by default).
        It is useful for consulting what precisely means the
        simplified parametrizations used to generate a
        ParamsHierarchical instance and possibly adjust its attributes to factorize again.


    Returns:
        F the Faust object result of the factorization.<br/>
        if ret_lambda == True (and ret_params == False), then the function
        returns the tuple (F,_lambda) (_lambda is the scale factor at the
        end of factorization).<br/>
        if ret_params == True (and ret_lambda == False), then the function
        returns the tuple (F, p) (p being the ParamsHierarchical
        instance really used by the algorithm).<br/>
        if ret_lambda == True and ret_params == True, then the function
        returns the tuple (F, _lambda, p).

    Examples:
        <b> 1. Fully Defined Parameters for a Random Matrix Factorization </b>
        >>> from pyfaust.fact import hierarchical
        >>> from pyfaust.factparams import ParamsHierarchical, ConstraintList, StoppingCriterion
        >>> import numpy as np
        >>> M = np.random.rand(500, 32)
        >>> fact_cons = ConstraintList('splin', 5, 500, 32, 'sp', 96, 32, 32, 'sp', 96, 32, 32)
        >>> res_cons = ConstraintList('normcol', 1, 32, 32, 'sp', 666, 32, 32, 'sp', 333, 32, 32)
        >>> stop_crit1 = StoppingCriterion(num_its=200)
        >>> stop_crit2 = StoppingCriterion(num_its=200)
        >>> param = ParamsHierarchical(fact_cons, res_cons, stop_crit1, stop_crit2)
        >>> F = hierarchical(M, param)
        Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorization
        1/3<br/>
        Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorization
        2/3<br/>
        Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorization
        3/3<br/>
        >>> F
        Faust size 500x32, density 0.189063, nnz_sum 3025, 4 factor(s):
            - FACTOR 0 (real) SPARSE, size 500x32, density 0.15625, nnz 2500
            - FACTOR 1 (real) SPARSE, size 32x32, density 0.09375, nnz 96
            - FACTOR 2 (real) SPARSE, size 32x32, density 0.09375, nnz 96
            - FACTOR 3 (real) SPARSE, size 32x32, density 0.325195, nnz 333

        <b>2. Simplified Parameters for Hadamard Factorization</b>

        >>> from pyfaust import wht
        >>> from pyfaust.fact import hierarchical
        >>> from numpy.linalg import norm
        >>> # generate a Hadamard Faust of size 32x32
        >>> FH = wht(32)
        >>> H = FH.toarray() # the full matrix version
        >>> # factorize it
        >>> FH2 = hierarchical(H, 'squaremat');
        >>> # test the relative error
        >>> (FH-FH2).norm('fro')/FH.norm('fro') # the result is 1.1015e-16, the factorization is accurate
        >>> FH
        Faust size 32x32, density 0.3125, nnz_sum 320, 5 factor(s):
        - FACTOR 0 (real) SPARSE, size 32x32, density 0.0625, nnz 64
        - FACTOR 1 (real) SPARSE, size 32x32, density 0.0625, nnz 64
        - FACTOR 2 (real) SPARSE, size 32x32, density 0.0625, nnz 64
        - FACTOR 3 (real) SPARSE, size 32x32, density 0.0625, nnz 64
        - FACTOR 4 (real) SPARSE, size 32x32, density 0.0625, nnz 64

        >>> FH2
        Faust size 32x32, density 0.3125, nnz_sum 320, 5 factor(s):
            - FACTOR 0 (real) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 1 (real) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 2 (real) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 3 (real) SPARSE, size 32x32, density 0.0625, nnz 64
            - FACTOR 4 (real) SPARSE, size 32x32, density 0.0625, nnz 64

       <b> 3. Simplified Parameters for a Rectangular Matrix Factorization
       (the BSL demo MEG matrix)</b>
       >>> from pyfaust import *
       >>> from scipy.io import loadmat
       >>> from pyfaust.demo import get_data_dirpath
       >>> d = loadmat(get_data_dirpath()+'/matrix_MEG.mat')
       >>> MEG = d['matrix']
       >>> num_facts = 9
       >>> k = 10
       >>> s = 8
       >>> MEG16 = hierarchical(MEG, ['rectmat', num_facts, k, s])
       >>> MEG16
       Faust size 204x8193, density 0.0631655, nnz_sum 105573, 9 factor(s):
           - FACTOR 0 (real) SPARSE, size 204x204, density 0.293613, nnz
           12219
           - FACTOR 1 (real) SPARSE, size 204x204, density 0.0392157, nnz
           1632
           - FACTOR 2 (real) SPARSE, size 204x204, density 0.0392157, nnz
           1632
           - FACTOR 3 (real) SPARSE, size 204x204, density 0.0392157, nnz
           1632
           - FACTOR 4 (real) SPARSE, size 204x204, density 0.0392157, nnz
           1632
           - FACTOR 5 (real) SPARSE, size 204x204, density 0.0392157, nnz
           1632
           - FACTOR 6 (real) SPARSE, size 204x204, density 0.0392157, nnz
           1632
           - FACTOR 7 (real) SPARSE, size 204x204, density 0.0392157, nnz
           1632
           - FACTOR 8 (real) SPARSE, size 204x8193, density 0.0490196, nnz
           81930

       >>> # verify the constraint k == 10, on column 4
       >>> count_nonzero(MEG16.factors(8)[:,4].toarray())
       10
       >>> # now verify the s constraint is respected on MEG16 factor 1
       >>> count_nonzero(MEG16.factors(1).toarray())/MEG16.shape[0]
       8.0

    """
    p = _prepare_hierarchical_fact(M,p, "hierarchical", ret_lambda,
                              ret_params)
    core_obj,_lambda = _FaustCorePy.FaustFact.fact_hierarchical(M, p)
    F = Faust(core_obj=core_obj)
    ret_list = [ F ]
    if(ret_lambda):
        ret_list += [ _lambda ]
    if(ret_params):
        ret_list += [ p ]
    if(ret_lambda or ret_params):
        return ret_list
    else:
        return F


def _prepare_hierarchical_fact(M, p, callee_name, ret_lambda, ret_params,
                               M_name='M'):
    """
    Utility func. for hierarchical() and fgft_palm().
    Among other checkings, it sets parameters from simplified ones.
    """
    from pyfaust.factparams import (ParamsHierarchical,
                                    ParamsFactFactory)
    if(not isinstance(p, ParamsHierarchical) and
       ParamsFactFactory.is_a_valid_simplification(p)):
        p = ParamsFactFactory.createParams(M, p)
    if(not isinstance(p, ParamsHierarchical)):
        raise TypeError("p must be a ParamsHierarchical object.")
    _check_fact_mat(''+callee_name+'()', M)
    if(not isinstance(ret_lambda, bool)):
        raise TypeError("ret_lambda must be a bool.")
    if(not isinstance(ret_params, bool)):
        raise TypeError("ret_params must be a bool.")
    if(not p.is_mat_consistent(M)):
        raise ValueError("M's number of columns must be consistent with "
                         "the last residuum constraint defined in p. "
                         "Likewise its number of rows must be consistent "
                         "with the first factor constraint defined in p.")
    return p

# experimental block start
def hierarchical_constends(M, p, A, B, ret_lambda=False, ret_params=False):
    """
    Tries to approximate M by \f$ A \prod_j S_j B \f$ using the hierarchical.

    It needs to add one additional residuum constraint and also one factor
    constraint into p relatively to the number of constraints you would
    have defined if you just factorized M into prod \f$ \prod_j S_j \f$.

    Examples:
        from pyfaust import rand
        from pyfaust import hierarchical
        import numpy as np
        from numpy.random import rand
        from pyfaust.factparams import (ParamsHierarchical, ConstraintList,
                                        StoppingCriterion)
        M = rand(32,32)
        A = rand(32,32)
        B = rand(32,32)
        stop_crit1 = StoppingCriterion(num_its=200)
        stop_crit2 = StoppingCriterion(num_its=200)
        fact_cons = ConstraintList('splin', 5, 32, 32, 'sp', 96, 32, 32, 'sp', 96,
                                   32, 32,
                                   'sp', 225, 32, 32)
        res_cons = ConstraintList('normcol', 1, 32, 32,
                                  'normcol', 1, 32, 32, 'sp', 666, 32, 32, 'sp',
                                  333, 32, 32)
        param = ParamsHierarchical(fact_cons, res_cons, stop_crit1, stop_crit2)
        F, _lambda = hierarchical_constends(M, param, A, B, ret_lambda=True)

        assert(np.allclose(F.factors(0).toarray(), A))
        assert(np.allclose(F.factors(5).toarray(), B))

    """
    from pyfaust.factparams import ConstraintList
    from pyfaust.factparams import ParamsHierarchical
    A = np.asfortranarray(A)
    B = np.asfortranarray(B)
    consA = ConstraintList('const', A, *A.shape)
    consB = ConstraintList('const', B, *B.shape)
    consts = p.constraints
    nconsts = len(consts)
    # consts : factor constraints + residuum constraints
    fac_cons = [ consts[i] for i in range(0, p.num_facts-1) ]
    res_cons = [ consts[i] for i in range(p.num_facts-1, len(consts))]
    assert(nconsts == len(fac_cons) + len(res_cons))
    # add two CONST(ants) factor constraint: A and B to the constraints
    # depending of the factorization direction we need to switch A and B
    # constraint positions
    if(p.is_fact_side_left):
        new_fact_cons = consB + ConstraintList(*fac_cons)
        new_res_cons = ConstraintList(*res_cons) + consA
    else:
        new_fact_cons = consA + ConstraintList(*fac_cons)
        new_res_cons = ConstraintList(*res_cons) + consB
    assert(nconsts == len(new_fact_cons) + len(new_res_cons) - 2)
    p = ParamsHierarchical(new_fact_cons, new_res_cons,
                               p.stop_crits[0], p.stop_crits[1],
                               p.is_update_way_R2L, p.init_lambda,
                               p.step_size, p.constant_step_size,
                               p.is_fact_side_left,
                               p.is_verbose)
    F, _lambda = hierarchical(M, p, ret_lambda=True,
                                          ret_params=False)
    FA_lambda = F.factors(0)
    FA = FA_lambda/_lambda
    new_F_factors = [ FA ]
    F_without_A = F.factors(range(1,F.numfactors()))
    F_without_A = F_without_A * _lambda
    new_F_factors += [ F_without_A.factors(i) for i in
                      range(0,F_without_A.numfactors()) ]
    new_F = Faust(new_F_factors)
    F = new_F
    ret_list = [F]
    if(ret_lambda):
        ret_list += [ _lambda ]
    if(ret_params):
        ret_list += [ p ]
    if(ret_lambda or ret_params):
        return ret_list
    else:
        return F
# experimental block end

# experimental block start
def palm4msa_constends(M, p, A, B=None, ret_lambda=False):
    """
    Tries to approximate M by \f$ A \prod_j S_j B\f$ using palm4msa (B being optional).


    Example:
        from pyfaust import rand
        from pyfaust.fact import palm4msa_constends
        from pyfaust.factparams import (ParamsPalm4MSA, ConstraintList,
                                        StoppingCriterion)
        import numpy as np
        from numpy.random import rand
        from random import randint
        n = 32
        k = 16
        M = rand(n, n)
        A = rand(n, k)
        B = rand(k, n)
        stop_crit = StoppingCriterion(num_its=100)
        consts = ConstraintList('spcol', randint(1, A.shape[1]), A.shape[1],
                                                 B.shape[0])
        param = ParamsPalm4MSA(consts, stop_crit)
        F = palm4msa_constends(M, param, A, B)

        assert(np.allclose(F.factors(0).toarray(), A))
        assert(np.allclose(F.factors(2).toarray(), B))

    """
    from pyfaust.factparams import ConstraintList
    from pyfaust.factparams import ParamsPalm4MSA
    A = np.asfortranarray(A)
    B = np.asfortranarray(B)
    consA = ConstraintList('const', A, *A.shape)
    consts = ConstraintList(*p.constraints)
    new_consts = consA + consts
    if(isinstance(B, np.matrix) or isinstance(B, np.ndarray)):
        consB = ConstraintList('const', B, *B.shape)
        new_consts = new_consts + consB
    p = ParamsPalm4MSA(new_consts, stop_crit=p.stop_crit, init_facts=p.init_facts,
                       is_update_way_R2L=p.is_update_way_R2L,
                       init_lambda=p.init_lambda, step_size=p.step_size,
                       constant_step_size = p.constant_step_size,
                       is_verbose = p.is_verbose)
    F, _lambda = palm4msa(M, p, ret_lambda=True)
    F = \
    Faust([F.factors(0)/_lambda]+[F.factors(1)*_lambda]+
          [F.factors(i)
           for i in
           range(2,F.numfactors())])
    if(ret_lambda):
        return F, _lambda
    else:
        return F
# experimental block end

def fgft_palm(U, Lap, p, init_D=None, ret_lambda=False, ret_params=False):
    """
    Computes the FGFT for the Fourier matrix U (it should be the eigenvectors of the Laplacian Lap).


    NOTE: this algorithm is a variant of hierarchical.


    Example:
        from pyfaust import FaustFactory
        from pyfaust.factparams import *
        from scipy.io import loadmat, savemat
        from pyfaust.demo import get_data_dirpath
        from os.path import sep
        from numpy.linalg import eig, eigh, norm
        from numpy import sort, argsort, log2, size, copy, diag

        d = loadmat(sep.join((get_data_dirpath(),'Laplacian_256_ring.mat')))
        Lap = d['Lap'].astype('float')


        D, U = eig(Lap)

        indices = argsort(D)
        D = D[indices].astype('float')
        U = U[:,indices].astype('float')

        print(D.shape, type(D))
        print("eig(Lap), U error:", norm(Lap.dot(U)-U.dot(diag(D))))
        dim = Lap.shape[0]
        nfacts = int(round(log2(dim))-3)
        over_sp = 1.5 # sparsity overhead
        dec_fact = .5 # decrease of the residum sparsity
        fact_cons, res_cons = [], []
        for j in range(1, nfacts):
            fact_cons += [ ConstraintInt('sp',dim,dim,
                            min(int(round(dec_fact**j*dim**2*over_sp)),size(Lap)))
                        ]
            res_cons += [
                ConstraintInt('sp',
                                dim,
                                dim,
                                min(int(round(2*dim*over_sp)),size(Lap)))
            ]

        params = ParamsHierarchical(fact_cons,
            res_cons,
            StoppingCriterion(num_its=50),
            StoppingCriterion(num_its=100),
            step_size=1.0000e-06,
            constant_step_size=True,
            init_lambda=1.0,
            is_fact_side_left=False)


        Lap = Lap.astype(float)
        Uhat,Dhat = fgft_palm(Lap, U, params, init_D=D)

        err_U = (Uhat-U).norm()/norm(U)
        err_Lap = norm(Uhat.todense()*diag(Dhat)*Uhat.T.todense()-Lap)/norm(Lap)
        print(norm(diag(Dhat), 2))

        print("err_U:", err_U)
        print("err_Lap:", err_Lap)


        #end of output:
        #	eig(Lap), U error: 1.36688422173e-13
        #	Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 1/3
        #	Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 2/3
        #	Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 3/3
        #	err_U: 1.00138959974
        #	err_Lap: 0.997230709335


    Args:
        Lap: (numpy.ndarray) The Laplacian matrix.
        U: (numpy.ndarray) The Fourier matrix.
        init_D: (numpy.ndarray) The initial diagonal vector. if None it will be the ones() vector by default.
        p: (ParamsHierarchical) The PALM hierarchical algorithm parameters.
        ret_lambda: (bool) True to return the lambda scale factor used in PALM algorithm.
        ret_params: (bool) True to return the parameters used by the hierarchical PALM algorithm (normally it's p except if p is a simplifed form of parameters -- not instance of ParamsHierarchical).

    Returns:
        FGFT: (Faust object) The Fourier transform.
        - if ret_lambda == True (and ret_params == False), then the function returns the
        tuple (FGFT,_lambda) (_lambda is the scale factor at the end of
        factorization).
        - if ret_params == True (and ret_lambda == False), then the function returns the
        tuple (FGFT, p) (p being the ParamsHierarchical instance really used by the
        algorithm).
        - if ret_lambda == True and ret_params == True, then the function returns the
        tuple (FGFT, _lambda, p).

    See also:
        hierarchical, fgft_givens, eigtj

    References:
        - [1]   Le Magoarou L., Gribonval R. and Tremblay N., "Approximate fast
        graph Fourier transforms via multi-layer sparse approximations",
        IEEE Transactions on Signal and Information Processing
        over Networks 2018, 4(2), pp 407-420
        <https://hal.inria.fr/hal-01416110>
        - [2] Le Magoarou L. and Gribonval R., "Are there approximate Fast
        Fourier Transforms on graphs ?", ICASSP, 2016.  <https://hal.inria.fr/hal-01254108>


    """
    from pyfaust.factparams import _init_init_D

    p = _prepare_hierarchical_fact(U, p, "fgft_palm", ret_lambda,
                              ret_params, 'U')
    if(init_D.dtype != Lap.dtype or Lap.dtype != U.dtype or U.dtype != init_D.dtype):
        raise ValueError("All the numpy arrays must be of the same dtype.")

    init_D = _init_init_D(init_D, U.shape[0])
    core_obj, _lambda, D = _FaustCorePy.FaustFact.fact_hierarchical_fft(U, Lap, p,
                                                                    init_D)
    F = Faust(core_obj=core_obj)
    ret_list = [ F, D ]
    if(ret_lambda):
        ret_list += [ _lambda ]
    if(ret_params):
        ret_list += [ p ]
    return ret_list
