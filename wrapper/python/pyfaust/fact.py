# -*- coding: utf-8 -*-
# @PYFAUST_LICENSE_HEADER@

## @package pyfaust.fact @brief The pyfaust factorization module
##
##    This module gives access to the main factorization algorithms of
##    FAuST. These algorithms can factorize a dense matrix into a sparse product
##    (i.e. a Faust object). A few of them are only available in experimental
##    packages.
##
##    There are several factorization algorithms.
##
##    - The first one is Palm4MSA :
##    which stands for Proximal Alternating Linearized Minimization for
##    Multi-layer Sparse Approximation. Note that Palm4MSA is not
##    intended to be used directly. You should rather rely on the second algorithm.
##
##    - The second one is the Hierarchical Factorization algorithm:
##    this is the central algorithm to factorize a dense matrix into a Faust.
##    It makes iterative use of Palm4MSA to proceed with the factorization of a given
##    dense matrix.
##
##    - The third group of algorithms is for approximate eigenvalue decomposition (eigtj) and singular value decomposition (svdtj).
##
##
##


import numpy as np, scipy
from scipy.io import loadmat
from scipy.sparse import csr_matrix, csc_matrix, eye as seye, kron as skron
import pyfaust
import pyfaust.factparams
from pyfaust import Faust
from pyfaust import bitrev_perm
import _FaustCorePy
import warnings
from pyfaust.tools import _sanitize_dtype

# experimental block start
def svdtj2(M, nGivens, tol=0, relerr=True,  nGivens_per_fac=None, verbosity=0,
          enable_large_Faust=False):
    """
        Performs an approximate singular value decomposition and returns the left and right singular vectors as Faust transforms.

        NOTE: this function is based on fact.eigtj which relies on the
        truncated Jacobi algorithm, hence the 'tj' in the name.

        Args:
            M: a real or complex, dense or sparse (csr) matrix.
            nGivens: see fact.eigtj
            tol: see fact.eigtj
            relerr: see fact.eigtj
            nGivens_per_fac: see fact.eigtj
            verbosity: see fact.eigtj

        Returns:
            The tuple U,S,V: U*S.toarray()@V' being the approximation of M.
                - (sparse diagonal matrix) S the singular values in
                descending order.
                - (Faust object) U the left-singular transform.
                - (Faust object) V the right-singular transform.

        Example:
            >>> from pyfaust.fact import svdtj
            >>> from numpy.random import rand
            >>> M = rand(128,128)
            >>> U,S,V = svdtj(M, 1024, nGivens_per_fac=64)

     See also:
        eigtj
    """
    from scipy.sparse import spdiags
    from scipy import diag
    from pyfaust import Faust
    from numpy import argsort,sign,eye, sort
    if verbosity > 1:
        print("input eigtj nGivens, tol, order, relerr, "
              "nGivens_per_fac, verbosity, enable_large_Faust:" , nGivens, tol,
              'descend', relerr, nGivens_per_fac, verbosity, enable_large_Faust)
    D1, W1 = eigtj(M.dot(M.T.conj()), nGivens, tol, 'descend', relerr,
                   nGivens_per_fac, verbosity, enable_large_Faust)
    D2, W2 = eigtj(M.T.conj().dot(M), nGivens, tol, 'descend', relerr,
                   nGivens_per_fac, verbosity, enable_large_Faust)
    S = diag(W1.T.conj() @ M @ W2)
    I = argsort(abs(S))[::-1]
    sign_S = spdiags(sign(S[I]), [0], W1.shape[1], W1.shape[1])
    S = spdiags(S[I] * sign(S[I]), [0], W1.shape[1], W2.shape[1])
    if len(I) == W1.shape[1]:
        PU = eye(len(I))[:, I]
    else:
        # W1.shape[1] > len(I)
        J = np.hstack((I, np.arange(len(I), W1.shape[1])))
        PU = eye((W1.shape[1]))[:,J]
    if len(I) == W2.shape[1]:
        PV = eye(len(I))[:, I]
    else:
        # W2.shape[0] > len(I)
        J = np.hstack((I, np.arange(len(I), W2.shape[1])))
        PV = eye((W2.shape[1]))[:,J]
    U = W1.left(len(W1)-2) @ Faust(W1.factors(len(W1)-1) @ PU @ sign_S.toarray())
    V = W2.left(len(W2)-2) @ Faust(W2.factors(len(W2)-1) @ PV)
    return U,S,V
# experimental block end


def svdtj(M, nGivens=None, tol=0, err_period=100, relerr=True,
          nGivens_per_fac=None, enable_large_Faust=False, **kwargs):
    """
        Performs a singular value decomposition and returns the left and right singular vectors as Faust transforms.

        NOTE: this function is based on fact.eigtj which relies on the truncated Jacobi algorithm, hence the 'tj' in the name. See below the examples for further details on how svdtj is defined using eigtj.

        Args:
            M: a real or complex matrix (np.ndarray or scipy.sparse.csr_matrix). The dtype must be float32, float64
            or complex128 (the dtype might have a large impact on performance).
            nGivens: (int or tuple(int, int)) defines the number of Givens rotations
            that will be used at most to compute U and V.
            If it is an integer, it will apply both to U and V.
            If it is a tuple of two integers as nGivens = (JU, JV), JU
            will be the limit number of rotations for U and JV the same for V.
            nGivens argument is optional if tol is set but becomes mandatory otherwise.
            tol: (float) this is the error target on U S V' against M.
            if error <= tol, the algorithm stops. See relerr below for the error formula.
            This argument is optional if nGivens is set, otherwise it becomes mandatory.
            err_period: (int) it defines the period, in number of factors of U
            or V, the error is compared to tol (reducing the period spares
            some factors but increases slightly the computational cost because the error
            is computed more often).
            relerr: (bool) defines the type of error used to evaluate the stopping
            criterion defined by tol. true for a relative error (\f$\| U' S V -
            M\|_F \over \| M \|_F\f$) otherwise this is an absolute error (\f$\| U' S V - M\|_F\f$).
            nGivens_per_fac: (int or tuple(int, int)) this argument is the number of Givens
            rotations to set at most by factor of U and V.
            If this is an integer it will be the same number of U and V.
            Otherwise, if it is a tuple of integers (tU, tV), tU will be the number
            of Givens rotations per factor for U and tV the same for V.
            By default, this parameter is maximized for U and V,
            i.e. tU = M.shape[0] / 2, tV = M.shape[1] / 2.
            enable_large_Faust: see eigtj.


        Returns:
            The tuple U,S,V: such that U*numpy.diag(S)*V.H is the approximate of M.
            - (np.array vector) S the singular values in descending order.
            - (Faust objects) U,V orthonormal Fausts.
        <br/>
        Examples:
            >>> from pyfaust.fact import svdtj
            >>> from numpy.random import rand, seed
            >>> import numpy as np
            >>> from scipy.sparse import spdiags
            >>> seed(42) # just for reproducibility
            >>> M = rand(16, 32)
            >>> # Factoring by specifying the number of Givens rotations
            >>> U1, S1, V1 = svdtj(M, 4096, enable_large_Faust=True)
            >>> S1_ = spdiags(S1, [0], U1.shape[0], V1.shape[0])
            >>> np.allclose(U1 @ S1_ @ V1.H, M)
            True
            >>> # Specifying a different number of rotations for U and V
            >>> # Because U is smaller it should need less rotations
            >>> U2, S2, V2 = svdtj(M, (2400, 3200), enable_large_Faust=True)
            >>> S2_ = spdiags(S2, [0], U2.shape[0], V2.shape[0])
            >>> np.allclose(U2 @ S2_ @ V2.H, M)
            True
            >>> # Factoring according to an approximate accuracy target
            >>> U3, S3, V3 = svdtj(M, tol=1e-12, enable_large_Faust=False)
            >>> S3_ = spdiags(S3, [0], U3.shape[0], V3.shape[0])
            >>> np.allclose(U3 @ S3_ @ V3.H, M)
            True
            >>> # verify the relative error is lower than 1e-12
            >>> np.linalg.norm(U3 @ S3_ @ V3.H - M) / np.linalg.norm(M)
            9.122446623891136e-13
            >>> # try with an absolute tolerance (the previous one was relative to M)
            >>> U4, S4, V4 = svdtj(M, tol=1e-6, relerr=False, enable_large_Faust=False)
            >>> S4_ = spdiags(S4, [0], U4.shape[0], V4.shape[0])
            >>> np.allclose(U4 @ S4_ @ V4.H, M)
            True
            >>> # verify the absolute error is lower than 1e-6
            >>> np.linalg.norm(U4 @ S4_ @ V4.H - M)
            1.2044207331117425e-11
            >>> # try a less accurate approximate to get less factors
            >>> U5, S5, V5 = svdtj(M, nGivens=(256, 512), tol=1e-1, relerr=True, enable_large_Faust=False)
            >>> S5_ = spdiags(S5, [0], U5.shape[0], V5.shape[0])
            >>> # verify the absolute error is lower than 1e-3
            >>> np.linalg.norm(U5 @ S5_ @ V5.H - M) / np.linalg.norm(M)
            0.09351811486725303
            >>> ### Let's see the lengths of the different U, V Fausts
            >>> len(V1) # it should be 4096 / nGivens_per_fac, which is (M.shape[1] // 2) = 256
            256
            >>> len(U1) # it should be 4096 / nGivens_per_fac, which is (M.shape[0] // 2) = 512
            100
            >>> # but it is not, svdtj stopped automatically extending U1 because the error stopped enhancing
            >>> # (it can be verified with verbosity=1)
            >>> (len(U2), len(V2))
            (100, 200)
            >>> (len(U3), len(V3))
            (64, 200)
            >>> (len(U4), len(V4))
            (64, 200)
            >>> # not surprisingly U5 and V5 use the smallest number of factors (nGivens and tol were the smallest)
            >>> (len(U5), len(V5))
            (32, 32)

        Explanations:

        If we call svdtj on the matrix M, it makes two internal calls to eigtj. In Python it would be:

        1.  D1, W1 = eigtj(M.dot(M.H), next_args...)
        2.  D2, W2 = eigtj(M.H.dot(M), next_args...)

        It gives the following equalities (ignoring the fact that eigtj computes approximations):
        \f[
            W_1 D_1 W_1^* = M M^*
        \f]
        \f[
                W_2 D_2 W_2^* = M^* M
        \f]
        But because of the SVD \f$ M = USV^* \f$ we also have:
        \f[MM^* = U S V^* V S U^* = U S^2 U^* = W_1 D_1 W_1^*\f]
        \f[M^* M = V S U^* U S V^* = V S^2 V^* = W_2 D_2  W_2^*\f]
        It allows to identify the left singular vectors of M to W1,
        and likewise the right singular vectors to W2.

        To compute a consistent approximation of S we observe that U and V are orthogonal hence \f$ S  = U^* M V \f$ so we ignore the off-diagonal coefficients of the approximation and take \f$ S = diag(U^* M V)  \approx diag(W_1^* M W_2)\f$

        The last step performed by svdtj() is to sort the singular values of S
        in descending order and build a signed permutation matrix to order the
        left singular vectors of W1 accordingly. The -1 elements of the signed
        permutation matrix allow to change the sign of each negative values of
        S by reporting it on the corresponding left singular vector (\f$
        \sigma_i v_i = (-\sigma_i) (-v_i )\f$).<br/>
        To sum up W1 is replaced by W1 P and W2 by W2 abs(P) (because W2 also
        needs to be ordered), with P the signed permutation resulting of the
        descending sort of S. The resulting transforms/Fausts W1 and W2 are
        returned by svdtj along with the ordered S. Note that the permutation
        factor P (resp. abs(P)) is fused with the rightmost factor of the Faust
        object W1 (resp. W2).

     See also:
        eigtj
    """

    # for nGivens_per_fac and nGivens checking
    def pair_of_posreals_or_posreal(v):
        return isinstance(v, (tuple, list)) and \
                len(v) == 2 and \
                (np.isreal(v[0]) and v[0] > 0) and \
                (np.isreal(v[1]) and v[1] > 0) or \
                np.isreal(v) and v > 0

    if(nGivens == None):
        if(tol == 0):
            raise Exception("You must specify nGivens or tol argument"
                            " (to define a stopping  criterion)")
        nGivens = 0
    else:
        # verify user value
        if not pair_of_posreals_or_posreal(nGivens):
            raise ValueError('nGivens must be a tuple of two positive integers '
                             'or an integer')

    if(nGivens_per_fac == None):
        # default nGivens_per_fac
        nGivens_per_fac = (int(M.shape[0] // 2),
                           int(M.shape[1] // 2))
    else:
        # verify user value
        if not pair_of_posreals_or_posreal(nGivens_per_fac):
            raise ValueError('nGivens_per_fac must be a tuple of two positive integers '
                             'or an integer')
    if('verbosity' in kwargs.keys()):
        verbosity = kwargs['verbosity']
        if(not isinstance(verbosity, int)): raise TypeError('verbosity must be'
                                                            ' a int')
    else:
        verbosity = 0

    is_real = np.empty((1,))
    M = _check_fact_mat('svdtj()', M, is_real)

    assert np.isreal(err_period)
    err_period = int(err_period)

    svdtj_args = [M, nGivens, nGivens_per_fac, verbosity, tol, relerr,
                  enable_large_Faust, err_period]
    if is_real:
        is_float = M.dtype == 'float32'
        if(isinstance(M, np.ndarray)):
            if is_float:
                Ucore, S, Vcore =  _FaustCorePy.FaustAlgoGenGivensFlt.svdtj(*svdtj_args)
            else:
                Ucore, S, Vcore =  _FaustCorePy.FaustAlgoGenGivensDbl.svdtj(*svdtj_args)
        elif(isinstance(M, csr_matrix)):
            if is_float:
                Ucore, S, Vcore =  _FaustCorePy.FaustAlgoGenGivensFlt.svdtj_sparse(*svdtj_args)
            else:
                Ucore, S, Vcore =  _FaustCorePy.FaustAlgoGenGivensDbl.svdtj_sparse(*svdtj_args)
        else:
            raise ValueError("invalid type for M (first argument): only np.ndarray "
                             "or scipy.sparse.csr_matrix are supported.")
    else:
        if(isinstance(M, np.ndarray)):
            Ucore, S, Vcore =  _FaustCorePy.FaustAlgoGenGivensCplxDbl.svdtj(*svdtj_args)
        elif(isinstance(M, csr_matrix)):
            Ucore, S, Vcore =  _FaustCorePy.FaustAlgoGenGivensCplxDbl.svdtj_sparse(*svdtj_args)
        else:
            raise ValueError("invalid type for M (first argument): only np.ndarray "
                             "or scipy.sparse.csr_matrix are supported.")

    U = Faust(core_obj=Ucore)
    V = Faust(core_obj=Vcore)
    return U, S, V

def pinvtj(M, nGivens=None, tol=0, err_period=100, relerr=True,
           nGivens_per_fac=None, enable_large_Faust=False, **kwargs):
    """
    Computes the pseudoinverse of M using svdtj.

    Args:
        M: the matrix to compute the pseudoinverse.
        nGivens: cf. svdtj
        tol:  cf. svdtj (here the error is computed on U.H S^+ V).
        err_period: cf. svdtj
        relerr: cf. svdtj
        nGivens_per_fac: cf. svdtj
        enable_large_Faust: see svdtj.

    Returns:
        The tuple V,Sp,Uh: such that V*numpy.diag(Sp)*Uh is the approximate of M^+.
        - (np.array vector) Sp the inverses of the min(m, n) nonzero singular values
        in ascending order. Note however that zeros might occur if M is rank r < min(*M.shape).
        - (Faust objects) V, Uh orthonormal Fausts.


    Example:
        >>> from pyfaust.fact import pinvtj
        >>> from scipy.sparse import spdiags
        >>> import numpy as np
        >>> from numpy.random import rand, seed
        >>> seed(42)
        >>> M = np.random.rand(128, 64)
        >>> V, Sp, Uh = pinvtj(M, tol=1e-3)
        >>> scipy_Sp = spdiags(Sp, [0], M.shape[1], M.shape[0])
        >>> np.linalg.norm(V @ scipy_Sp @ Uh @ M - np.eye(64, 64)) / np.linalg.norm(np.eye(64, 64))
        0.00012007583711484639

    """
    from os import environ
    environ['PINVTJ_ERR'] = '1'
    U, S, V = svdtj(M, nGivens=nGivens, tol=tol, err_period=err_period,
          relerr=relerr, nGivens_per_fac=nGivens_per_fac,
          enable_large_Faust=enable_large_Faust, **kwargs)
    environ['PINVTJ_ERR'] = '0'
    Sp = 1/S
    Sp[Sp == np.inf] = 0
    return V, Sp, U.H

def eigtj(M, nGivens=None, tol=0, err_period=100, order='ascend', relerr=True,
          nGivens_per_fac=None, verbosity=0, enable_large_Faust=False):
    """
    Performs an approximate eigendecomposition of M and returns the eigenvalues in W along with the corresponding normalized right eigenvectors (as the columns of the Faust object V).

    The output is such that V*numpy.diag(W)*V.H approximates M. V is a product
    of Givens rotations obtained by truncating the Jacobi algorithm.

    The trade-off between accuracy and complexity of V can be set through the
    parameters nGivens and tol that define the targeted number of Givens rotations
    and targeted error.

    Args:
        M: (numpy.ndarray or csr_matrix) the matrix to diagonalize. Must be
        real and symmetric, or complex hermitian. Can be in dense or sparse
        format. The dtype must be float32, float64
        or complex128 (the dtype might have a large impact on performance).
        nGivens: (int) targeted number of Givens rotations (this argument is optional
        only if tol is set).
        tol: (float) the tolerance error at which the algorithm stops. The
        default value is zero so that stopping is based on reaching the
        targeted nGivens (this argument is optional only if nGivens is set).
        err_period: (int) it defines the period, in number of factors of V
        the error is compared to tol (reducing the period spares
        some factors but increases slightly the computational cost because the error
        is computed more often).
        order: order of eigenvalues, possible choices are ‘ascend,
        'descend' or 'undef' (to avoid a sorting operation and save some time).
        nGivens_per_fac: (int) targeted number of Givens rotations per factor
        of V. Must be an integer between 1 to floor(M.shape[0]/2) (the default
        value).
        relErr: (bool) the type of error used as stopping criterion.  True
        for the relative error norm(V*D*V'-M, 'fro')/norm(M, 'fro'), False
        for the absolute error norm(V*D*V'-M, 'fro').
        verbosity: (int) the level of verbosity. The greater the value the more
        info is displayed. It can be helpful to understand for example why the
        algorithm stopped before reaching the tol error or the number of Givens
        (nGivens).
        enable_large_Faust: (bool)  if true, it allows to compute a transform
        that doesn't worth it regarding its complexity compared to the matrix
        M. Otherwise by default, an exception is raised before the algorithm starts.


    Returns:
        The tuple (W,V):
           - W (numpy.ndarray) the vector of the approximate eigenvalues of M
           (in ascending order by default).
           - V the Faust object representing the approximate eigenvector
            transform. The column V[:, i] is the eigenvector
            corresponding to the eigenvalue W[i].<br/>

    Remarks:
        - When  ‘nGivens’ and ‘tol’ are used simultaneously, the number of Givens
        rotations in V may be smaller than specified by ‘nGivens’ if the error
        criterion is met first, and the achieved error may be larger than specified
        if ‘nGivens’ is reached first during the iterations of the truncated Jacobi
        algorithm.
        - When nGivens_per_fac > 1, all factors have exactly
        nGivens_per_fac except the leftmost one which may have fewer if the
        total number of Givens rotations is not a multiple of
        nGivens_per_fact

    References:
    [1]   Le Magoarou L., Gribonval R. and Tremblay N., "Approximate fast
    graph Fourier transforms via multi-layer sparse approximations",
    IEEE Transactions on Signal and Information Processing
    over Networks 2018, 4(2), pp 407-420
    <https://hal.inria.fr/hal-01416110>

    Example:
        >>> import numpy as np
        >>> from pyfaust.fact import eigtj
        >>> from scipy.io import loadmat
        >>> from os.path import sep
        >>> from pyfaust.demo import get_data_dirpath
        >>> from numpy.linalg import norm
        >>> # get a graph Laplacian to diagonalize
        >>> demo_path = sep.join((get_data_dirpath(),'Laplacian_256_community.mat'))
        >>> data_dict = loadmat(demo_path)
        >>> Lap = data_dict['Lap'].astype('float64')
        >>> Dhat, Uhat = eigtj(Lap, nGivens=Lap.shape[0]*100, enable_large_Faust=True)
        >>> # Uhat is the Fourier matrix/eigenvectors approximation as a Faust
        >>> # (200 factors)
        >>> # Dhat the eigenvalues diagonal matrix approx.
        >>> print("err: ", norm(Lap-Uhat@np.diag(Dhat)@Uhat.H)/norm(Lap)) # about 6.5e-3
        err:  0.006494660956003688
        >>> Dhat2, Uhat2 = eigtj(Lap, tol=0.01)
        >>> assert(norm(Lap-Uhat2@np.diag(Dhat2)@Uhat2.H)/norm(Lap) < .011)
        >>> # and then asking for an absolute error
        >>> Dhat3, Uhat3 = eigtj(Lap, tol=0.1, relerr=False)
        >>> assert(norm(Lap-Uhat3@np.diag(Dhat3)@Uhat3.H) < .11)
        >>> # now recompute Uhat2, Dhat2 but asking a descending order of eigenvalues
        >>> Dhat4, Uhat4 = eigtj(Lap, tol=0.01, order='descend')
        >>> assert((Dhat4[::-1] == Dhat2[::]).all())
        >>> # and now with no sort
        >>> Dhat5, Uhat5 = eigtj(Lap, tol=0.01, order='undef')
        >>> assert((np.sort(Dhat5) == Dhat2).all())

    See also:
        svdtj
    """
    if not isinstance(M, scipy.sparse.csr_matrix):
        M = _check_fact_mat('eigtj()', M)

    M_dtype = _sanitize_dtype(M.dtype)

    eigtj_args = [M, nGivens, tol, relerr,
                  nGivens_per_fac, verbosity, order,
                  enable_large_Faust,
                  err_period]
    if M_dtype == 'float32':
        D, core_obj = _FaustCorePy.FaustAlgoGenGivensFlt.eigtj(*eigtj_args)

    elif M_dtype == 'float64':
        D, core_obj = _FaustCorePy.FaustAlgoGenGivensDbl.eigtj(*eigtj_args)
    else: # M_dtype == 'complex'
        D, core_obj = _FaustCorePy.FaustAlgoGenGivensCplxDbl.eigtj(*eigtj_args)
    return D, Faust(core_obj=core_obj)

def _check_fact_mat(funcname, M, is_real=None):
    if not isinstance(M, np.ndarray):
        raise Exception(funcname+" 1st argument must be a numpy ndarray.")

    use_is_real = isinstance(is_real, (list, np.ndarray))
    M_dtype = _sanitize_dtype(M.dtype)
    if M_dtype in ['float64', 'float32']:
        if use_is_real:
            is_real[0] = True
    else: #M_dtype == 'complex':
        if use_is_real:
            is_real[0] = False

    ndim_M=M.ndim;

    if (ndim_M > 2) or (ndim_M < 1):
        raise ValueError(funcname+': input matrix/array invalid number of dimensions')

    if not M.flags['F_CONTIGUOUS']:
#        raise ValueError(funcname+' input array must be Fortran contiguous (Colmajor)')
        M = np.asfortranarray(M)
    return M



# experimental block start
from numpy.linalg import norm
def hierarchical_py(A, J, N, res_proxs, fac_proxs, is_update_way_R2L=False,
                    is_fact_side_left=False, compute_2norm_on_arrays=False,
                    norm2_max_iter=100, norm2_threshold=1e-6, use_csr=True,
                    dev='cpu'):
    """
    PALM4MSA-hierarchical factorization.

    Args:
        A: (np.ndarray) the matrix to factorize. The dtype can be float32 or float64 (it might change the performance).
        J: number of factors.
        N: number of iterations of PALM4MSA calls.
        res_proxs: the residual factor proximity operators.
        fac_proxs: the main factor proximity operators.
        is_update_way_R2L: True to update factors from the right to left in PALM4MSA, False for the inverse order.
        is_fact_side_left: True to slit the left-most factor in the hierarchical factorization, False to split right-most factor.
        compute_2norm_on_arrays: True to compute left/right factors of the current one updated as arrays instead of computing the Faust norms.
        norm2_max_iter: the maximum number of iterations of the 2-norm computation algorithm.
        norm2_threshold: the threshold stopping criterion of the 2-norm computation algorithm.
        use_csr: True to update factors as CSR matrix in PALM4MSA.
        dev: 'cpu' or 'gpu', to use CPU or GPU Faust in the algorithm.

    Returns:
        the Faust resulting from the factorization of A.
    """
    S = Faust([A], dev=dev, dtype=A.dtype)
    l2_ = 1
    compute_2norm_on_arrays_ = compute_2norm_on_arrays
    for i in range(J-1):
        if(isinstance(compute_2norm_on_arrays, list)):
            compute_2norm_on_arrays_ = compute_2norm_on_arrays[i]
        print("hierarchical_py factor", i+1)
        if(is_fact_side_left):
            Si = S.factors(0)
            split_proxs = [res_proxs[i], fac_proxs[i]]
        else:
            Si = S.factors(i)
            split_proxs = [fac_proxs[i], res_proxs[i]]
        Si, l_ = palm4msa_py(Si, 2, N, split_proxs, is_update_way_R2L,
                             S='zero_and_ids', _lambda=1,
                             compute_2norm_on_arrays=compute_2norm_on_arrays_,
                             norm2_max_iter=norm2_max_iter,
                             norm2_threshold=norm2_threshold, use_csr=use_csr,
                             dev=dev)
        l2_ *= l_
        if i > 1:
            if is_fact_side_left:
                S = Si@S.right(1)
            else:
                S = S.left(i-1)@Si
        elif i > 0:
            if is_fact_side_left:
                S = Si@Faust(S.right(1), dev=dev)
            else:
                S = Faust(S.left(i-1), dev=dev)@Si
        else:  # i == 0
            S = Si
        S = S*1/l_
        if is_fact_side_left:
            fp = [*fac_proxs[0:i+1]]
            fp = list(reversed(fp))
            n_proxs = [res_proxs[i], *fp]
        else:
            n_proxs = [*fac_proxs[0:i+1],
                       res_proxs[i]]
            S, l2_ = palm4msa_py(A, S.numfactors(), N, n_proxs,
                                 is_update_way_R2L, S=S, _lambda=l2_,
                                 compute_2norm_on_arrays=compute_2norm_on_arrays_,
                                 norm2_max_iter=norm2_max_iter,
                                 norm2_threshold=norm2_threshold, use_csr=use_csr,
                                 dev=dev)
        S = S*1/l2_
    S = S*l2_
    return S


def palm4msa_py(A, J, N, proxs, is_update_way_R2L=False, S=None, _lambda=1,
                compute_2norm_on_arrays=False, norm2_max_iter=100,
                norm2_threshold=1e-6, use_csr=True, dev='cpu'):
    """
    PALM4MSA factorization.

    Args:
        A: (np.ndarray) the matrix to factorize. The dtype can be float32 or float64 (it might change the performance).
        J: number of factors.
        N: number of iterations of PALM4MSA.
        proxs: the factor proximity operators.
        is_update_way_R2L: True to update factors from the right to left in PALM4MSA, False for the inverse order.
        S: The Faust (sequence of factors) to initialize the PALM4MSA. By
        default, the first factor to be updated is zero, the other are the
        identity/eye matrix.
        compute_2norm_on_arrays: True to compute left/right factors of the current one updated as arrays instead of computing the Faust norms.
        norm2_max_iter: the maximum number of iterations of the 2-norm computation algorithm.
        norm2_threshold: the threshold stopping criterion of the 2-norm computation algorithm.
        use_csr: True to update factors as CSR matrix in PALM4MSA.
        dev: 'cpu' or 'gpu', to use CPU or GPU Faust in the algorithm.

    Returns:
        the Faust resulting from the factorization of A and its scale factor
        lambda (note that lambda is already applied to the output Faust. It is
        returned only for information, which is useful in hierarchical_py).
    """

    dims = [(prox.constraint._num_rows, prox.constraint._num_cols) for prox in
            proxs]
    A_H = A.T.conj()
    if not isinstance(A_H, np.ndarray):
        A_H = A_H.toarray()
    if S == 'zero_and_ids' or S is None:
        # start Faust, identity factors and one zero
        if is_update_way_R2L:
            S = Faust([np.eye(dims[i][0], dims[i][1], dtype=A.dtype) for i in
                       range(J-1)]+[np.zeros((dims[J-1][0], dims[J-1][1]), dtype=A.dtype)],
                      dev=dev)
        else:
            S = Faust([np.zeros((dims[0][0], dims[0][1]), dtype=A.dtype)] +
                      [np.eye(dims[i+1][0],
                              dims[i+1][1], dtype=A.dtype)
                       for i in
                       range(J-1)],
                      dev=dev)
    lipschitz_multiplicator = 1.001
    for i in range(N):
        if is_update_way_R2L:
            iter_ = reversed(range(J))
        else:
            iter_ = range(J)
        for j in iter_:
            if j == 0:
                S_j = S.factors(j)
                R = S.right(j+1)
                L = np.eye(S_j.shape[0], S_j.shape[0], dtype=A.dtype)
            elif(j == S.numfactors()-1):
                S_j = S.factors(j)
                R = np.eye(S_j.shape[1], S_j.shape[1], dtype=A.dtype)
                L = S.left(j-1)
            else:
                L = S.left(j-1)
                R = S.right(j+1)
                S_j = S.factors(j)
            if not pyfaust.isFaust(L):
                L = Faust(L, dev=dev, dtype=A.dtype)
            if not pyfaust.isFaust(R):
                R = Faust(R, dev=dev, dtype=A.dtype)
            if compute_2norm_on_arrays:
                c = \
                        lipschitz_multiplicator*_lambda**2*norm(R.toarray(), 2)**2 * \
                        norm(L.toarray(), 2)**2

            else:
                c = \
                        lipschitz_multiplicator*_lambda**2*R.norm(2, max_num_its=norm2_max_iter,
                                                                  threshold=norm2_threshold)**2 * \
                        L.norm(2, max_num_its=norm2_max_iter, threshold=norm2_threshold)**2
            if np.isnan(c) or c == 0:
                raise Exception("Failed to compute c (inverse of descent step),"
                                "it could be because of the Faust 2-norm error,"
                                "try option compute_2norm_on_arrays=True")
            if not isinstance(S_j, np.ndarray):
                S_j = S_j.toarray()
            D = S_j-_lambda*(L.H@(_lambda*L@(S_j@R)-A)@R.H)*1/c
            if(not isinstance(D, np.ndarray)):
                D = D.toarray()
            S_j = proxs[j](D)
            if use_csr:
                S_j = csr_matrix(S_j)
            if S.numfactors() > 2 and j > 0 and j < S.numfactors()-1:
                S = L@Faust(S_j, dev=dev, dtype=A.dtype)@R
            elif j == 0:
                S = Faust(S_j, dev=dev, dtype=A.dtype)@R
            else:
                S = L@Faust(S_j, dev=dev, dtype=A.dtype)
        _lambda = np.trace(A_H@S).real/S.norm()**2
    S = _lambda*S
    return S, _lambda
# experimental block end

# experimental block start
def hierarchical2020(M, nites, constraints, is_update_way_R2L,
                     is_fact_side_left, factor_format, packing_RL, norm2_threshold,
                     norm2_max_iter):
    factor_format = ParamsFact.factor_format_str2int(factor_format)
    core_obj,_lambda = \
            _FaustCorePy.FaustAlgoGenDbl.hierarchical2020(M, nites, constraints, is_update_way_R2L,
                                                    is_fact_side_left,
                                                    factor_format,
                                                    packing_RL, norm2_threshold,
                                                    norm2_max_iter)
    F = Faust(core_obj=core_obj)
    return F, _lambda

# experimental block end

def palm4msa(M, p, ret_lambda=False, backend=2016, on_gpu=False):
    """
    Factorizes the matrix M with Palm4MSA algorithm using the parameters set in p.

    Args:
        M: the numpy array to factorize. The dtype must be float32, float64
        or complex128 (the dtype might have a large impact on performance).
        p: the The pyfaust.factparams.ParamsPalm4MSA instance to define the algorithm parameters.
        ret_lambda: set to True to ask the function to return the scale factor (False by default).
        backend: the C++ implementation to use (default to 2016, 2020 backend
        should be faster for most of the factorizations).
        on_gpu: if True the GPU implementation is executed (this option applies only to 2020 backend).

    NOTE: If backend=2020 and independently to the StoppingCriterion defined in p,
    it is possible to stop the algorithm manually at any iteration by the key
    combination CTRL-C.
    The last Faust computed in the factorization process will be returned.
    A typical use case is when the verbose mode is enabled and you see that the error
    doesn't change anymore or only slightly, you might stop iterations by typing CTRL-C.

    Returns:
        The Faust object resulting of the factorization.
        if ret_lambda == True then the function returns a tuple (Faust, lambda).

    Examples:
        >>> from pyfaust.fact import palm4msa
        >>> from pyfaust.factparams import ParamsPalm4MSA, ConstraintList, StoppingCriterion
        >>> import numpy as np
        >>> M = np.random.rand(500, 32)
        >>> cons = ConstraintList('splin', 5, 500, 32, 'normcol', 1.0, 32, 32)
        >>> # or alternatively using pyfaust.proj
        >>> # from pyfaust.proj import splin, normcol
        >>> # cons = [ splin((500,32), 5), normcol((32,32), 1.0)]
        >>> stop_crit = StoppingCriterion(num_its=200)
        >>> param = ParamsPalm4MSA(cons, stop_crit)
        >>> F = palm4msa(M, param)
        >>> F
        Faust size 500x32, density 0.22025, nnz_sum 3524, 2 factor(s):
        - FACTOR 0 (double) SPARSE, size 500x32, density 0.15625, nnz 2500
        - FACTOR 1 (double) DENSE, size 32x32, density 1, nnz 1024

    See also pyfaust.factparams.ParamsPalm4msaWHT to factorize a Hadamard matrix using the SKPERM projector.
    """
    if(not isinstance(p, pyfaust.factparams.ParamsPalm4MSA)):
        raise TypeError("p must be a ParamsPalm4MSA object.")
    is_real = np.empty((1,))
    M = _check_fact_mat('palm4msa()', M, is_real)
    p.are_constraints_consistent(M)
    if(backend == 2016):
        if on_gpu: raise ValueError("on_gpu applies only on 2020 backend.")
        if is_real:
            is_float = M.dtype == 'float32'
            if is_float:
                core_obj, _lambda = _FaustCorePy.FaustAlgoGenFlt.fact_palm4msa(M, p)
            else:
                core_obj, _lambda = _FaustCorePy.FaustAlgoGenDbl.fact_palm4msa(M, p)
        else:
            core_obj, _lambda = _FaustCorePy.FaustAlgoGenCplxDbl.fact_palm4msa(M, p)
    elif(backend == 2020):
        if is_real:
            if M.dtype == 'float64':
                core_obj, _lambda = _FaustCorePy.FaustAlgoGenDbl.palm4msa2020(M, p, on_gpu)
            else: # M.dtype == 'float32': ensured by _check_fact_mat
                core_obj, _lambda = _FaustCorePy.FaustAlgoGenFlt.palm4msa2020(M, p, on_gpu)
        else:
            if M.dtype == np.complex128 and p.factor_format != 'dense':
                p.factor_format = 'dense'
                warnings.warn("forcing the factor_format parameter to 'dense'")
            core_obj, _lambda = _FaustCorePy.FaustAlgoGenCplxDbl.palm4msa2020(M, p)
    else:
        raise ValueError("Unknown backend (only 2016 and 2020 are available).")
    F = Faust(core_obj=core_obj)
    if(ret_lambda):
        return F, _lambda
    else:
        return F

def palm4msa_mhtp(M, palm4msa_p, mhtp_p, ret_lambda=False, on_gpu=False):
    """
    Runs the MHTP-PALM4MSA algorithm to factorize the matrix M.

    MHTP stands for Multilinear Hard Tresholding Pursuit.
    This is a generalization of the Bilinear HTP algorithm describe in [1].

    [1] Quoc-Tung Le, Rémi Gribonval. Structured Support Exploration For
    Multilayer Sparse Matrix Factorization. ICASSP 2021 - IEEE International Conference on Acoustics,
    Speech and Signal Processing,
    Jun 2021, Toronto, Ontario, Canada. pp.1-5. <a
    href="https://hal.inria.fr/hal-03132013/document">hal-03132013</a>



    Args:
        M: the numpy array to factorize. The dtype must be float32, float64
        or complex128 (the dtype might have a large impact on performance).
        palm4msa_p: the The pyfaust.factparams.ParamsPalm4MSA instance to
        define the PALM4MSA algorithm parameters.
        mhtp_p: the pyfaust.factparams.MHTPParams instance to define the MHTP algorithm parameters.
        ret_lambda: set to True to ask the function to return the scale factor (False by default).
        on_gpu: if True the GPU implementation is executed.

    Returns:
        The Faust object resulting of the factorization.
        if ret_lambda == True then the function returns a tuple (Faust, lambda).

    Example:
		>>> from pyfaust.fact import palm4msa_mhtp
		>>> from pyfaust.factparams import ParamsPalm4MSA, StoppingCriterion, MHTPParams
		>>> import numpy as np
		>>> from pyfaust.proj import splin, normcol
        >>> np.random.seed(42) # just for reproducibility
		>>> M = np.random.rand(500, 32)
		>>> cons = [ splin((500,32), 5), normcol((32,32), 1.0)]
		>>> stop_crit = StoppingCriterion(num_its=200)
		>>> param = ParamsPalm4MSA(cons, stop_crit)
		>>> # MHTP will run every 100 iterations of PALM4MSA (that is 2 times) for 50 iterations on each factor
		>>> mhtp_param = MHTPParams(num_its=50, palm4msa_period=100)
		>>> G = palm4msa_mhtp(M, param, mhtp_param)
        >>> G
        Faust size 500x32, density 0.22025, nnz_sum 3524, 2 factor(s):
        - FACTOR 0 (double) SPARSE, size 500x32, density 0.15625, nnz 2500
        - FACTOR 1 (double) DENSE, size 32x32, density 1, nnz 1024

    <b/> See also pyfaust.factparams.MHTPParams
    <b/> See also pyfaust.fact.palm4msa
    <b/> See also pyfaust.fact.hierarchical_mhtp
    """
    palm4msa_p.use_MHTP = mhtp_p
    return palm4msa(M, palm4msa_p, ret_lambda=ret_lambda, backend=2020, on_gpu=on_gpu)

def hierarchical_mhtp(M, hierar_p, mhtp_p, ret_lambda=False, ret_params=False,
                      on_gpu=False):
    """
    Runs the MHTP-PALM4MSA hierarchical factorization algorithm on the matrix M.

    This algorithm uses the MHTP-PALM4MSA (pyfaust.fact.palm4msa_mhtp) instead of only PALM4MSA as
    pyfaust.fact.hierarchical.

    Args:
        M: the numpy array to factorize. The dtype must be float32, float64
        or complex128 (the dtype might have a large impact on performance).
        p: is a set of hierarchical factorization parameters. See pyfaust.fact.hierarchical.
        mhtp_p: the pyfaust.factparams.MHTPParams instance to define the MHTP algorithm parameters.
        on_gpu: if True the GPU implementation is executed.
        ret_lambda: set to True to ask the function to return the scale factor (False by default).
        ret_params: set to True to ask the function to return the
        ParamsHierarchical instance used (False by default).


    Example:
		>>> from pyfaust.fact import hierarchical_mhtp
		>>> from pyfaust.factparams import ParamsHierarchical, StoppingCriterion
		>>> from pyfaust.factparams import MHTPParams
		>>> from pyfaust.proj import sp, normcol, splin
		>>> import numpy as np
		>>> M = np.random.rand(500, 32)
		>>> fact_cons = [splin((500, 32), 5), sp((32,32), 96), sp((32,32), 96)]
		>>> res_cons = [normcol((32,32), 1), sp((32,32), 666), sp((32,32), 333)]
		>>> stop_crit1 = StoppingCriterion(num_its=200)
		>>> stop_crit2 = StoppingCriterion(num_its=200)
		>>> # 50 iterations of MHTP will run every 100 iterations of PALM4MSA (each time PALM4MSA is called by the hierarchical algorithm)
		>>> mhtp_param = MHTPParams(num_its=50, palm4msa_period=100)
		>>> param = ParamsHierarchical(fact_cons, res_cons, stop_crit1, stop_crit2)
		>>> param.is_verbose = True
		>>> F = hierarchical_mhtp(M, param, mhtp_param)

    <b/> See also pyfaust.fact.hierarchical, pyfaust.fact.palm4msa_mhtp
    """
    hierar_p.use_MHTP = mhtp_p
    return hierarchical(M, hierar_p, ret_lambda=False, ret_params=False, backend=2020,
                        on_gpu=False)


# experimental block start
def _palm4msa_fgft(Lap, p, ret_lambda=False):
    """
    """
    if(not isinstance(p, pyfaust.factparams.ParamsPalm4MSAFGFT)):
        raise TypeError("p must be a ParamsPalm4MSAFGFT object.")
    is_real = np.empty((1,))
    Lap = _check_fact_mat('_palm4msa_fgft()', Lap, is_real)
    if((Lap.T != Lap).any() or Lap.shape[0] != Lap.shape[1]):
        raise ValueError("Laplacian matrix must be symmetric.")
    if(not p.is_mat_consistent(Lap)):
        raise ValueError("M's number of columns must be consistent with "
                         "the last residual constraint defined in p. "
                         "Likewise its number of rows must be consistent "
                         "with the first factor constraint defined in p.")
    if is_real:
        core_obj, _lambda, D = _FaustCorePy.FaustAlgoGenDbl.fact_palm4msa_fft(Lap, p)
    else:
        core_obj, _lambda, D = _FaustCorePy.FaustAlgoGenCplxDbl.fact_palm4msa_fft(Lap, p)

    F = Faust(core_obj=core_obj)
    if(ret_lambda):
        return F, D, _lambda
    else:
        return F, D

# experimental block end

def hierarchical(M, p, ret_lambda=False, ret_params=False, backend=2016,
                 on_gpu=False):
    """
    Factorizes the matrix M using the hierarchical PALM4MSA algorithm according to the parameters set in p.

    @note This function has its shorthand pyfaust.faust_fact(). For
    convenience you might use it like this:<br/>
            <code>
            from pyfaust import *;
            F = faust_fact(M, p) # equiv. to hierarchical(M, p)
            </code>

    Basically, the hierarchical factorization works in a sequence of splits in two of the
    last factor.
    The first split consists to split the matrix M in two, when
    p.is_fact_side_left is False (which is the case by default), the
    factorization works toward right, so M is factorized as follows:
        \f[M \approx S_1 R_1\f]
    We call \f$S_1\f$ the main factor and \f$R_1\f$ the residual factor.
    On step 2, \f$R_1\f$ is split in two such as \f$R_1 \approx S_2 R_2\f$, which gives:
        \f[M \approx S_1 S_2 R_2 \f]
    And so on until the algorithm reaches the targeted number of factors:
        \f[M \approx S_1 S_2 ... S_{N-1} R_N \f]
    If p.is_fact_side_left is False, the residual is factorized toward left, so
    it gives rather :
        \f[M \approx R_1 S_1 \\ \vdots \\ M \approx R_2 S_2 S_1 \\ \vdots \\ M \approx R_N S_{N-1}  S_2 S_1 ... \f]

    Args:
        M: the numpy array to factorize. The dtype must be float32, float64
        or complex128 (the dtype might have a large impact on performance).
        p: is a set of hierarchical factorization parameters. It might be a fully defined instance of parameters (pyfaust.factparams.ParamsHierarchical) or a simplified expression which designates a pre-defined parametrization:
            - 'hadamard' to use pre-defined parameters typically used to factorize a Hadamard matrix of order a power of two (see pyfaust.demo.hadamard).
            - ['rectmat', j, k, s] to use pre-defined parameters used for
            instance in factorization of the MEG matrix which is a rectangular
            matrix of size m*n such that m < n (see pyfaust.demo.bsl); j is the
            number of factors, k the sparsity of the main factor's columns, and
            s the sparsity of rows for all other factors except the residual factor
            (that is the first factor here because the factorization is made
            toward the left -- is_side_fact_left == true, cf.
            pyfaust.factparams.ParamsHierarchical and pyfaust.factparams.ParamsHierarchicalRectMat).
            <br/>The residual factor has a sparsity of P*rho^(num_facts-1). <br/> By default, rho == .8 and P = 1.4. It's possible to set custom values with for example p == ( ['rectmat', j, k, s], {'rho':.4, 'P':.7 }). <br/>The sparsity is here the number of non-zero elements.
        backend: the C++ implementation to use (default to 2016, 2020 backend
        should be faster for most of the factorizations).
        on_gpu: if True the GPU implementation is executed (this option applies only to 2020 backend).

        ret_lambda: set to True to ask the function to return the scale factor (False by default).
        ret_params: set to True to ask the function to return the
        ParamsHierarchical instance used (False by default).
        It is useful for consulting what precisely means the
        simplified parametrizations used to generate a
        ParamsHierarchical instance and possibly adjust its attributes to factorize again.

    NOTE: If backend=2020 and independently to the StoppingCriterion-s defined in p,
    it is possible to stop any internal call to PALM4MSA manually at any iteration
    by the key combination CTRL-C.
    The last Faust computed in the PALM4MSA instance will be used to continue
    the hierarchical factorization.
    A typical use case is when the verbose mode is enabled and you see that the error
    doesn't change anymore or only slightly, you might stop iterations by typing CTRL-C.

    Returns:
        F the Faust object result of the factorization: Faust\f$([S_1, S_2, ...
        ,S_{N-1}, R_N]])\f$
        if p.is_fact_side_left == False, Faust\f$([R_N, S_{N-1}, ... , S_2,
        S_1])\f$
        otherwise.<br/>
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
        >>> # or alternatively using pyfaust.proj
        >>> # from pyfaust.proj import *
        >>> # res_cons = [normcol((32,32), 1), sp((32,32), 666), sp((32,32), 333)]
        >>> stop_crit1 = StoppingCriterion(num_its=200)
        >>> stop_crit2 = StoppingCriterion(num_its=200)
        >>> param = ParamsHierarchical(fact_cons, res_cons, stop_crit1, stop_crit2)
        >>> F = hierarchical(M, param)
        >>> F
        Faust size 500x32, density 0.189063, nnz_sum 3025, 4 factor(s):
        - FACTOR 0 (double) SPARSE, size 500x32, density 0.15625, nnz 2500
        - FACTOR 1 (double) SPARSE, size 32x32, density 0.09375, nnz 96
        - FACTOR 2 (double) SPARSE, size 32x32, density 0.09375, nnz 96
        - FACTOR 3 (double) SPARSE, size 32x32, density 0.325195, nnz 333

        <b>2. Simplified Parameters for Hadamard Factorization</b>
        >>> from pyfaust import wht
        >>> from pyfaust.fact import hierarchical
        >>> from numpy.linalg import norm
        >>> # generate a Hadamard Faust of size 32x32
        >>> FH = wht(32)
        >>> H = FH.toarray() # the full matrix version
        >>> # factorize it
        >>> FH2 = hierarchical(H, 'hadamard');
        >>> # test the relative error
        >>> (FH-FH2).norm('fro')/FH.norm('fro') # the result is about 1e-16, the factorization is accurate
        1.6883057247219393e-16
        >>> FH
        Faust size 32x32, density 0.3125, nnz_sum 320, 5 factor(s):
        - FACTOR 0 (double) SPARSE, size 32x32, density 0.0625, nnz 64
        - FACTOR 1 (double) SPARSE, size 32x32, density 0.0625, nnz 64
        - FACTOR 2 (double) SPARSE, size 32x32, density 0.0625, nnz 64
        - FACTOR 3 (double) SPARSE, size 32x32, density 0.0625, nnz 64
        - FACTOR 4 (double) SPARSE, size 32x32, density 0.0625, nnz 64
        >>> FH2
        Faust size 32x32, density 0.3125, nnz_sum 320, 5 factor(s):
        - FACTOR 0 (double) SPARSE, size 32x32, density 0.0625, nnz 64
        - FACTOR 1 (double) SPARSE, size 32x32, density 0.0625, nnz 64
        - FACTOR 2 (double) SPARSE, size 32x32, density 0.0625, nnz 64
        - FACTOR 3 (double) SPARSE, size 32x32, density 0.0625, nnz 64
        - FACTOR 4 (double) SPARSE, size 32x32, density 0.0625, nnz 64


       <br/>
       See also pyfaust.factparams.ParamsHierarchicalWHT
       <br/>

       <b> 3. Simplified Parameters for a Rectangular Matrix Factorization
       (the BSL demo MEG matrix)</b>
       >>> from pyfaust import *
       >>> from pyfaust.fact import hierarchical
       >>> from scipy.io import loadmat
       >>> from pyfaust.demo import get_data_dirpath
       >>> import numpy as np
       >>> d = loadmat(get_data_dirpath()+'/matrix_MEG.mat')
       >>> MEG = d['matrix'].T
       >>> num_facts = 9
       >>> k = 10
       >>> s = 8
       >>> MEG16 = hierarchical(MEG, ['rectmat', num_facts, k, s])
       >>> MEG16
       Faust size 204x8193, density 0.0631655, nnz_sum 105573, 9 factor(s):
       - FACTOR 0 (double) SPARSE, size 204x204, density 0.293613, nnz 12219
       - FACTOR 1 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
       - FACTOR 2 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
       - FACTOR 3 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
       - FACTOR 4 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
       - FACTOR 5 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
       - FACTOR 6 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
       - FACTOR 7 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
       - FACTOR 8 (double) SPARSE, size 204x8193, density 0.0490196, nnz 81930
       >>> # verify the constraint k == 10, on column 4
       >>> np.count_nonzero(MEG16.factors(8)[:,4].toarray())
       10
       >>> # now verify the s constraint is respected on MEG16 factor 1
       >>> np.count_nonzero(MEG16.factors(1).toarray())/MEG16.shape[0]
       8.0


       <br/>
       See also pyfaust.factparams.ParamsHierarchicalRectMat
       <br/>

       <b>4. Simplified Parameters for the Discrete Fourier Transform Factorization</b>
       >>> from pyfaust import dft
       >>> from pyfaust.fact import hierarchical
       >>> from numpy.linalg import norm
       >>> # generate a DFT matrix of size 32x32
       >>> FDFT = dft(32)
       >>> DFT = FDFT.toarray()
       >>> # factorize it
       >>> FDFT2 = hierarchical(DFT, 'dft');
       >>> # test the relative error
       >>> (FDFT-FDFT2).norm('fro')/FDFT.norm('fro') # the result is about 1e-16, the factorization is accurate
       1.5867087530782683e-16
       >>> FDFT
       Faust size 32x32, density 0.34375, nnz_sum 352, 6 factor(s):
       - FACTOR 0 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
       - FACTOR 1 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
       - FACTOR 2 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
       - FACTOR 3 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
       - FACTOR 4 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
       - FACTOR 5 (complex) SPARSE, size 32x32, density 0.03125, nnz 32
       >>> FDFT2
       Faust size 32x32, density 0.34375, nnz_sum 352, 6 factor(s):
       - FACTOR 0 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
       - FACTOR 1 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
       - FACTOR 2 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
       - FACTOR 3 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
       - FACTOR 4 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
       - FACTOR 5 (complex) SPARSE, size 32x32, density 0.03125, nnz 32

       See also pyfaust.factparams.ParamsHierarchicalDFT
       <br/>

        <b>5. Simplified Parameters for Hadamard Factorization without residual
        constraints</b>

        This factorization parameterization is the same as the one shown in 2.
        except that there is no constraints at all on residual factors. See
        pyfaust.factparams.ParamsHierarchicalNoResCons and
        pyfaust.factparams.ParamsHierarchicalWHTNoResCons for more details.

        >>> from pyfaust import wht
        >>> from pyfaust.fact import hierarchical
        >>> from numpy.linalg import norm
        >>> # generate a Hadamard Faust of size 32x32
        >>> FH = wht(32)
        >>> H = FH.toarray() # the full matrix version
        >>> # factorize it
        >>> FH2 = hierarchical(H, 'hadamard_simple', backend=2020);
        >>> # test the relative error
        >>> (FH-FH2).norm('fro')/FH.norm('fro') # the result is about 1e-16, the factorization is accurate
        7.85046215906392e-16
        >>> FH2
        Faust size 32x32, density 0.3125, nnz_sum 320, 5 factor(s):
        - FACTOR 0 (double) SPARSE, size 32x32, density 0.0625, nnz 64
        - FACTOR 1 (double) SPARSE, size 32x32, density 0.0625, nnz 64
        - FACTOR 2 (double) SPARSE, size 32x32, density 0.0625, nnz 64
        - FACTOR 3 (double) SPARSE, size 32x32, density 0.0625, nnz 64
        - FACTOR 4 (double) SPARSE, size 32x32, density 0.0625, nnz 64

       <br/>
       See also pyfaust.factparams.ParamsHierarchicalWHTNoResCons
       <br/>

       <b> 6. Simplified Parameters for a Rectangular Matrix Factorization (the
       BSL demo MEG matrix) without residual constraints</b>

       The factorization parameterization shown here is the same as in 3.
       except that there is no constraint at all on residual factors. See
       pyfaust.factparams.ParamsHierarchicalNoResCons and
       pyfaust.factparams.ParamsHierarchicalRectMatNoResCons for more details.
       In the example below the MEG matrix is factorized according to the
       parameterization shown in 3. (aka "MEG") and on the other hand with
       the parameterization of interest here (aka "MEG_SIMPLE", with no
       residual constraints), the approximate accuracy is quite the same so we
       can conclude that on this case (as in 5.) removing residual constraints can
       not only simplify the parameterization of hierarchical PALM4MSA but also
       be as efficient.

        >>> from scipy.io import loadmat
        >>> from pyfaust import Faust
        >>> from pyfaust.fact import hierarchical
        >>> from pyfaust.demo import get_data_dirpath
        >>> MEG = loadmat(get_data_dirpath()+'/matrix_MEG.mat')['matrix'].T
        >>> F1 = hierarchical(MEG, ['MEG', 5, 10, 8], backend=2020)
        >>> F2 = hierarchical(MEG, ['MEG_SIMPLE', 5, 10, 8], backend=2020)
        >>> # compare the errors:
        >>> (F2 - MEG).norm() / Faust(MEG).norm()
        0.13033595653237676
        >>> (F1 - MEG).norm() / Faust(MEG).norm()
        0.12601709513741005

       <br/>
       <b>See also</b>  pyfaust.factparams.ParamsHierarchicalRectMatNoResCons
       <br/>

    <b/> See also pyfaust.factparams.ParamsHierarchicalRectMat, pyfaust.factparams.ParamsHierarchicalSquareMat, pyfaust.factparams.ParamsHierarchicalDFT
    <b/> See also pyfaust.factparams.ParamsHierarchicalRectMatNoResCons, pyfaust.factparams.ParamsHierarchicalWHTNoResCons

    [1] <a href="https://hal.archives-ouvertes.fr/hal-01167948">Le Magoarou L. and Gribonval R., “Flexible multi-layer
    sparse approximations of matrices and applications”, Journal of Selected
    Topics in Signal Processing, 2016.</a>

    """
    is_real = np.empty((1,))
    p, M = _prepare_hierarchical_fact(M, p, "hierarchical", ret_lambda,
                              ret_params, is_real=is_real)
    if(backend == 2016):
        if on_gpu: raise ValueError("on_gpu applies only on 2020 backend.")
        if is_real:
            if M.dtype == 'float64':
                core_obj,_lambda = _FaustCorePy.FaustAlgoGenDbl.fact_hierarchical(M, p)
            else:
                core_obj,_lambda = _FaustCorePy.FaustAlgoGenFlt.fact_hierarchical(M, p)
        else:
            core_obj,_lambda = _FaustCorePy.FaustAlgoGenCplxDbl.fact_hierarchical(M, p)
    elif(backend == 2020):
        if is_real:
            if M.dtype == 'float64':
                core_obj, _lambda = _FaustCorePy.FaustAlgoGenDbl.hierarchical2020(M, p,
                                                                                  on_gpu)
            else: #M.dtype = 'float32': # after _prepare_hierarchical_fact it is ensured
                core_obj, _lambda = _FaustCorePy.FaustAlgoGenFlt.hierarchical2020(M, p,
                                                                                  on_gpu)
        else:
#            if M.dtype == np.complex128 and p.factor_format != 'dense':
#                p.factor_format = 'dense'
#                warnings.warn("forcing the factor_format parameter to 'dense'")
            core_obj, _lambda = _FaustCorePy.FaustAlgoGenCplxDbl.hierarchical2020(M, p,
                                                                        on_gpu)
    else:
        raise ValueError("backend must be 2016 or 2020")
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
                               M_name='M', is_real=None):
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
    M = _check_fact_mat(callee_name+'()', M, is_real=is_real)
    if(not isinstance(ret_lambda, bool)):
        raise TypeError("ret_lambda must be a bool.")
    if(not isinstance(ret_params, bool)):
        raise TypeError("ret_params must be a bool.")
    p.are_constraints_consistent(M)
    return p, M

# experimental block start
def hierarchical_constends(M, p, A, B, ret_lambda=False, ret_params=False):
    """
    Tries to approximate M by \f$ A \prod_j S_j B \f$ using the hierarchical.

    It needs to add one additional residual constraint and also one factor
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
    # consts : factor constraints + residual constraints
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

# experimental block start
def fgft_palm(U, Lap, p, init_D=None, ret_lambda=False, ret_params=False):
    """
    Computes the FGFT for the Fourier matrix U (it should be the eigenvectors of the Laplacian Lap).


    NOTE: this algorithm is a variant of hierarchical.


    Example:
        from pyfaust.fact import fgft_palm
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
        err_Lap = norm(Uhat.toarray()*diag(Dhat)*Uhat.T.toarray()-Lap)/norm(Lap)
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
    from pyfaust.factparams import ParamsPalm4MSAFGFT

    is_real = np.empty((1,))
    p, M = _prepare_hierarchical_fact(U, p, "fgft_palm", ret_lambda,
                              ret_params, 'U', is_real=is_real)
    if(init_D.dtype != Lap.dtype or Lap.dtype != U.dtype or U.dtype != init_D.dtype):
        raise ValueError("All the numpy arrays must be of the same dtype.")

    init_D = ParamsPalm4MSAFGFT._set_init_D(init_D, U.shape[0])
    if is_real:
        core_obj, _lambda, D = _FaustCorePy.FaustAlgoGenDbl.fact_hierarchical_fft(U, Lap, p,
                                                                            init_D)
    else:
        core_obj, _lambda, D = _FaustCorePy.FaustAlgoGenCplxDbl.fact_hierarchical_fft(U, Lap, p,
                                                                                init_D)
    F = Faust(core_obj=core_obj)
    ret_list = [ F, D ]
    if(ret_lambda):
        ret_list += [ _lambda ]
    if(ret_params):
        ret_list += [ p ]
    return ret_list
# experimental block end


def butterfly(M, type="bbtree", perm=None, diag_opt=False):
    """
    Factorizes M according to a butterfly support and optionally a permutation using the algorithms described in [1].

    The result is a Faust F of the form BP where B has a butterfly structure
    and P is a permutation matrix determined by the optional parameter perm.

    Args:
        M: the numpy ndarray. The dtype must be float32, float64
        or complex128 (the dtype might have a large impact on performance). M
        must be square and its dimension must be a power of two.
        type: (str) the type of factorization 'right'ward, 'left'ward or
        'bbtree'. More precisely: if 'left' (resp. 'right') is used then at each stage of the
        factorization the most left factor (resp. the most right factor) is split in two.
        If 'bbtree' is used then the matrix is factorized according to a Balanced
        Binary Tree (which is faster as it allows parallelization).
        perm: five kinds of values are possible for this argument.
            1. perm is a list of column indices of the permutation matrix P which is such that
            the returned Faust is F = B@P where B is the Faust butterfly
            approximation of M@P.T.
            If the list of indices is not a valid permutation the behaviour
            is undefined (however an invalid size or an out of bound index raise
            an exception).
            2. perm is a list of lists of permutation column indices as defined
            in 1. In that case, all permutations passed to the function are
            used as explained in 1, each one producing a Faust, the best one
            (that is the best approximation of M in the Frobenius norm) is kept and returned by butterfly.
            3. perm is 'default_8', this is a particular case of 2. Eight
            default permutations are used. For the definition of those
            permutations please refer to [2].
            4. perm is 'bitrev': in that case the permutation is the
            bit-reversal permutation (cf. pyfaust.bitrev_perm).
            5. By default this argument is None, no permutation is used (this
            is equivalent to using the identity permutation matrix in 1).
        diag_opt: if True then the returned Faust is optimized using
        pyfaust.opt_butterfly_faust.

    Note: Below is an example of how to create a permutation scipy CSR matrix from a permutation list
    of indices (as defined by the perm argument) and conversely how to convert
    a permutation matrix to a list of indices of permutation.

    >>> from scipy.sparse import random, csr_matrix
    >>> from numpy.random import permutation
    >>> import numpy as np
    >>> np.random.seed(42)
    >>> I = permutation(8)  # random permutation as a list of indices
    >>> I
    array([1, 5, 0, 7, 2, 4, 3, 6])
    >>> n = len(I)
    >>> # convert a permutation as indices to a csr_matrix
    >>> P = csr_matrix((np.ones(n), (I, np.arange(n))))
    >>> P.toarray()
    array([[0., 0., 1., 0., 0., 0., 0., 0.],
           [1., 0., 0., 0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 1., 0., 0., 0.],
           [0., 0., 0., 0., 0., 0., 1., 0.],
           [0., 0., 0., 0., 0., 1., 0., 0.],
           [0., 1., 0., 0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 0., 0., 0., 1.],
           [0., 0., 0., 1., 0., 0., 0., 0.]])
    >>> # convert a permutation as a list of indices to a permutation matrix P as a csr_matrix
    >>> I_ = P.T.nonzero()[1]
    >>> I_
    array([1, 5, 0, 7, 2, 4, 3, 6], dtype=int32)
    >>> np.allclose(I_, I)
    True

    Returns:
        The Faust F which is an approximation of M according to a butterfly support.

    Examples:
        >>> import numpy as np
        >>> from random import randint
        >>> from pyfaust.fact import butterfly
        >>> from pyfaust import Faust, wht, dft
        >>> H = wht(8).toarray()
        >>> F = butterfly(H, type='bbtree')
        >>> # compute the error
        >>> (F-H).norm()/Faust(H).norm()
        1.2560739454502295e-15
        >>> # the same can be done with dft in place of wht
        >>> # all you need is to specify the bit-reversal permutation
        >>> # since the Discrete Fourier Transform is the product of a butterfly factors with this particular permutation
        >>> DFT = dft(8).toarray()
        >>> F = butterfly(DFT, type='bbtree', perm='bitrev')
        >>> # compute the error
        >>> (F-DFT).norm()/Faust(DFT).norm()
        1.1427230601405052e-15


        Use simple permutations:
        >>> import numpy as np
        >>> from random import randint
        >>> from pyfaust.fact import butterfly
        >>> M = np.random.rand(4, 4)
        >>> # without any permutation
        >>> F1 = butterfly(M, type='bbtree')
        >>> # which is equivalent to using the identity permutation
        >>> p = np.arange(0, 4)
        >>> F2 = butterfly(M, type='bbtree', perm=p)
        >>> # compute the relative diff
        >>> (F2-F1).norm()/F1.norm()
        0.0
        >>> # then try another permutation
        >>> p2 = [1, 0, 3, 2]
        >>> F3 = butterfly(M, type='bbtree', perm=p2)

        Use butterfly with a permutation defined by a list of indices J:
        >>> import numpy as np
        >>> from pyfaust.fact import butterfly
        >>> from pyfaust import Faust, wht, dft
        >>> H = wht(8).toarray()
        >>> J = np.arange(7, -1, -1)
        >>> F = butterfly(H, type='bbtree', perm=J)
        >>> # this is equivalent to passing a list containing a single permutation :
        >>> # F = butterfly(H, type='bbtree', perm=[J])
        # use butterfly with successive permutations J1 and J2
        # and keep the best approximation
        >>> J1 = J
        >>> from numpy.random import permutation
        >>> J2 = permutation(len(J)) # another random permutation
        >>> F = butterfly(H, type='bbtree', perm=[J1, J2])
        >>> # or to to use the 8 default permutations (keeping the best approximation resulting Faust)
        >>> F = butterfly(H, type='bbtree', perm='default_8')

    References:
        [1] Quoc-Tung Le, Léon Zheng, Elisa Riccietti, Rémi Gribonval. Fast
        learning of fast transforms, with guarantees. ICASSP, IEEE
        International Conference on Acoustics, Speech and Signal Processing,
        May 2022, Singapore, Singapore. (<a href="https://hal.inria.fr/hal-03438881">hal-03438881</a>) <br/>
		[2] T. Dao, A. Gu, M. Eichhorn, A. Rudra, and C. Re,
		“Learning Fast Algorithms for Linear Transforms Using
		Butterfly Factorizations,” in Proceedings of the 36th
		International Conference on Machine Learning. June
		2019, pp. 1517–1527, PMLR

     See also:
         pyfaust.wht, pyfaust.dft, pyfaust.rand_butterfly
    """
    def perm2indices(P):
            return P.T.nonzero()[1]
    is_real = np.empty((1,))
    M = _check_fact_mat('butterfly()', M, is_real)
    if isinstance(perm, str):
        if perm == 'bitrev':
            P = bitrev_perm(M.shape[1])
            return butterfly(M, type, perm=perm2indices(P), diag_opt=diag_opt)
        elif perm == 'default_8':
            # the three modified functions below were originally extracted from the 3 clause-BSD code hosted here: https://github.com/leonzheng2/butterfly
            # please look the header license here https://github.com/leonzheng2/butterfly/blob/main/src/utils.py
            def perm_type(i, type):
                """
                Type 0 is c in paper. Type 1 is b in paper. Type 2 is a in paper.
                :param i:
                :param type:
                :return:
                """
                size = 2 ** i
                if type == 0:
                    row_inds = np.hstack((np.arange(size//2), size//2 + np.arange(size//2)))
                    col_inds = np.hstack((np.arange(size//2), size - 1 - np.arange(size//2)))
                elif type == 1:
                    row_inds = np.hstack((size // 2 - 1 - np.arange(size//2), size // 2 + np.arange(size//2)))
                    col_inds = np.hstack((np.arange(size//2), size//2 + np.arange(size//2)))
                else:
                    row_inds = np.hstack((np.arange(size//2), size//2 + np.arange(size//2)))
                    col_inds = np.hstack((np.arange(size//2) * 2, np.arange(size//2) * 2 + 1))
                result = csr_matrix((np.ones(row_inds.size), (row_inds, col_inds)))
                return result

            def shared_logits_permutation(num_factors, choices):
                """
                :param num_factors:
                :param choices: array of three bool
                :return:
                """
                permutations = []
                for i in range(2, num_factors + 1):
                    block = seye(2 ** i)
                    if choices[0]:
                        block = block @ perm_type(i, 0)
                    if choices[1]:
                        block = block @ perm_type(i, 1)
                    if choices[2]:
                        block = block @ perm_type(i, 2)
                    perm = skron(seye(2 ** (num_factors - i)), block)
                    permutations.append(perm)
                return permutations

            def get_permutation_matrix(num_factors, perm_name):
                """
                :param num_factors:
                :param perm_name: str, 000, 001, ..., 111
                :return:
                """
                if perm_name.isnumeric():
                    choices = [int(char) for char in perm_name]
                    p_list = shared_logits_permutation(num_factors, choices)
                    p = csr_matrix(Faust(p_list).toarray()) # TODO: keep csr along the whole product
                else:
                    raise TypeError("perm_name must be numeric")
                return p

    #        print(list(perm2indices(get_permutation_matrix(int(np.log2(M.shape[0])),
    #                                               perm_name))+1 \
    #                        for perm_name in  ["000", "001", "010", "011", "100",
    #                                           "101", "110", "111"]))
            permutations = [perm2indices(get_permutation_matrix(int(np.log2(M.shape[0])),
                                                   perm_name)) \
                            for perm_name in  ["000", "001", "010", "011", "100", "101", "110", "111"]]
            return butterfly(M, type, permutations, diag_opt=diag_opt)
    elif isinstance(perm, (list, tuple)) and isinstance(perm[0], (list, tuple,
                                                                 np.ndarray)):
        # loop on each perm and keep the best approximation
        best_err = np.inf
        best_F = None
        for p in perm:
            # print(p)
            row_inds = np.arange(len(p))
            P = csr_matrix((np.ones(row_inds.size), (row_inds, p)))
            F = butterfly(M, type, p)
            # compute error
            error = np.linalg.norm(F.toarray()-M)/Faust(M).norm()
            # print(error)
            if error < best_err:
                best_err = error
                best_F = F
        return best_F
    args = (M, type, perm, not diag_opt)
    if is_real:
        is_float = M.dtype == 'float32'
        if is_float:
            F = Faust(core_obj=_FaustCorePy.FaustAlgoGenFlt.butterfly_hierarchical(*args))
        else:
            F = Faust(core_obj=_FaustCorePy.FaustAlgoGenDbl.butterfly_hierarchical(*args))
    else:
        F = Faust(core_obj=_FaustCorePy.FaustAlgoGenCplxDbl.butterfly_hierarchical(*args))
    if diag_opt:
        F = pyfaust.opt_butterfly_faust(F)
    return F
