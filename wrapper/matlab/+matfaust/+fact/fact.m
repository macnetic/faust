%> @package matfaust.factparams @brief The module for the parametrization of FAuST's algorithms (Palm4MSA and Hierarchical Factorization)
%> @package matfaust.fact @brief The matfaust factorization module.
%>
%>
%>    This module gives access to the main factorization algorithms of
%>    FAuST. These algorithms can factorize a dense matrix to a sparse product
%>    (i.e. a Faust object).
%>
%>    There are several factorization algorithms.
%>
%>    - The first one is PALM4MSA:
%>    which stands for Proximal Alternating Linearized Minimization for
%>    Multi-layer Sparse Approximation. Note that PALM4MSA is not
%>    intended to be used directly. You should rather rely on the second algorithm.
%>
%>    - The second one is the Hierarchical Factorization algorithm:
%>    this is the central algorithm to factorize a dense matrix to a Faust.
%>    It makes iterative use of Palm4MSA to proceed with the factorization of a given
%>    dense matrix.
%>
%>    - The third group of algorithms is for approximate eigenvalue decomposition (eigtj) and singular value decomposition (svdtj).
%>    - The fourth algorithm is fact.butterfly.
%>
%>
%>
%======================================================================

