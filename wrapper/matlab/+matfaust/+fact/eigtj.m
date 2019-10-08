%==========================================================================================
%> @brief Computes the eigenvalues and the eigenvectors transform (as a Faust object) using the truncated Jacobi algorithm.
%>
%> The eigenvalues and the eigenvectors are approximate. The trade-off between accuracy and sparsity can be set through the parameters J and t.
%>
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b eigtj(M, J) calls the non-parallel Givens algorithm.<br/>
%> &nbsp;&nbsp;&nbsp; @b eigtj(M, J, 0) or eigtj(M, J, 1) do the same as in previous line.<br/>
%> &nbsp;&nbsp;&nbsp; @b eigtj(M, J, t) calls the parallel Givens algorithm (if t > 1, otherwise it calls basic Givens algorithm)<br/>
%> &nbsp;&nbsp;&nbsp; @b eigtj(M, J, t, 'verbosity', 2) same as above with a level of verbosity of 2 in output. <br/>
%>
%> @param M the matrix to diagonalize. Must be real and symmetric.
%> @param J defines the number of factors in the eigenvector transform V. The number of factors is round(J/t). Note that the last permutation factor is not in count here (in fact, the total number of factors in V is rather round(J/t)+1).
%> @param t the number of Givens rotations per factor. Note that t is forced to the value min(J,t). Besides, a value of t such that t > size(M,1)/2 won't lead to the desired effect because the maximum number of rotation matrices per factor is anyway size(M,1)/2. The parameter t is meaningful in the parallel version of the truncated Jacobi algorithm (cf. references below). If t <= 1 (by default) then the function runs the non-parallel algorithm.
%> @param verbosity the level of verbosity, the greater the value the more info. is displayed.
%>
%> @retval [V,D]
%> - V the Faust object representing the approximate eigenvector transform. V has its last factor being a permutation matrix, the goal of this factor is to apply to the columns of V the same order as eigenvalues set in D.
%> - D the approximate sparse diagonal matrix of the eigenvalues (in ascendant order along the rows/columns).
%>
%> @b Example
%> @code
%> import matfaust.fact.eigtj
%>
%> % get a Laplacian to diagonalize
%> load('Laplacian_256_community.mat')
%> % do it
%> [Uhat, Dhat] = eigtj(Lap, size(Lap,1)*100, size(Lap, 1)/2, 'verbosity', 2)
%> % Uhat is the Fourier matrix/eigenvectors approximattion as a Faust (200 factors + permutation mat.)
%> % Dhat the eigenvalues diagonal matrix approx.
%> @endcode
%>
%>
%>
%> @b References:
%> - [1]   Le Magoarou L., Gribonval R. and Tremblay N., "Approximate fast
%> graph Fourier transforms via multi-layer sparse approximations",
%> IEEE Transactions on Signal and Information Processing
%> over Networks 2018, 4(2), pp 407-420
%>
%> <p> @b See @b also fact.fgft_givens, fact.fgft_palm
%>
%==========================================================================================
function [V,D] = eigtj(M, J, varargin)
	[V, D] = matfaust.fact.fgft_givens(M, J, varargin{:});
	V = factors(V,1:numfactors(V))
	% copy seems unecessary but it's to workaround a bug (temporarily)
end
