%==========================================================================================
%> @brief Runs the truncated Jacobi algorithm to compute the eigenvalues of M (returned in D) and the corresponding transform of eigenvectors (in Faust V columns).
%>
%> The output is such that V*D*V' approximates M.
%>
%> The trade-off between accuracy and sparsity can be set through the parameters maxiter and nGivens_per_fac.
%>
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b eigtj(M, J) computes only one Givens rotation per factor of V.<br/>
%> &nbsp;&nbsp;&nbsp; @b eigtj(M, J, 'nGivens_per_fac', 1) do the same as in previous line.<br/>
%> &nbsp;&nbsp;&nbsp; @b eigtj(M, J, 'nGivens_per_fac', t) as above but with t Givens rotations per V factor.<br/>
%> &nbsp;&nbsp;&nbsp; @b eigtj(M, J, 'nGivens_per_fac', t, 'verbosity', 2) same as above with a level of verbosity of 2 in output. <br/>
%> &nbsp;&nbsp;&nbsp; @b eigtj(M, J, 'nGivens_per_fac', t, 'tol', 0.01, 'relerr', true) Uses a stopping criterion based on relative squared error norm(V*D*V'-M, 'fro')^2/norm(M, 'fro')^2. This criterion is concurrent to maxiter (here J).<br/>
%> &nbsp;&nbsp;&nbsp; @b eigtj(M, J, 'nGivens_per_fac', t, 'tol', 0.01, 'relerr', true) Uses a stopping criterion based on absolute squared error norm(V*D*V'-M, 'fro')^2. This criterion is concurrent to maxiter (here J).<br/>
%>
%> @param M the matrix to diagonalize. Must be real and symmetric.
%> @param maxiter defines the number of Givens rotations that are computed in eigenvector transform V. The number of rotations per factor of V is defined by nGivens_per_fac.
%> @param 'nGivens_per_fac', integer the number of Givens rotations per factor of V, must be an integer between 1 to floor(size(M, 1)/2) which is the default value.
%> @param 'tol', number (optional) the tolerance error under what the algorithm stops. By default, it's zero for not stopping on error criterion.
%> @param 'relerr', true (optional) For a stopping criterion based on the relative squared error (this is the default error).
%> @param 'relerr', false (optional) For a stopping criterion based on the absolute squared error.
%> @param 'verbosity', integer (optional) the level of verbosity, the greater the value the more info. is displayed.
%> @param 'order', char (optional) 'descend' for a descending order of eigenvalues, 'ascend' for an ascending order (default value) or 'undef' for no sort.
%>
%> @retval [V,D]
%> - V the Faust object representing the approximate eigenvector transform. The column V(:, i) is the eigenvector corresponding to the eigenvalue D(i,i).
%> The last factor of V is a permutation matrix. The goal of this factor is to apply to the columns of V the same order as eigenvalues set in D.
%> - D the approximate sparse diagonal matrix of the eigenvalues (by default in ascendant order along the rows/columns).
%>
%> @b Example
%> @code
%> import matfaust.fact.eigtj
%>
%> % get a Laplacian to diagonalize
%> load('Laplacian_256_community.mat')
%> % do it
%> [Uhat, Dhat] = eigtj(Lap, size(Lap,1)*100, 'nGivens_per_fac', size(Lap, 1)/2, 'verbosity', 2)
%> % Uhat is the Fourier matrix/eigenvectors approximation as a Faust (200 factors + permutation mat.)
%> % Dhat the eigenvalues diagonal matrix approx.
%> @endcode
%>
%>
%>
%> @b References:
%> - [1]   Le Magoarou L., Gribonval R. and Tremblay N., "Approximate fast
%> graph Fourier transforms via multi-layer sparse approximations",
%> IEEE Transactions on Signal and Information Processing
%> over Networks 2018, 4(2), pp 407-420 <https://hal.inria.fr/hal-01416110>
%>
%> <p> @b See @b also fact.fgft_givens, fact.fgft_palm
%>
%==========================================================================================
function [V,D] = eigtj(M, maxiter, varargin)
	[V, D] = matfaust.fact.fgft_givens(M, maxiter, varargin{:});
	%V = factors(V,1:numfactors(V)) % tmp workaround
end
