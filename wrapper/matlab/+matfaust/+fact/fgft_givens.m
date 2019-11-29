% experimental block start
%==========================================================================================
%> @brief Computes the FGFT of the Laplacian matrix Lap (using fact.eigtj).
%>
%>
%> @b Usage
%>  &nbsp;&nbsp;&nbsp; @b See fact.eigtj
%>
%> @param Lap the Laplacian matrix (which is symmetric or hermitian).
%> @param maxiter see fact.eigtj
%> @param 'nGivens_per_fac', integer see fact.eigtj
%> @param 'tol', number see fact.eigtj
%> @param 'relerr', bool see fact.eigtj
%> @param 'verbosity', integer see fact.eigtj
%> @param 'order', integer see fact.eigtj
%>
%> @retval [FGFT,D]:
%> - with FGFT being the Faust object representing the Fourier transform and,
%> -  D as a sparse diagonal matrix of the eigenvalues by default in ascendant order along the rows/columns.
%>
%>
%> @b References:
%> - [1]   Le Magoarou L., Gribonval R. and Tremblay N., "Approximate fast
%> graph Fourier transforms via multi-layer sparse approximations",
%> IEEE Transactions on Signal and Information Processing
%> over Networks 2018, 4(2), pp 407-420 <https://hal.inria.fr/hal-01416110>
%>
%>
%> <p> @b See @b also fact.eigtj, fact.fgft_palm
%>
%==========================================================================================
function [FGFT,D] = fgft_givens(Lap, maxiter, varargin)
	[FGFT, D] = matfaust.fact.eigtj(Lap, maxiter, varargin{:});
end
% experimental block end
