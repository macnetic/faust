% experimental block start
%======================================================================
%> @brief Gives the next Faust basis of dimension (n+1) from the Faust F polynomial basis of dimension n.
%>
%> @param F The polynomial basis Faust (must have been generated with matfaust.poly.basis).
%> @reval G the basis of dimension (n+1) with one additional factor added to those of F.
%>
%> @note The identical factors between F and G are just views of each other (they are not duplicated in memory).
%>
%> @b Example
%> @code
%> % in a matlab terminal
%> >> import matfaust.poly.*
%> >> L = sprand(50, 50, .2);
%> >> L = L*L';
%> >> K = 3;
%> >> F = basis(L, K, 'chebyshev');
%> >> G = next(F);
%> >> F
%> @endcode
%> F =
%>
%>Faust size 200x50, density 0.6672, nnz_sum 6672, 4 factor(s):
%>- FACTOR 0 (real) SPARSE, size 200x150, density 0.0758, nnz 2274
%>- FACTOR 1 (real) SPARSE, size 150x100, density 0.148267, nnz 2224
%>- FACTOR 2 (real) SPARSE, size 100x50, density 0.4248, nnz 2124
%>- FACTOR 3 (real) SPARSE, size 50x50, density 0.02, nnz 50
%>
%> @code
%> >> G
%> @endcode
%> G =
%>
%> Faust size 250x50, density 0.71968, nnz_sum 8996, 5 factor(s):
%> - FACTOR 0 (real) SPARSE, size 250x200, density 0.04648, nnz 2324
%> - FACTOR 1 (real) SPARSE, size 200x150, density 0.0758, nnz 2274
%> - FACTOR 2 (real) SPARSE, size 150x100, density 0.148267, nnz 2224
%> - FACTOR 3 (real) SPARSE, size 100x50, density 0.4248, nnz 2124
%> - FACTOR 4 (real) SPARSE, size 50x50, density 0.02, nnz 50
%>  identity matrix flag
%> @endcode
%>
%======================================================================
function G = next(F)
	G = next(F);
end
% experimental block end
