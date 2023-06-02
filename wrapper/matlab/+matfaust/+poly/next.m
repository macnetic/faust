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
%> >> rng(42)
%> >> L = sprand(50, 50, .2);
%> >> L = L*L';
%> >> K = 3;
%> >> F = basis(L, K, 'chebyshev');
%> >> G = next(F);
%> >> F
%>
%> F =
%>
%> Faust size 200x50, density 0.6624, nnz_sum 6624, 4 factor(s):
%> - FACTOR 0 (double) SPARSE, size 200x150, density 0.0752667, nnz 2258
%> - FACTOR 1 (double) SPARSE, size 150x100, density 0.1472, nnz 2208
%> - FACTOR 2 (double) SPARSE, size 100x50, density 0.4216, nnz 2108
%> - FACTOR 3 (double) SPARSE, size 50x50, density 0.02, nnz 50
%>  identity matrix flag
%>
%> >> G
%>
%> G =
%>
%> Faust size 250x50, density 0.71456, nnz_sum 8932, 5 factor(s):
%> - FACTOR 0 (double) SPARSE, size 250x200, density 0.04616, nnz 2308
%> - FACTOR 1 (double) SPARSE, size 200x150, density 0.0752667, nnz 2258
%> - FACTOR 2 (double) SPARSE, size 150x100, density 0.1472, nnz 2208
%> - FACTOR 3 (double) SPARSE, size 100x50, density 0.4216, nnz 2108
%> - FACTOR 4 (double) SPARSE, size 50x50, density 0.02, nnz 50
%>  identity matrix flag
%>
%> @endcode
%>
%======================================================================
function G = next(F)
	G = next(F);
end
