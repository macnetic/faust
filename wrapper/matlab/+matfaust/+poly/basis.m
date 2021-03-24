% experimental block start
%======================================================================
%> @brief  Builds the Faust of the polynomial basis defined on the sparse matrix L.
%>
%> @param L the sparse square matrix.
%> @param K	the degree of the last polynomial, i.e. the K+1 first polynomials are built.
%> @param basis_name 'chebyshev', and others yet to come.
%>
%> @retval F the Faust of the basis composed of the K+1 orthogonal polynomials.
%>
%> @b Example
%> @code
%> % in a matlab terminal
%> >> import matfaust.poly.*
%> >> L = sprand(50, 50, .2);
%> >> L = L*L';
%> >> K = 3;
%> >> F = basis(L, K, 'chebyshev')
%> @endcode
%> F =
%>
%> Faust size 200x50, density 0.6612, nnz_sum 6612, 4 factor(s):
%> - FACTOR 0 (real) SPARSE, size 200x150, density 0.0751333, nnz 2254
%> - FACTOR 1 (real) SPARSE, size 150x100, density 0.146933, nnz 2204
%> - FACTOR 2 (real) SPARSE, size 100x50, density 0.4208, nnz 2104
%> - FACTOR 3 (real) SPARSE, size 50x50, density 0.02, nnz 50
%>
%=======================================================================
function F = basis(L, K, basis_name, varargin)
	% TODO: T0 in varargin
	if(~ ismatrix(L) || ~ issparse(L) || size(L, 1) ~= size(L, 2))
		error("L must be a square matrix.")
	end

	is_real = isreal(L);

	if(~ isstr(basis_name) && ~ ischar(basis_name))
		error('basis_name must be a character string/array.')
	end


	if(strcmp(basis_name, 'chebyshev'))
		if(is_real)
			core_obj = mexPolyReal('chebyshev', L, K);
		else
			core_obj = mexPolyCplx('chebyshev', L, K);
		end
	else
		error(['unknown basis name: ' basis_name])
	end

	F = matfaust.Faust(core_obj, is_real);
end
% experimental block end
