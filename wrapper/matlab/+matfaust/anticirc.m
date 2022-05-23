%==========================================================================================
%> @brief Returns an anticirculant Faust A defined by the vector c (which is the last column of full(A)).
%>
%> @b Example
%>
%> @code
%> % in a matlab terminal
%> >> import matfaust.anticirc
%> >> c = rand(1, 8);
%> >> A = anticirc(c)
%> @endcode
%>
%> A =
%>
%> Faust size 8x8, density 1.75, nnz_sum 112, 8 factor(s):
%> - FACTOR 0 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 1 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 2 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 3 (complex) SPARSE, size 8x8, density 0.125, nnz 8
%> - FACTOR 4 (complex) SPARSE, size 8x8, density 0.125, nnz 8
%> - FACTOR 5 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 6 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 7 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%>
%> @code
%> >> full_A = full(A);
%> >> all(full_A(:, end).' - c < 1e-15)
%>
%> ans =
%>
%>   logical
%>
%>    1
%>
%> >> c
%>
%> c =
%>
%>     0.0046    0.7749    0.8173    0.8687    0.0844    0.3998    0.2599    0.8001
%>
%> >> real(full_A)
%>
%> ans =
%>
%>     0.7749    0.8173    0.8687    0.0844    0.3998    0.2599    0.8001    0.0046
%>     0.8173    0.8687    0.0844    0.3998    0.2599    0.8001    0.0046    0.7749
%>     0.8687    0.0844    0.3998    0.2599    0.8001    0.0046    0.7749    0.8173
%>     0.0844    0.3998    0.2599    0.8001    0.0046    0.7749    0.8173    0.8687
%>     0.3998    0.2599    0.8001    0.0046    0.7749    0.8173    0.8687    0.0844
%>     0.2599    0.8001    0.0046    0.7749    0.8173    0.8687    0.0844    0.3998
%>     0.8001    0.0046    0.7749    0.8173    0.8687    0.0844    0.3998    0.2599
%>     0.0046    0.7749    0.8173    0.8687    0.0844    0.3998    0.2599    0.8001
%> @endcode
%>
%> @b See also matfaust.circ, matfaust.toeplitz
%==========================================================================================
function A = anticirc(c)
	C = matfaust.circ(c);
	n = numel(c);
	I = n:-1:1;
	J = 1:n;
	P = sparse(I, J, 1);
	N = numfactors(C);
	A = left(C, N-1) * matfaust.Faust(factors(C, N) * P);
end
