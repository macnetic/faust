%==========================================================================================
%> @brief Returns an anticirculant Faust A defined by the vector c (which is the last column of full(A)).
%>
%> @b Example
%>
%> @code
%> % in a matlab terminal
%> >> import matfaust.anticirc
%> >> c = 1:8;
%> >> A = anticirc(c)
%> @endcode
%>
%> A =
%>
%> Faust size 8x8, density 1.5, nnz_sum 96, 6 factor(s):
%> - FACTOR 0 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 1 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 2 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 3 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 4 (complex) SPARSE, size 8x8, density 0.25, nnz 16
%> - FACTOR 5 (complex) SPARSE, size 8x8, density 0.25, nnz 16
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
%> >> real(full_A)
%>
%> ans =
%>
%>    2.0000    3.0000    4.0000    5.0000    6.0000    7.0000    8.0000    1.0000
%>    3.0000    4.0000    5.0000    6.0000    7.0000    8.0000    1.0000    2.0000
%>    4.0000    5.0000    6.0000    7.0000    8.0000    1.0000    2.0000    3.0000
%>    5.0000    6.0000    7.0000    8.0000    1.0000    2.0000    3.0000    4.0000
%>    6.0000    7.0000    8.0000    1.0000    2.0000    3.0000    4.0000    5.0000
%>    7.0000    8.0000    1.0000    2.0000    3.0000    4.0000    5.0000    6.0000
%>    8.0000    1.0000    2.0000    3.0000    4.0000    5.0000    6.0000    7.0000
%>    1.0000    2.0000    3.0000    4.0000    5.0000    6.0000    7.0000    8.0000
%> >> % Look at the density of a larger anticirculant Faust
%> >> % it indicates a speedup of the Faust-matrix/vector product
%> >> density(anticirc(rand(1, 1024)))
%> 0.0391

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
