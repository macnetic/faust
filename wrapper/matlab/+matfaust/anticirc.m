%==========================================================================================
%> @brief Returns a anticirculant Faust C defined by the vector c (which is the last column of the full(C)).
%>
%> @b See also matfaust.circ, matfaust.toeplitz
%==========================================================================================
function C = anticirc(c)
	C = matfaust.circ(c)
	P = zeros(numel(c))
	n = numel(c)
	I = n:-1:1
	J = 1:n
	for k=1:n
		P(I(k), J(k)) = 1
	end
	N = numfactors(C)
	C = left(C, N-1) * matfaust.Faust(factors(C, N) * P)
end
