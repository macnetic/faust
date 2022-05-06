%==========================================================================================
%> @brief Returns a circulant Faust C defined by the vector c (which is the first column of the full(C)).
%>
%> @b See also matfaust.anticirc
%==========================================================================================
function C = circ(c)
	log2c = log2(numel(c));
	if(log2c ~= floor(log2c))
		error('Only power of two length vectors are supported')
	end
	n = numel(c);
	F = matfaust.dft(n, 'normed', false);
	FH = F';
	if (size(c, 1) < size(c, 2))
		c = c.';
	end
	S = sparse(diag(FH*(c/n)));
	C = F * matfaust.Faust(S*factors(FH, 1)) * right(FH, 2);
end
