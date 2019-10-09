%==========================================================================================
%> @brief Constructs a Faust implementing the Walsh-Hadamard Transform of order n.
%>
%> The resulting Faust has log2(n) sparse factors of order n, each one having 2 non-zero elements
%> per row and per column.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b H = wht(n) <br/>
%> &nbsp;&nbsp;&nbsp; @b H = wht(n, normed)
%>
%> @param n the power of two exponent for a Hadamard matrix of order n and a factorization into log2(n) factors.
%> @param normed: (optional) true (by default) to normalize the returned Faust as if Faust.normalize() was called, false otherwise.
%>
%> @retval H the Faust implementing the Hadamard transform of dimension n.
%>
%> @b Example
%> @code
%> % in a matlab terminal
%> >> import matfaust.*
%> >> H = wht(1024) % is equal to
%> >> H = normalize(wht(1024, false))
%> @endcode
%>
%>
%>H =
%>
%>Faust size 1024x1024, density 0.0195312, nnz_sum 20480, 10 factor(s):
%>- FACTOR 0 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>- FACTOR 1 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>- FACTOR 2 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>- FACTOR 3 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>- FACTOR 4 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>- FACTOR 5 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>- FACTOR 6 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>- FACTOR 7 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>- FACTOR 8 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>- FACTOR 9 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>
%==========================================================================================
function H = wht(n, varargin)
	% check n (must be integer > 0)
	if(~ isreal(n) || n < 0 || abs(n-floor(n)) > 0)
		error('n must be an integer greater than zero')
	end
	log2n = floor(log2(n));
	if(2^log2n < n)
		error('n must be a power of 2')
	end
	if(log2n > 31)
		error('Can''t handle a Hadamard Faust of order larger than 2^31')
	end
	if(length(varargin) > 0)
		if(~ islogical(varargin{1}))
			error('wht optional second argument must be a boolean');
		end
		normed = varargin{1};
	else
		normed = true; % normalization by default
	end
	core_obj = mexFaustReal('hadamard', log2n, normed);
	is_real = true;
	e = MException('FAUST:OOM', 'Out of Memory');
	if(core_obj == 0)
		throw(e)
	end
	H = matfaust.Faust(core_obj, is_real);
end
