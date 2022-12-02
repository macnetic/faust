%> @package matfaust.lazylinop @brief The matfaust module for lazy linear operators.
%=============================================================
%> @brief Returns the LazyLinearOp for eye.
%>
%> @param m: (int) the number of rows.
%> @param n: (int, optional) the number of columns.
%> @param k: (int, optional) diagonal to place ones on. Default is 0 (main diagonal). Negative integer for a diagonal below the main diagonal, strictly positive integer for a diagonal above.
%> @param 'dtype', str: data type of the LazyLinearOp ('double', 'single', 'complex').
%>
%> @b Example:
%> @code
%>
%> >> import matfaust.lazylinop.eye
%> >> le1 = eye(5)
%>
%> le1 =
%>
%>   5x5 LazyLinearOp array with no properties.
%>
%> >> full(le1)
%>
%> ans =
%>
%>      1     0     0     0     0
%>      0     1     0     0     0
%>      0     0     1     0     0
%>      0     0     0     1     0
%>      0     0     0     0     1
%>
%> >> le2 = eye(5, 2)
%>
%> le2 =
%>
%>   5x2 LazyLinearOp array with no properties.
%>
%> >> full(le2)
%>
%> ans =
%>
%>      1     0
%>      0     1
%>      0     0
%>      0     0
%>      0     0
%>
%> >> le3 = eye(5, 3, 1)
%>
%> le3 =
%>
%>   5x3 LazyLinearOp array with no properties.
%>
%> >> full(le3)
%>
%> ans =
%>
%>      0     1     0
%>      0     0     1
%>      0     0     0
%>      0     0     0
%>      0     0     0
%>
%> >> le4 = eye(5, 3, -1)
%>
%> le4 =
%>
%>   5x3 LazyLinearOp array with no properties.
%>
%> >> full(le4)
%>
%> ans =
%>
%>      0     0     0
%>      1     0     0
%>      0     1     0
%>      0     0     1
%>      0     0     0
%> @endcode
%>
%=============================================================
function EL = eye(m, varargin)
	import matfaust.lazylinop.*
	p = inputParser;


	validK = @(k) isscalar(k) && 0 == k - floor(k);
	addOptional(p, 'n', 'undefined', validK)

	addOptional(p, 'k', 0, validK)

	%addParameter(p, 'n', 'undefined', validK)

	%addParameter(p, 'k', 0, validK)

	validDtype = @(dtype) any(strcmp(dtype, {'complex', 'double', 'single', 'undefined'}));
	% TODO: use dtype
	addParameter(p, 'dtype', 'double', validDtype)



	parse(p, varargin{:})
	dtype = p.Results.dtype;
	k = p.Results.k;
	n = p.Results.n;

	if strcmp(n, 'undefined')
		n = m;
	end

	function P = matmat(x, m, n, k)
		if n ~= size(x, 1)
			error('Dimensions must agree')
		end
		minmn = min(m, n);
		x_islop = isLazyLinearOp(x);
		if k < 0
			neg_k = true;
			nz = zeros(abs(k), size(x, 2));
			if x_islop
				nz = aslazylinearoperator(nz);
			end
			limk = min(minmn, m - abs(k));
			k = 0;
		else
			limk = min(minmn, n -k);
			neg_k = false;
		end
		mul = x(k+1: k + limk, :);
		if neg_k
			mul = [ nz ; mul ];
		end
		if size(mul, 1) < m
			z = zeros(m - size(mul,1), size(mul, 2));
			if x_islop
				z = azlazylinearoperator(z);
			end
			mul = [mul ; z];
		end
		P = mul;
	end
	EL = LazyLinearOperator([m, n], 'matmat', @(x) matmat(x, m, n, k), 'rmatmat', @(x) matmat(x, n, m, -k));
end
