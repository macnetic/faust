function EL = eye(m, varargin)
	import matfaust.lazylinop.*
	p = inputParser;


	validK = @(k) isscalar(k) && 0 == k - floor(k);
	addParameter(p, 'n', 'undefined', validK)

	addParameter(p, 'k', 0, validK)

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
