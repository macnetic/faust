% experimental block start
%> @package matfaust.poly @brief The matfaust module for polynomial basis as Faust objects.

%======================================================================
%> @brief Computes the linear combination of the polynomials defined by basis.
%>
%> @param coeffs the linear combination coefficients (vector).
%> @param basis either the name of the polynomial basis to build on L or the basis if already built externally (as a Faust or an equivalent full array).
%> @param 'L', matrix the sparse matrix on which the polynomial basis is built if basis is not already a Faust or a full array.
%> @param 'X', matrix if X is set, the linear combination of basis*X is computed (note that the memory space is optimized compared to the manual way of doing first B = basis*X and then calling poly on B witout X set).
%>
%> @retval LC The linear combination Faust or full array depending on if basis is itself a Faust or a np.ndarray.
%>
%> @b Example
%> @code
%> >> import matfaust.poly.*
%> >> L = sprand(50, 50, .02);
%> >> L = L*L';
%> >> coeffs = [.5 1 2 3];
%> >> G = poly(coeffs, 'chebyshev', 'L', L)
%> @endcode
%>
%> G =
%>
%> Faust size 50x50, density 0.3524, nnz_sum 881, 5 factor(s):
%> - FACTOR 0 (real) SPARSE, size 50x200, density 0.02, nnz 200
%> - FACTOR 1 (real) SPARSE, size 200x150, density 0.00923333, nnz 277
%> - FACTOR 2 (real) SPARSE, size 150x100, density 0.0151333, nnz 227
%> - FACTOR 3 (real) SPARSE, size 100x50, density 0.0254, nnz 127
%> - FACTOR 4 (real) SPARSE, size 50x50, density 0.02, nnz 50
%>
%> Which is equivalent to do as below (in two times):
%>
%> @code
%> >> K = 3;
%> >> F = basis(L, K, 'chebyshev');
%> >> G = poly(coeffs, F)
%> @endcode
%>
%> G =
%>
%> Faust size 50x50, density 0.3524, nnz_sum 881, 5 factor(s):
%> - FACTOR 0 (real) SPARSE, size 50x200, density 0.02, nnz 200
%> - FACTOR 1 (real) SPARSE, size 200x150, density 0.00923333, nnz 277
%> - FACTOR 2 (real) SPARSE, size 150x100, density 0.0151333, nnz 227
%> - FACTOR 3 (real) SPARSE, size 100x50, density 0.0254, nnz 127
%> - FACTOR 4 (real) SPARSE, size 50x50, density 0.02, nnz 50
%>
%> Above G is a Faust because F is too. Below the full array of the Faust F is passed, so an array is returned into GA.
%>
%> @code
%> GA = poly(coeffs, full(F));
%> ismatrix(GA)
%> @endcode
%>
%>ans =
%>
%>	logical
%>
%>	1
%>
%> But of course they are equal:
%>
%> @code
%> >> norm(GA-full(G))/(norm(GA))
%> @endcode
%>
%>ans =
%>
%>   8.8455e-17
%>
%======================================================================
function LC = poly(coeffs, basis, varargin)
	dev = 'cpu';
	argc = length(varargin);
	X = {}; % by default no X argument is passed, set it as a cell (see why in matfaust.Faust.poly)
	if(argc > 0)
		for i=1:2:argc
			if(argc > i)
				% next arg (value corresponding to the key varargin{i})
				tmparg = varargin{i+1};
			end
			switch(varargin{i})
				case 'L'
					if(argc == i || ~ ismatrix(tmparg) || ~ issparse(tmparg))
						error('L argument must be followed by a square and sparse matrix.')
					else
						L = tmparg;
					end
				case 'dev' % not used yet
					if(argc == i || ~ strcmp(tmparg, 'cpu') && ~ startsWith(tmparg, 'gpu'))
						error('dev keyword argument is not followed by a valid value: cpu, gpu*.')
					else
						dev = tmparg;
					end
				case 'X'
					if(argc == i || ~ ismatrix(tmparg) || issparse(tmparg))
						error('X argument must be followed by a dense matrix.')
					else
						X = tmparg;
					end
				otherwise
					if((isstr(varargin{i}) || ischar(varargin{i}))  && ~ strcmp(tmparg, 'cpu') && ~ startsWith(tmparg, 'gpu'))
						error([ tmparg ' unrecognized argument'])
					end
			end
		end
	end

	if(~ ismatrix(coeffs) || size(coeffs, 1) ~= 1 && size(coeffs, 2) ~= 1)
		error('coeffs must be a scalar vector')
		if(isreal(coeffs) ~= isreal(basis))
			error('coeffs and basis must be of the same scalar type (real or complex)')
		end
	end

	K = numel(coeffs)-1;

	if(isstr(basis) || ischar(basis))
		if(exist('L') ~= 1)
			error('L key-value pair argument is missing in the argument list.')
		end
		basis = matfaust.poly.basis(L, K, basis, 'dev', dev);
	end


	is_real = isreal(basis);
	if(matfaust.isFaust(basis))
		if(numfactors(basis) ~= numel(coeffs))
			error('coeffs and basis dimensions must agree.')
		end
		LC = matfaust.Faust.poly(basis, coeffs, 0, X);
		% LC is a Faust iff X is not set, otherwise it's a matrix
	elseif(ismatrix(basis))
		if(mod(size(basis, 1), K+1) ~= 0)
			error('coeffs and basis dimensions must agree.')
		end
		d = floor(size(basis,1) / (K+1));
		if(is_real)
			LC = mexPolyReal('polyMatrix', d, K, size(basis,2), coeffs, basis);
		else
			LC = mexPolyCplx('polyMatrix', d, K, size(basis,2), coeffs, basis);
		end
		% LC is a matrix
	end
end
%> experimental block end
