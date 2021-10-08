%> @package matfaust.poly @brief The matfaust module for polynomial basis as Faust objects.
%> @note This module is still in BETA status.

%======================================================================
%> @brief Computes the linear combination of the polynomials defined by basis.
%>
%> @param coeffs the linear combination coefficients (vector).
%> @param basis either the name of the polynomial basis to build on L or the basis if already built externally (as a Faust or an equivalent full array).
%> @param 'L', matrix the sparse matrix on which the polynomial basis is built if basis is not already a Faust or a full array.
%> @param 'X', matrix if X is set, the linear combination of basis*X is computed (note that the memory space is optimized compared to the manual way of doing first B = basis*X and then calling poly on B without X set).
%> @param 'dev', str (optional): the computating device ('cpu' or 'gpu').
%> @param 'dtype', str (optional): to decide in which data type the resulting Faust or array will be encoded ('float' or 'double' by default). If basis is a Faust or an array its dtype/class is prioritary over this parameter which is in fact useful only if basis is the name of the basis (a str/char array).
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
	X = {}; % by default no X argument is passed, set it as a cell (see why in matfaust.Faust.poly)
	argc = length(varargin);
	dev = 'cpu';
	dtype = 'double';
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
				case 'dtype'
					if(argc == i || ~ strcmp(tmparg, 'float') && ~ startsWith(tmparg, 'double'))
						error('dtype keyword argument is not followed by a valid value: float or double.')
					else
						dtype = tmparg;
					end
				otherwise
					if((isstr(varargin{i}) || ischar(varargin{i}))  && ~ strcmp(tmparg, 'cpu') && ~ startsWith(tmparg, 'gpu'))
						error([ tmparg ' unrecognized argument'])
					end
			end
		end
	end

	on_gpu = startsWith(dev , 'gpu');

	if(~ ismatrix(coeffs)) %|| size(coeffs, 1) ~= 1 && size(coeffs, 2) ~= 1)
		error('coeffs must be a scalar matrix')
	end

	if(isreal(coeffs) ~= isreal(basis))
		error('coeffs and basis must be of the same scalar type (real or complex)')
	end

	if isvector(coeffs)
		K = numel(coeffs)-1;
	else
		K = size(coeffs, 1)-1;
	end

	if(isstr(basis) || ischar(basis))
		if(exist('L') ~= 1)
			error('L key-value pair argument is missing in the argument list.')
		end
		basis = matfaust.poly.basis(L, K, basis, 'dev', dev, 'dtype', dtype);
	end


	is_real = isreal(basis);
	is_float = strcmp(class(basis), 'single') || strcmp('dtype', 'float');
	if(matfaust.isFaust(basis))
		if(numfactors(basis) ~= numel(coeffs))
			error('coeffs and basis dimensions must agree.')
		end
		LC = matfaust.poly.FaustPoly.poly_(basis, coeffs, X);
		% LC is a Faust iff X is not set, otherwise it's a matrix
	elseif(ismatrix(basis))
		if(mod(size(basis, 1), K+1) ~= 0)
			error('coeffs and basis dimensions must agree.')
		end
		d = floor(size(basis,1) / (K+1));
		if(size(coeffs, 2) == 1)
			if(is_real)
				if(is_float)
					LC = mexPolyRealFloat('polyMatrix', d, K, size(basis,2), coeffs, basis, on_gpu);
				else
					LC = mexPolyReal('polyMatrix', d, K, size(basis,2), coeffs, basis, on_gpu);
				end
			else
				LC = mexPolyCplx('polyMatrix', d, K, size(basis,2), coeffs, basis, on_gpu);
			end
		else
			if(is_real)
				if(is_float)
					LC = mexPolyRealFloat('polyMatrixGroupCoeffs', d, K, size(basis,2), size(coeffs, 2), coeffs, basis, on_gpu);
				else
					LC = mexPolyReal('polyMatrixGroupCoeffs', d, K, size(basis,2), size(coeffs, 2), coeffs, basis, on_gpu);
				end
			else
				LC = mexPolyCplx('polyMatrixGroupCoeffs', d, K, size(basis,2), size(coeffs, 2), coeffs, basis, on_gpu);
			end

		end
		% LC is a matrix
	end
end
