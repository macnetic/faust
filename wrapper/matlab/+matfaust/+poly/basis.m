%======================================================================
%> @brief  Builds the Faust of the polynomial basis defined on the sparse matrix L.
%>
%> @param L the sparse square matrix.
%> @param K	the degree of the last polynomial, i.e. the K+1 first polynomials are built.
%> @param basis_name 'chebyshev', and others yet to come.
%> @param 'T0', matrix (optional): a sparse matrix to replace the identity as a 0-degree polynomial of the basis.
%> @param 'dev', str (optional): the computing device ('cpu' or 'gpu').
%> @param 'dtype', str (optional): to decide in which data type the resulting Faust will be encoded ('float' or 'double' by default).
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
%>
%> F =
%>
%> Faust size 200x50, density 0.6612, nnz_sum 6612, 4 factor(s):
%> - FACTOR 0 (double) SPARSE, size 200x150, density 0.0751333, nnz 2254
%> - FACTOR 1 (double) SPARSE, size 150x100, density 0.146933, nnz 2204
%> - FACTOR 2 (double) SPARSE, size 100x50, density 0.4208, nnz 2104
%> - FACTOR 3 (double) SPARSE, size 50x50, density 0.02, nnz 50
%>  identity matrix flag
%> @endcode
%>
%> By default, the 0-degree polynomial is the identity.
%> However it is possible to replace the corresponding matrix by any sparse matrix T0 of your choice (with the only constraint that size(T0,1) == size(L, 1)). In that purpose, do as follows:
%>
%> @code
%> >> rng(42)
%> >> T0 = sprand(50, 2, .3);
%> >> F = basis(L, 3, 'chebyshev', 'T0', T0)
%> F =
%>
%> Faust size 200x2, density 16.53, nnz_sum 6612, 4 factor(s):
%> - FACTOR 0 (double) SPARSE, size 200x150, density 0.0751333, nnz 2254
%> - FACTOR 1 (double) SPARSE, size 150x100, density 0.146933, nnz 2204
%> - FACTOR 2 (double) SPARSE, size 100x50, density 0.4208, nnz 2104
%> - FACTOR 3 (double) SPARSE, size 50x2, density 0.5, nnz 50
%> @endcode
%>
%>
%=======================================================================
function F = basis(L, K, basis_name, varargin)

	if(~ ismatrix(L) || ~ issparse(L) || size(L, 1) ~= size(L, 2))
		error("L must be a square matrix.")
	end

	is_real = isreal(L);

	if(~ isstr(basis_name) && ~ ischar(basis_name))
		error('basis_name must be a character string/array.')
	end

	T0_is_set = false;
	T0 = []; % no T0 by default
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
				case 'T0'
					if(argc == i || ~ ismatrix(tmparg))
						error('T0 argument must be followed by a matrix.')
					else
						T0 = tmparg;
						T0_is_set = true;
					end
				case 'dev' % not used yet
					if(argc == i || ~ strcmp(tmparg, 'cpu') && ~ startsWith(tmparg, 'gpu'))
						error('dev keyword argument is not followed by a valid value: cpu, gpu*.')
					else
						dev = tmparg;
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

	is_float = strcmp(class(L), 'single') || strcmp(dtype, 'float'); % since it's impossible to create a float sparse matrix in matlab use a dtype argument
	mex_args = {basis_name, L, K, startsWith(dev, 'gpu')};

	if(T0_is_set)
		if(size(T0,1) ~= size(L,2))
			error('T0 number of rows must be equal to L number of columns/rows')
		end
		if(~ issparse(T0))
			error('T0 must be a sparse matrix')
		end
		mex_args = [mex_args {T0}];
	end

	if(strcmp(basis_name, 'chebyshev'))
		if(is_real)
			if(is_float)
				core_obj = mexPolyRealFloat(mex_args{:});
			else
				core_obj = mexPolyReal(mex_args{:});
			end
		else
			core_obj = mexPolyCplx(mex_args{:});
		end
	else
		error(['unknown basis name: ' basis_name])
	end

	F = matfaust.poly.FaustPoly(core_obj, is_real, 'cpu', dtype);
end
