%==========================================================================================
%> @brief Generates a random Faust.
%>
%> @warning if this function is imported through 'import matfaust.rand' or 'import matfaust.*' the Matlab builtin function rand will be unreachable (matfaust.rand will be called instead). So, generally, it is not advisable to import this function, rather directly call matfaust.rand without any import.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b rand(@b M,@b N) generates a M-by-N Faust object. The numbers of rows and columns of intermediary factors are all randomly chosen between M and N (included). The number of factors is 5. The factors are sparse and real. The nnz per row of factors is 5.
%>
%> &nbsp;&nbsp;&nbsp; @b rand(@b M, @b N, @b 'num_factors', @b NF, @b 'dim_sizes', @b S) with NF and S two integers, generates a M-by-N Faust of NF factors. The factor size is S for both dimensions (except the first factor number of rows and the last factor number of columns which are respectively M and N). The factors are sparse and real.
%>
%> &nbsp;&nbsp;&nbsp; @b rand(@b M, @b N, @b 'num_factors', [@b N1, @b N2], @b 'dim_sizes', @b S) same as above except that here the number of factors is randomly chosen between N1 and N2 inclusively.
%>
%> &nbsp;&nbsp;&nbsp; @b rand(@b M, @b N, @b 'num_factors', [@b N1, N@b 2], [@b S1, @b S2]) or @b rand(@b M, N@b , @b 'num_factors', @b NF, @b 'dim_sizes', [@b S1, @b S2]) same as above except that here the intermediary factor dimension sizes are random; the number of rows and columns are both randomly chosen between S1 and S2 inclusively.
%>
%> &nbsp;&nbsp;&nbsp; @b rand(@b M, @b N, @b 'num_factors', @b NF, @b 'dim_sizes', @b S, @b 'density', @b D) or @b rand(@b M, @b N, @b 'num_factors', [@b N1, @b N2], @b 'dim_sizes', [@b S1, @b S2], @b 'density', @b D) same as above but specifying D the approximate density of each factor.
%>
%> &nbsp;&nbsp;&nbsp; @b rand(@b M, @b N, @b 'num_factors', @b NF, @b 'dim_sizes', @b S, @b 'density', @b D, @b 'per_row', @b true) or @b rand(@b M, @b N, @b 'num_factors', [@b N1, @b N2], @b 'dim_sizes', [@b S1, @b S2], @b 'density',@b D, @b 'per_row', true) same as above but specifying D, the density of each factor per row ('per_row', true) or per column ('per_row', false).
%>
%> &nbsp;&nbsp;&nbsp; @b @b rand(@b M, @b N, @b 'num_factors', @b NF, @b 'dim_sizes', @b S, @b 'density', @b D, @b 'fac_type', @b 'dense') or @b rand(@b M, @b N, @b 'num_factors', @b [@b  N1, @b N2], @b 'dim_sizes', [@b S1, @b S2], @b 'density', @b D, @b 'fac_type', @b 'dense') same as above but generating only dense matrices as factors.
%>
%> &nbsp;&nbsp;&nbsp; @b rand(@b M, @b N, @b 'num_factors', @b NF, @b 'dim_sizes', @b S, @b 'density', @b D, @b 'mixed') or @b rand(@b M, @b N, @b 'num_factors', [@b N1, @b N2], @b 'dim_sizes', [@b S1, @b S2], @b 'density', @b D, @b 'fac_type', @b 'sparse') same as above but generating either sparse or dense matrices as factors.
%>
%> &nbsp;&nbsp;&nbsp; @b rand(@b M, @b N, @b 'num_factors', @b NF, @b 'dim_sizes', @b S, @b 'density', @b D, @b 'fac_type', @b 'sparse', @b 'field' , @b 'complex'), @b rand(@b M, @b N, @b 'num_factors', [@b N1, @b N2], @b 'dim_sizes', [@b S1, @b S2], @b 'density', @b D, @b 'fac_type', @b 'sparse', @b 'per_row', @b false), rand(@b M, @b N, @b 'num_factors', @b NF, @b 'dim_sizes', @b S, @b 'density', @b D, @b 'fac_type', @b 'dense', @b 'field' , @b 'complex') or @b rand(@b M, @b N, @b 'num_factors', [@b N1, @b N2], @b 'dim_sizes', [@b S1, @b S2], @b D, @b 'fac_type', @b 'dense', @b 'field', @b 'complex') same as above but generating a complex Faust, that is, matrices defined over the complex field.
%>
%>
%>
%>
%>
%>
%>
%> @param M the Faust number of rows.
%> @param N the Faust number of columns.
%> @param 'num_factors', NF if it's an integer it is the number of random factors to set in the Faust.
%>                    If NF is a vector of 2 integers then the
%>                    number of factors is set randomly between
%>                    NF(1) and NF(2) (inclusively). Defaultly, a 5 factors long Faust is generated.
%> @param 'dim_sizes',S if it's an integer all Faust factors are square matrices (except maybe the first
%>			and last ones, depending on M and N). The size of the intermediary square factors is S**2.
%>			If it's a vector of 2 integers then the number of rows and columns are both a random number between size_dims(1) and
%> 			size_dims(2) (inclusively). Defaultly, dim_sizes is [M, N].
%> @param 'density',D the approximate density of generated factors.
%>		      D must be a floating point number greater than 0 and lower of equal to 1.
%> 	              The default value is such that each factor gets 5 non-zero
%> 	              elements per row, if per_row is true, or per column otherwise.
%> 	              A density of zero is equivalent to the default case.
%> @param 'per_row',bool this argument is to specify the density per row or per column.
%>		         By default the density is set per row and is such that the Faust's factors will have 5 non-zero elements per row.
%> @param fac_type, str the storage representation of factors. str must be 'sparse', 'dense' or 'mixed'.
%> 			The latter designates a mix of dense and sparse matrices in the generated Faust (the choice is made according
%>                      to a uniform distribution).
%>			The default value is 'sparse'.
%> @param 'field', str	str is either 'real' or 'complex' (the Faust field).
%>                      The default value is 'real'.
%> @param 'dev', str 'gpu or 'cpu' to create the random Faust on CPU or GPU (by default on CPU).
%> @param 'class', str 'double' (by default) or 'single' to select the scalar type used for the Faust generated.
%> @param 'seed', int  seed to initialize the PRNG used to generate the Faust factors. If 0 (default) the seed will be initialized with a random value depending of the time clock.
%>
%>
%> @retval F the random Faust.
%>
%> @b Example @b 1
%> @code
%> >> F = matfaust.rand(10,5)
%>
%> F =
%>
%> Faust size 10x5, density 4.32, nnz_sum 216, 5 factor(s):
%> - FACTOR 0 (real) SPARSE, size 10x9, density 0.555556, nnz 50
%> - FACTOR 1 (real) SPARSE, size 9x6, density 0.666667, nnz 36
%> - FACTOR 2 (real) SPARSE, size 6x10, density 0.5, nnz 30
%> - FACTOR 3 (real) SPARSE, size 10x10, density 0.5, nnz 50
%> - FACTOR 4 (real) SPARSE, size 10x5, density 1, nnz 50
%> @endcode
%>
%> @b Example @b 2
%> @code
%> % in a matlab terminal
%> >> G = matfaust.rand(10, 10, 'num_factors', 2, 'dim_sizes', 10, 'density', .5, 'mixed', 'complex')
%>
%> G =
%>
%> Faust size 10x10, density 1, nnz_sum 100, 2 factor(s):
%> - FACTOR 0 (complex) SPARSE, size 10x10, density 0.5, nnz 50
%> - FACTOR 1 (complex) DENSE, size 10x10, density 0.5, nnz 50
%> @endcode
%>
%> @b Example @b 3
%> @code
%> >> H = matfaust.rand(10, 10, 'num_factors', [2, 5], 'dim_sizes', [10, 20], 'density', .5, 'fac_type', 'dense')
%>
%> H =
%>
%> Faust size 5x10, density 4.64, nnz_sum 232, 3 factor(s):
%> - FACTOR 0 (real) DENSE,  size 5x14, density 0.5, nnz 35
%> - FACTOR 1 (real) DENSE,  size 14x17, density 0.470588, nnz 112
%> - FACTOR 2 (real) DENSE,  size 17x10, density 0.5, nnz 85
%>
%> @endcode
%>
%> <p>@b See @b also Faust.Faust, matfaust.rand_bsr.
%==========================================================================================
function F = rand(M, N, varargin)
	import matfaust.*
	%> Identifies a complex Faust.
	COMPLEX=3;
	%> Identifies a real Faust.
	REAL=4;
	% Constants to identify kind of factors to generate
	%> Designates a dense factor matrix
	DENSE=0;
	%> Designates a dense factor matrix
	SPARSE=1;
	%> Means DENSE or SPARSE
	MIXED=2;

	if(nargin < 2)
		error('matfaust.rand(): the number of arguments must be at least 2.')
	end
	num_rows = M;
	num_cols = N;
	% default values
	field = 'real';
	fac_type = 'sparse';
	per_row = true;
	density = -1; % default density: 5 elements per row or per column for each factor
	min_num_factors = 5;
	max_num_factors = 5;
	min_dim_size = min(num_rows, num_cols);
	max_dim_size = max(num_rows, num_cols);
	argc = length(varargin);
	dev = 'cpu';
	class = 'double';
	seed = 0;
	if(argc > 0)
		for i=1:2:argc
			if(argc > i)
				% next arg (value corresponding to the key varargin{i})
				tmparg = varargin{i+1};
			end
			switch(varargin{i})
				case 'num_factors'
					if(argc == i || ~ ismatrix(tmparg) || numel(tmparg) ~= 1 && any(size(tmparg) ~= [1 2]) || ~ any(isnumeric(tmparg)) || any(tmparg-floor(tmparg)) > 0 || any(tmparg <= 0))
						error('num_factors keyword argument is not followed by an integer or an array of positive integers')
					else
						if(isscalar(tmparg))
							min_num_factors = tmparg;
							max_num_factors = tmparg;
						else % matrix
							min_num_factors = tmparg(1);
							max_num_factors = tmparg(2);
						end
					end
				case 'dim_sizes'
					if(argc == i || ~ ismatrix(tmparg) || numel(tmparg) ~= 1 && any(size(tmparg) ~= [1 2])|| ~ any(isnumeric(tmparg)) || any(tmparg-floor(tmparg)) > 0 || any(tmparg <= 0))
						error('dim_sizes keyword argument is not followed by an integer or an array of positive integers')
					else
						if(isscalar(tmparg))
							min_dim_size = tmparg;
							max_dim_size = tmparg;
						else % matrix
							min_dim_size = tmparg(1);
							max_dim_size = tmparg(2);
						end
					end
				case 'density'
					if(argc == i || ~ isscalar(tmparg) || ~ (tmparg >= 0))
						error('density keyword argument is not followed by a positive number')
					else
						density = tmparg;
					end
				case 'fac_type'
					if(argc == i || ~ strcmp(tmparg, 'dense') && ~ strcmp(tmparg, 'sparse') && ~ strcmp(tmparg, 'mixed'))
						error('fac_type keyword argument is not followed by a valid value: dense, sparse or mixed.')
					else
						fac_type = tmparg;
					end
				case 'field'
					if(argc == i || ~ strcmp(tmparg, 'real') && ~ strcmp(tmparg, 'complex'))
						error('field keyword argument is not followed by a valid value: real, complex.')
					else
						field = tmparg;
					end
				case 'per_row'
					if(argc == i || ~ islogical(tmparg))
						error('per_row keyword argument is not followed by a logical')
					else
						per_row = tmparg;
					end
				case 'seed'
					if(argc == i || ~ isscalar(tmparg) || ~ isreal(tmparg))
						error('isscalar keyword argument is not followed by a real scalar')
					else
						seed = uint64(tmparg);
					end
				case 'dev'
					if(argc == i || ~ strcmp(tmparg, 'cpu') && ~ startsWith(tmparg, 'gpu'))
						error('dev keyword argument is not followed by a valid value: cpu, gpu*.')
					else
						dev = tmparg;
					end
				case 'class'
					if(argc == i || ~ strcmp(tmparg, 'double') && ~ strcmp(tmparg, 'single'))
						error('class keyword argument is not followed by a valid value: single, double.')
					else
						class = tmparg;
					end
				case 'seed'
				otherwise
					if((isstr(varargin{i}) || ischar(varargin{i}))  && ~ strcmp(tmparg, 'cpu') && ~ startsWith(tmparg, 'gpu') && ~ strcmp(tmparg, 'dense') && ~ strcmp(tmparg, 'sparse') && ~ strcmp(tmparg, 'mixed') && ~ strcmp(tmparg, 'real') && ~ strcmp(tmparg, 'complex'))
						error([ tmparg ' unrecognized argument'])
					end
			end
		end
	end
	% replace char array values integer values (for the C++ backend)
	if(strcmp(fac_type, 'sparse'))
		fac_type = SPARSE;
	elseif(strcmp(fac_type,'dense'))
		fac_type = DENSE;
	elseif(strcmp(fac_type,'mixed'))
		fac_type = MIXED;
	else
		error('fac_type has an unknown char array value.')
	end
	if(strcmp(field, 'real'))
		field = REAL;
	elseif(strcmp(field,'complex'))
		field = COMPLEX;
	else
		error('field has an unknown char array value.')
	end
	e = MException('FAUST:OOM', 'Out of Memory');
	if(field == COMPLEX && strcmp(class, 'single'))
		error('Complex Faust-s are not available in single precision (only double precision is possible).')
	end
	if(strcmp(dev, 'cpu'))
		if(field == COMPLEX)
			core_obj = mexFaustCplx('rand', num_rows, num_cols, fac_type, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row, seed);
			is_real = false;
		else %if(field == REAL)
			if(strcmp(class, 'double'))
				core_obj = mexFaustReal('rand', num_rows, num_cols, fac_type, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row, seed);
			else % single
				core_obj = mexFaustRealFloat('rand', num_rows, num_cols, fac_type, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row, seed);
			end
			is_real = true;
		end
	else
		if(field == COMPLEX)
			core_obj = mexFaustGPUCplx('rand', num_rows, num_cols, fac_type, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row, seed);
			is_real = false;
		else %if(field == REAL)
			if(strcmp(class, 'double'))
				core_obj = mexFaustGPUReal('rand', num_rows, num_cols, fac_type, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row, seed);
			else % float
				core_obj = mexFaustGPURealFloat('rand', num_rows, num_cols, fac_type, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row, seed);
			end
			is_real = true;
		end
	end
	if(core_obj == 0)
		throw(e)
	end
	F = matfaust.Faust(core_obj, is_real, dev, class);
end
