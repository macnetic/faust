%==========================================================================================
%> @brief Generates a random Faust composed only of BSR matrices.
%> @param M (arg. 1) The number of rows of the random Faust.
%> @param N (arg. 2) The number of columns of the random Faust.
%> @param BM (arg. 3) The nonzero block number of rows (must divide M).
%> @param BN (arg. 4)  The nonzero block number of columns (must divide N).
%> @param 'num_factors', NF (optional) If it's an integer it will be the number of random factors to set in the Faust.
%>                    If NF is a vector of 2 integers then the
%>                    number of factors will be set randomly between
%>                    NF(1) and NF(2) (inclusively).
%> @param 'density',D	(optional) the approximate density of generated factors.
%> 				   D must be a floating point number between 0 and 1.
%> @param 'field', str	(optional) str is either 'real' or 'complex' to set the Faust field.
%>                  The default value is 'real'.
%> @param 'dev', 'gpu or 'cpu' (optional) to create the random Faust on CPU or GPU (by default on CPU).
%> @param 'class', 'double' (by default) or 'single' (optional) to select the scalar type used for the Faust generated.
%>
%> @retval F the random Faust.
%>
%> @b Example @b 1
%> @code
%> >> matfaust.rand_bsr(10, 10, 2, 2)
%> ans =
%>
%> Faust size 10x10, density 0.6, nnz_sum 60, 5 factor(s):
%> - FACTOR 0 (double) BSR, size 10x10, density 0.12, nnz 12
%> - FACTOR 1 (double) BSR, size 10x10, density 0.12, nnz 12
%> - FACTOR 2 (double) BSR, size 10x10, density 0.12, nnz 12
%> - FACTOR 3 (double) BSR, size 10x10, density 0.12, nnz 12
%> - FACTOR 4 (double) BSR, size 10x10, density 0.12, nnz 12
%> @endcode
%>
%> @b Example @b 2
%> @code
%>>> matfaust.rand_bsr(128, 128, 32, 32, 'num_factors', [6, 10], 'field', 'real', 'class', 'double', 'density', .8)
%>
%>ans = 
%>
%>Faust size 128x128, density 8.125, nnz_sum 133120, 10 factor(s): 
%>- FACTOR 0 (double) BSR, size 128x128, density 0.8125, nnz 13312
%>- FACTOR 1 (double) BSR, size 128x128, density 0.8125, nnz 13312
%>- FACTOR 2 (double) BSR, size 128x128, density 0.8125, nnz 13312
%>- FACTOR 3 (double) BSR, size 128x128, density 0.8125, nnz 13312
%>- FACTOR 4 (double) BSR, size 128x128, density 0.8125, nnz 13312
%>- FACTOR 5 (double) BSR, size 128x128, density 0.8125, nnz 13312
%>- FACTOR 6 (double) BSR, size 128x128, density 0.8125, nnz 13312
%>- FACTOR 7 (double) BSR, size 128x128, density 0.8125, nnz 13312
%>- FACTOR 8 (double) BSR, size 128x128, density 0.8125, nnz 13312
%>- FACTOR 9 (double) BSR, size 128x128, density 0.8125, nnz 13312
%>
%>
%> @endcode
function F = rand_bsr(M, N, BM, BN, varargin)
	argc = length(varargin);
	dev = 'cpu';
	class = 'double';
	field = 'real';
	min_num_factors = 5;
	max_num_factors = 5;
	density = .1;
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
				case 'field'
					if(argc == i || ~ strcmp(tmparg, 'real') && ~ strcmp(tmparg, 'complex'))
						error('field keyword argument is not followed by a valid value: real, complex.')
					else
						field = tmparg;
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
				otherwise
					if((isstr(varargin{i}) || ischar(varargin{i}))  && ~ strcmp(tmparg, 'cpu') && ~ startsWith(tmparg, 'gpu') && ~ strcmp(tmparg, 'real') && ~ strcmp(tmparg, 'complex'))
						error([ tmparg ' unrecognized argument'])
					end
			end
		end
	end
	if(strcmp(dev, 'cpu'))
		if(strcmp(field, 'complex'))
			core_obj = mexFaustCplx('rand_bsr', M, N, min_num_factors, max_num_factors, BM, BN, density);
			is_real = false;
		else %if(field == REAL)
			if(strcmp(class, 'double'))
				core_obj = mexFaustReal('rand_bsr', M, N, min_num_factors, max_num_factors, BM, BN, density);
			else % float/single
				core_obj = mexFaustRealFloat('rand_bsr', M, N, min_num_factors, max_num_factors, BM, BN, density);
			end
			is_real = true;
		end
	else
		% dev is 'gpu'
		if(strcmp(field, 'complex'))
			core_obj = mexFaustGPUCplx('rand_bsr', M, N, min_num_factors, max_num_factors, BM, BN, density);
			is_real = false;
		else %if(field == REAL)
			if(strcmp(class, 'double'))
				core_obj = mexFaustGPUReal('rand_bsr', M, N, min_num_factors, max_num_factors, BM, BN, density);
			else % float/single
				core_obj = mexFaustGPURealFloat('rand_bsr', M, N, min_num_factors, max_num_factors, BM, BN, density);
			end
			is_real = true;
		end
	end
	if(core_obj == 0)
		throw(MException('FAUST:OOM', 'Out of Memory'))
	end
	if(strcmp('single', class))
		class = 'float';
	end
	F = matfaust.Faust(core_obj, is_real, dev, class);
end
