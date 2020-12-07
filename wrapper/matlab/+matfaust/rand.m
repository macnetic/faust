%==========================================================================================
%> @brief Generates a random Faust.
%>
%> @warning if this function is imported through 'import matfaust.rand' or 'import matfaust.*' the Matlab builtin function rand will be unreachable (matfaust.rand will be called instead). So it is not advisable to import this function, rather directly call matfaust.rand without any import.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b rand(M,N) generates a M-by-N Faust object. The numbers of rows and columns of intermediary factors are all randomly chosen between M and N (included). The number of factors is 5. The factors are sparse and real. The nnz per row of factors is 5.
%>
%> &nbsp;&nbsp;&nbsp; @b rand(M, N, NF, S) with NF and S two integers, generates a M-by-N Faust of NF factors. The factor size is S for both dimensions (except the first factor number of rows and the last factor number of columns which are respectively M and N). The factors are sparse and real.
%>
%> &nbsp;&nbsp;&nbsp; @b rand(M, N, [N1, N2], S) same as above except that here the number of factors is randomly chosen between N1 and N2 inclusively.
%>
%> &nbsp;&nbsp;&nbsp; @b rand(M, N, [N1, N2], [S1, S2]) or @b rand(M, N, NF, [S1, S2]) same as above except that here the intermediary factor dimension sizes are random; the number of rows and columns are both randomly chosen between S1 and S2 inclusively.
%>
%> &nbsp;&nbsp;&nbsp; @b rand(M, N, NF, S, D) or @b rand(M, N, [N1, N2], [S1, S2], D) same as above but specifying D the approximate density of each factor.
%>
%> &nbsp;&nbsp;&nbsp; @b rand(@b M, @b N, @b NF, @b S, @b {D, @b 'per_row'}) or @b rand(@b M, @b N, [@b N1, @b N2], [@b S1, @b S2], {@b D, @b 'per_row'}) same as above but specifying D, the density of each factor per row ('per_row') or per column ('per_col').
%>
%> &nbsp;&nbsp;&nbsp; @b @b rand(@b M, @b N, NF, @b S, @b D, @b 'dense') or @b rand(@b M, @b N, @b [@b  N1, @b N2], [@b S1, @b S2], @b D, @b 'dense') same as above but generating only dense matrices as factors.
%>
%> &nbsp;&nbsp;&nbsp; @b rand(@b M, @b N, @b NF, @b S, @b D, @b 'mixed') or @b rand(@b M, @b N, [@b N1, @b N2], [@b S1, @b S2], @b D, @b 'sparse') same as above but generating either sparse or dense matrices as factors.
%>
%> &nbsp;&nbsp;&nbsp; @b rand(@b M, @b N, @b NF, @b S, @b D, @b 'sparse', @b 'complex'), @b rand(@b M, @b N,[@b N1, @b N2], [@b S1, @b S2], @b D, @b 'sparse', @b false), rand(@b M, @b N, @b NF, @b S, @b D, @b 'dense', @b 'complex') or @b rand(@b M, @b N, [@b N1, @b N2], [@b S1, @b S2], @b D, @b 'dense', @b 'complex') same as above but generating a complex Faust, that is, matrices defined over the complex field.
%>
%>
%>
%>
%>
%>
%>
%> @param num_rows (arg. 1) The number of rows of the random Faust.
%> @param num_cols (arg. 2) The number of columns of the random Faust.
%> @param num_factors (arg. 3, optional) If it's an integer it will be the number of random factors to set in the Faust.
%>                    If num_factors is a vector of 2 integers then the
%>                    number of factors will be set randomly between
%>                    num_factors(1) and num_factors(2) (inclusively).
%> @param dim_sizes (arg. 4, optional) if it's an integer it will be the order of the square
%> 					matrix factors (of size size_dims^2).
%> 					If it's a vector of 2 integers then the
%> 					number of rows and columns will
%> 					be a random number between size_dims(1) and
%> 					size_dims(2) (inclusively).
%> @param density	(arg. 5, optional) the approximate density of factors generated.
%> 					It should be a floating point number between 0 and 1.
%>					This argument can also be a cell array {D, 'per_row'} or {D, 'per_col'} to specify the density per row or per column.
%>					By default the density is set per row and is such that the Faust's factors will have 5 non-zero elements per row.
%> @param fac_type	(arg. 6 or 7, optional) the type of factors. Must be
%>                 	'sparse', 'dense' or 'mixed' if you want a mix of dense and
%>                  sparse matrices in the generated Faust (choice's done according
%>                  to an uniform distribution).
%>                  The default value is 'sparse'.
%> @param field	(arg. 6 or 7, optional) 'real' or 'complex' to set the Faust field.
%>                  The default value is 'real'.
%>
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
%> >> G = matfaust.rand(10, 10, 2, 10, .5, 'mixed', 'complex')
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
%> >> H = matfaust.rand(10, 10, [2, 5], [10, 20], .5, 'dense')
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
%> <p>@b See @b also Faust.Faust.
%==========================================================================================
function F = rand(varargin)
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
	num_rows = varargin{1};
	num_cols = varargin{2};
	field = REAL;
	fac_type = SPARSE;
	per_row = true;
	density = -1; % default density: 5 elements per row or per column for each factor
	err_dens_not_num = 'matfaust.rand(): the argument 3 (density) must be a real number in [0;1] or a cell array of length 2 with density at first and ''per_row'' or ''per_col'' char array in second cell.';
	min_num_factors = 5;
	max_num_factors = 5;
	min_dim_size = min(num_rows, num_cols);
	max_dim_size = max(num_rows, num_cols);
	if(nargin >= 3)
		% set num of factors
		num_factors = varargin{3};
		if(isscalar(num_factors) && mod(num_factors,1) == 0)
			min_num_factors = num_factors;
			max_num_factors = num_factors;
		elseif(ismatrix(num_factors) && size(num_factors, 1) == 1 && size(num_factors, 2) == 2)
			min_num_factors = num_factors(1);
			max_num_factors = num_factors(2);
		else
			error('matfaust.rand(): the argument 3 (num_factors) must be an integer or a vector of two integers.')
		end
		if(nargin >= 4)
			% set sizes of factors
			dim_sizes = varargin{4};
			if(isscalar(dim_sizes) && mod(dim_sizes, 1) == 0)
				min_dim_size = dim_sizes;
				max_dim_size = dim_sizes;
			elseif(ismatrix(dim_sizes) && size(dim_sizes,1) == 1 && size(dim_sizes,2) == 2)
				min_dim_size = dim_sizes(1);
				max_dim_size = dim_sizes(2);
			else
				error('matfaust.rand(): the argument 4 (dim_sizes) must be an integer or a vector of two integers.')
			end

			if(nargin >= 5)
				if(isscalar(varargin{5}) && isreal(varargin{5}))
					density = varargin{5};
				elseif(iscell(varargin{5}))
					density_cell = varargin{5};
					if(isreal(density_cell{1}))
						density = density_cell{1};
						if(length(density_cell) >= 2)
							if(strcmp(density_cell{2}, 'per_row'))
								per_row = true;
							elseif(strcmp(density_cell{2}, 'per_col'))
								per_row = false;
							else
								error('matfaust.rand(): when argument 5 (density) is a cell the first cell element must be a real number into [0;1] and the second a char array ''per_row'' or ''per_row''.')
							end
						end
					else
						error(err_dens_not_num)
					end
				else
					error(err_dens_not_num)
				end
				if(nargin >= 6)
					for i=6:nargin
						% set repr. type of factors ('sparse', 'dense', 'mixed') and field ('real' or 'complex')
						if(nargin >= i)
							err_unknown_arg6or5 = ['matfaust.rand(): the argument ' int2str(i) ' (fac_type) must be among a character array among ''sparse'', ''dense'', ''mixed'', ''real'', or ''complex''.'];
							if(ischar(varargin{i}))
								if(strcmp(varargin{i}, 'sparse'))
									fac_type = SPARSE;
								elseif(strcmp(varargin{i},'dense'))
									fac_type = DENSE;
								elseif(strcmp(varargin{i},'mixed'))
									fac_type = MIXED;
								elseif(strcmp(varargin{i}, 'real'))
									field = REAL;
								elseif(strcmp(varargin{i},'complex'))
									field = COMPLEX;
								else
									error(err_unknown_arg6or5)
								end
							else
								error(err_unknown_arg6or5)
							end
						end
					end
				end
			end
		end
	end
	e = MException('FAUST:OOM', 'Out of Memory');
	if(field == COMPLEX)
		core_obj = mexFaustCplx('rand', num_rows, num_cols, fac_type, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row);
		is_real = false;
	else %if(field == REAL)
		core_obj = mexFaustReal('rand', num_rows, num_cols, fac_type, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row);
		is_real = true;
	end
	if(core_obj == 0)
		throw(e)
	end
	F = matfaust.Faust(core_obj, is_real);
end
