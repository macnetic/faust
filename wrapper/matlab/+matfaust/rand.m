%==========================================================================================
%> @brief Generates a random Faust.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b rand(N,S) with N and S two integers, generates a Faust of N factors. All factors are square matrices of order S. The type of factors (dense or sparse) is a random choice.
%>
%> &nbsp;&nbsp;&nbsp; @b rand([N1,N2],S) same as above except that here the number of factors is randomly chosen between N1 and N2 inclusively.
%>
%> &nbsp;&nbsp;&nbsp; @b rand([N1,N2],[S1, S2]) or @b rand(N, [S1, S2]) same as above except that here the factor matrices have random sizes; the number of rows and columns are both randomly chosen between S1 and S2 inclusively.
%>
%> &nbsp;&nbsp;&nbsp; @b rand(N, S, D) or @b rand([N1, N2], [S1, S2], D) same as above but specifying D the approximate density of each factor.
%>
%> &nbsp;&nbsp;&nbsp; @b rand(@b N, @b S, @b {D, @b 'per_row'}) or @b rand([@b N1, @b N2], [@b S1, @b S2], {@b D, @b 'per_row'}) same as above but specifying D, the density of each factor per row ('per_row') or per column ('per_col').
%>
%> &nbsp;&nbsp;&nbsp; @b @b rand(N, @b S, @b D, @b 'dense') or @b rand(@b [@b  N1, @b N2], [@b S1, @b S2], @b D, @b 'dense') same as above but generating only dense matrices as factors.
%>
%> &nbsp;&nbsp;&nbsp; @b rand(@b N, @b S, @b D, @b 'sparse') or @b rand([@b N1, @b N2], [@b S1, @b S2], @b D, @b 'sparse') same as above but generating only sparse matrices as factors.
%>
%> &nbsp;&nbsp;&nbsp; @b rand(@b N, @b S, @b D, @b 'sparse', @b 'complex'), @b rand([@b N1, @b N2], [@b S1, @b S2], @b D, @b 'sparse', @b false), rand(@b N, @b S, @b D, @b 'dense', @b 'complex') or @b rand([@b N1, @b N2], [@b S1, @b S2], @b D, @b 'dense', @b 'complex') same as above but generating a complex Faust, that is, matrices defined over the complex field.
%>
%>
%>
%>
%>
%>
%>
%> @param num_factors (arg. 1) If it's an integer it will be the number of random factors to set in the Faust.
%>                    If num_factors is a vector of 2 integers then the
%>                    number of factors will be set randomly between
%>                    num_factors(1) and num_factors(2) (inclusively).
%> @param dim_sizes (arg. 2) if it's an integer it will be the order of the square
%> 					matrix factors (of size size_dims^2).
%> 					If it's a vector of 2 integers then the
%> 					number of rows and columns will
%> 					be a random number between size_dims(1) and
%> 					size_dims(2) (inclusively).
%> @param density	(arg. 3, optional) the approximate density of factors generated.
%> 					It should be a floating point number between 0 and 1.
%>					This argument can also be a cell array {D, 'per_row'} or {D, 'per_col'} to specify the density per row or per column.
%>					By default the density is set per row and is such that the Faust's factors will have 5 non-zero elements per row.
%> @param fac_type	(arg. 4 or 5, optional) the type of factors. Must be
%>                 	'sparse', 'dense' or 'mixed' if you want a mix of dense and
%>                  sparse matrices in the generated Faust (choice's done according
%>                  to an uniform distribution).
%>                  The default value is 'mixed'.
%> @param field	(arg. 4 or 5, optional) 'real' or 'complex' to set the Faust field.
%>                  The default value is 'real'.
%>
%>
%>
%> @retval F the random Faust.
%>
%> @b Example @b 1
%> @code
%> % in a matlab terminal
%> >> F = matfaust.rand(2, 10, .5, 'mixed', 'complex')
%>
%> F =
%>
%> Faust size 10x10, density 1, nnz_sum 100, 2 factor(s):
%> - FACTOR 0 (complex) SPARSE, size 10x10, density 0.5, nnz 50
%> - FACTOR 1 (complex) DENSE, size 10x10, density 0.5, nnz 50
%> @endcode
%> @b Example @b 2
%> @code
%> >> G = matfaust.rand([2, 5], [10, 20], .5, 'dense')
%>
%> G =
%>
%> Faust size 19x18, density 0.973684, nnz_sum 333, 3 factor(s):
%> - FACTOR 0 (real) DENSE, size 19x12, density 0.5, nnz 114
%> - FACTOR 1 (real) DENSE, size 12x15, density 0.466667, nnz 84
%> - FACTOR 2 (real) DENSE, size 15x18, density 0.5, nnz 135
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
	% set num of factors
	num_factors = varargin{1};
	dim_sizes = varargin{2};
	if(isscalar(num_factors) && mod(num_factors,1) == 0)
		min_num_factors = num_factors;
		max_num_factors = num_factors;
	elseif(ismatrix(num_factors) && size(num_factors, 1) == 1 && size(num_factors, 2) == 2)
		min_num_factors = num_factors(1);
		max_num_factors = num_factors(2);
	else
		error('matfaust.rand(): the argument 1 (num_factors) must be an integer or a vector of two integers.')
	end
	% set sizes of factors
	if(isscalar(dim_sizes) && mod(dim_sizes, 1) == 0)
		min_dim_size = dim_sizes;
		max_dim_size = dim_sizes;
	elseif(ismatrix(dim_sizes) && size(dim_sizes,1) == 1 && size(dim_sizes,2) == 2)
		min_dim_size = dim_sizes(1);
		max_dim_size = dim_sizes(2);
	else
		error('matfaust.rand(): the argument 2 (dim_sizes) must be an integer or a vector of two integers.')
	end
	field = REAL;
	fac_type = MIXED;
	per_row = true;
	density = -1; % default density: 5 elements per row or per column for each factor
	err_dens_not_num = 'matfaust.rand(): the argument 3 (density) must be a real number in [0;1] or a cell array of length 2 with density at first and ''per_row'' or ''per_col'' char array in second cell.';
	if(nargin >= 3)
		if(isscalar(varargin{3}) && isreal(varargin{3}))
			density = varargin{3};
		elseif(iscell(varargin{3}))
			density_cell = varargin{3};
			if(isreal(density_cell{1}))
				density = density_cell{1};
				if(length(density_cell) >= 2)
					if(strcmp(density_cell{2}, 'per_row'))
						per_row = true;
					elseif(strcmp(density_cell{2}, 'per_col'))
						per_row = false;
					else
						error('matfaust.rand(): when argument 3 (density) is a cell the first cell element must be a real number into [0;1] and the second a char array ''per_row'' or ''per_row''.')
					end
				end
			else
				error(err_dens_not_num)
			end
		else
			error(err_dens_not_num)
		end
		if(nargin >= 4)
			for i=4:nargin
				% set repr. type of factors ('sparse', 'dense', 'mixed') and field ('real' or 'complex')
				if(nargin >= i)
					err_unknown_arg4or5 = ['matfaust.rand(): the argument ' int2str(i) ' (fac_type) must be among a character array among ''sparse'', ''dense'', ''mixed'', ''real'', or ''complex''.'];
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
							error(err_unknown_arg4or5)
						end
					else
						error(err_unknown_arg4or5)
					end
				end
			end
		end
	end
	e = MException('FAUST:OOM', 'Out of Memory');
	if(field == COMPLEX)
		core_obj = mexFaustCplx('rand', fac_type, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row);
		is_real = false;
	else %if(field == REAL)
		core_obj = mexFaustReal('rand', fac_type, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row);
		is_real = true;
	end
	if(core_obj == 0)
		throw(e)
	end
	F = matfaust.Faust(core_obj, is_real);
end
