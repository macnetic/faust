% ======================================================================
%> @brief     This factory class provides methods for generating a Faust especially by factorization of a dense matrix.
%>
%>    This class gives access to the main factorization algorithms of
%>    FAµST. Those algorithms can factorize a dense matrix to a sparse product
%>    (i.e. a Faust object).
%>
%>    There are two algorithms for factorization.
%>
%>    The first one is Palm4MSA :
%>    which stands for Proximal Alternating Linearized Minimization for
%>    Multi-layer Sparse Approximation. Note that Palm4MSA is not
%>    intended to be used directly. You should rather rely on the second algorithm.
%>
%>    The second one is the Hierarchical Factorization algorithm:
%>    this is the central algorithm to factorize a dense matrix to a Faust.
%>    It makes iterative use of Palm4MSA to proceed with the factorization of a given
%>    dense matrix.
%>
%>    A more secondary functionality of this class is the pseudo-random generation of a
%>    Faust with FaustFactory.rand().
%>
% ======================================================================
classdef FaustFactory
	properties (SetAccess = public)

	end
	methods(Static)

		%==========================================================================================
		%> @brief Factorizes the matrix M with Palm4MSA algorithm using the parameters set in p.
		%>
		%>
		%> @param M the dense matrix to factorize.
		%> @param p the ParamsPalm4MSA instance to define the algorithm parameters.
		%>
		%> @retval F the Faust object result of the factorization.
		%>
		%> @b Example
		%>
		%> @code
		%>  import matfaust.*
		%>  num_facts = 2
		%>  is_update_way_R2L = false
		%>  init_lambda = 1.0
		%>  M = rand(500, 32)
		%>  cons = cell(2,1)
		%>  cons{1} = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5);
		%>  cons{2} = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32, 32, 1.0);
		%>  stop_crit = StoppingCriterion(200);
		%>  params = ParamsPalm4MSA(num_facts, is_update_way_R2L, init_lambda, cons, stop_crit);
		%>  F = FaustFactory.fact_palm4msa(M, params)
		%> @endcode
		%>
		%> F =
		%>
		%> Faust size 500x32, density 0.22025, nnz_sum 3524, 2 factor(s):
		%> - FACTOR 0 (real) DENSE, size 500x32, density 0.15625, nnz 2500
		%> - FACTOR 1 (real) DENSE, size 32x32, density 1, nnz 1024
		%>
		%>
		%==========================================================================================
		function  F = fact_palm4msa(M, p)
			import matfaust.Faust
			mex_constraints = cell(1, length(p.constraints));
			matfaust.FaustFactory.check_fact_mat('FaustFactory.fact_palm4msa', M)
			for i=1:length(p.constraints)
				cur_cell = cell(1, 4);
				cur_cell{1} = p.constraints{i}.name.conv2str();
				cur_cell{2} = p.constraints{i}.param;
				cur_cell{3} = p.constraints{i}.num_rows;
				cur_cell{4} = p.constraints{i}.num_cols;
				mex_constraints{i} = cur_cell;
			end
			% put mex_constraints in a cell array again because mex eats one level of array
			mex_params = struct('data', M, 'nfacts', p.num_facts, 'cons', {mex_constraints}, 'init_facts', {p.init_facts}, 'niter', p.stop_crit.num_its, 'sc_is_criterion_error', p.stop_crit.is_criterion_error, 'sc_error_treshold', p.stop_crit.error_treshold, 'sc_max_num_its', p.stop_crit.max_num_its);
			[lambda, cell_facts] = mexPalm4MSA(mex_params);
			cell_facts{1} = cell_facts{1}*lambda;
			F = Faust(cell_facts);
		end

		%==========================================================================================
		%> @brief Factorizes the matrix M with Hierarchical Factorization using the parameters set in p.
		%>
		%>
		%> @param M the dense matrix to factorize.
		%> @param p the ParamsHierarchicalFact instance to define the algorithm parameters.
		%>
		%> @retval F The Faust object result of the factorization.
		%>
		%> @Example
		%> @code
		%>  import matfaust.*;
		%>  num_facts = 4;
		%>  is_update_way_R2L = false;
		%>  init_lambda = 1.0;
		%>  M = rand(500, 32);
		%>  fact_cons = cell(3, 1);
		%>  res_cons = cell(3, 1);
		%>  fact_cons{1} = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5);
		%>  fact_cons{2} = ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 96);
		%>  fact_cons{3} = ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 96);
		%>  res_cons{1} = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32, 32, 1);
		%>  res_cons{2} =  ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 666);
		%>  res_cons{3} =  ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 333);
		%>  stop_crit = StoppingCriterion(200);
		%>  stop_crit2 = StoppingCriterion(200);
		%>  params = ParamsHierarchicalFact(num_facts, is_update_way_R2L, init_lambda, fact_cons, res_cons, size(M,1), size(M,2), {stop_crit, stop_crit2});
		%>  F = FaustFactory.fact_hierarchical(M, params)
		%>  @endcode
		%>  Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorisation 1/3<br/>
		%>  Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorisation 2/3<br/>
		%>  Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorisation 3/3<br/>
		%>
		%>  F = 
		%>
		%>  Faust size 500x32, density 0.189063, nnz_sum 3025, 4 factor(s): 
		%>  - FACTOR 0 (real) DENSE, size 500x32, density 0.15625, nnz 2500
		%>  - FACTOR 1 (real) DENSE, size 32x32, density 0.09375, nnz 96
		%>  - FACTOR 2 (real) DENSE, size 32x32, density 0.09375, nnz 96
		%>  - FACTOR 3 (real) DENSE, size 32x32, density 0.325195, nnz 333
		%>
		%==========================================================================================
		function F = fact_hierarchical(M, p)
			import matfaust.Faust
			mex_constraints = cell(2, p.num_facts-1);
			matfaust.FaustFactory.check_fact_mat('FaustFactory.fact_hierarchical', M)
			%mex_fact_constraints = cell(1, p.num_facts-1)
			for i=1:p.num_facts-1
				cur_cell = cell(1, 4);
				cur_cell{1} = p.constraints{i}.name.conv2str();
				cur_cell{2} = p.constraints{i}.param;
				cur_cell{3} = p.constraints{i}.num_rows;
				cur_cell{4} = p.constraints{i}.num_cols;
				%mex_fact_constraints{i} = cur_cell;
				mex_constraints{1,i} = cur_cell;
			end
			%mex_residuum_constraints = cell(1, p.num_facts-1)
			for i=1:p.num_facts-1
				cur_cell = cell(1, 4);
				cur_cell{1} = p.constraints{i+p.num_facts-1}.name.conv2str();
				cur_cell{2} = p.constraints{i+p.num_facts-1}.param;
				cur_cell{3} = p.constraints{i+p.num_facts-1}.num_rows;
				cur_cell{4} = p.constraints{i+p.num_facts-1}.num_cols;
				%mex_residuum_constraints{i} = cur_cell;
				mex_constraints{2,i} = cur_cell;
			end
			mex_params = struct('data', M, 'nfacts', p.num_facts, 'cons', {mex_constraints}, 'niter1', p.stop_crits{1}.num_its,'niter2', p.stop_crits{2}.num_its, 'sc_is_criterion_error', p.stop_crits{1}.is_criterion_error, 'sc_error_treshold', p.stop_crits{1}.error_treshold, 'sc_max_num_its', p.stop_crits{1}.max_num_its, 'sc_is_criterion_error2', p.stop_crits{2}.is_criterion_error, 'sc_error_treshold2', p.stop_crits{2}.error_treshold, 'sc_max_num_its2', p.stop_crits{2}.max_num_its, 'nrow', p.data_num_rows, 'ncol', p.data_num_cols, 'fact_side', p.is_fact_side_left);
			[lambda, cell_facts] = mexHierarchical_fact(M, mex_params);
			cell_facts{1} = cell_facts{1}*lambda;
			F = Faust(cell_facts);
		end

		%==========================================================================================
		%> @brief Generates a random Faust.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b Faust.rand(N,S) with N and S two integers, generates a Faust of N factors. All factors are square matrices of order S. The type of factors (dense or sparse) is a random choice.
		%>
		%> &nbsp;&nbsp;&nbsp; @b Faust.rand([N1,N2],S) same as above except that here the number of factors is randomly chosen between N1 and N2 inclusively.
		%>
		%> &nbsp;&nbsp;&nbsp; @b Faust.rand([N1,N2],[S1, S2]) or @b Faust.rand(N, [S1, S2]) same as above except that here the factor matrices have random sizes; the number of rows and columns are both randomly chosen between S1 and S2 inclusively.
		%>
		%> &nbsp;&nbsp;&nbsp; @b Faust.rand(N, S, D) or @b Faust.rand([N1, N2], [S1, S2], D) same as above but specifying D the approximate density of each factor.
		%>
		%> &nbsp;&nbsp;&nbsp; @b @b Faust\.rand(N, @b S, @b D, @b 'dense') or @b Faust\.rand(@b [@b  N1, @b N2], [@b S1, @b S2], @b D, @b 'dense') same as above but generating only dense matrices as factors.
		%>
		%> &nbsp;&nbsp;&nbsp; @b Faust\.rand(@b N, @b S, @b D, @b 'sparse') or @b Faust\.rand([@b N1, @b N2], [@b S1, @b S2], @b D, @b 'sparse') same as above but generating only sparse matrices as factors.
		%>
		%> &nbsp;&nbsp;&nbsp; @b Faust\.rand(@b N, @b S, @b D, @b 'sparse', @b false), @b Faust\.rand([@b N1, @b N2], [@b S1, @b S2], @b D, @b 'sparse', @b false), Faust\.rand(@b N, @b S, @b D, @b 'dense', @b false) or @b Faust\.rand([@b N1, @b N2], [@b S1, @b S2], @b D, @b 'dense', @b false) same as above but generating a complex Faust, that is, matrices defined over a the complex field.
		%>
		%>
		%>
		%>
		%>
		%>
		%>
		%> @param num_factors (varargin{1}) If it's an integer it will be the number of random factors to set in the Faust.
		%>                    If num_factors is a vector of 2 integers then the
		%>                    number of factors will be set randomly between
		%>                    num_factors(1) and num_factors(2) (inclusively).
		%> @param dim_sizes (varargin{2}) if it's an integer it will be the order of the square
		%> 					matrix factors (of size size_dims^2).
		%> 					If it's a vector of 2 integers then the
		%> 					number of rows and columns will
		%> 					be a random number between size_dims(1) and
		%> 					size_dims(2) (inclusively).
		%> @param density	(varargin{3}, optional) the approximate density of factors generated. The default value is 0.1.
		%> 					It should be a floating point number between 0 and 1.
		%> @param fac_type	(varargin{4}, optional) the type of factors. Must be
		%>                 	'sparse', 'dense' or 'mixed' if you want a mix of dense and
		%>                  sparse matrices in the generated Faust (choice's done according
		%>                  to an uniform distribution).
		%>                  The default value is 'mixed'.
		%> @param is_real	(varargin{5}, optional) a boolean set to true to generate a real Faust,
		%>                  set to false to generate a complex Faust.
		%>                  The default value is true.
		%>
		%>
		%>
		%> @retval F the random Faust.
		%>
		%> @b Example @b 1
		%> @code
		%> % in a matlab terminal
		%> >> import matfaust.Faust
		%> >> F = Faust.rand(2, 10, .5, 'mixed', false)
		%>
		%> F =
		%>
		%> Faust size 10x10, density 0.99, nnz_sum 99, 2 factor(s):
		%> - FACTOR 0 (complex) SPARSE, size 10x10, density 0.4, nnz 40
		%> - FACTOR 1 (complex) DENSE, size 10x10, density 0.59, nnz 59
		%> @endcode
		%> @b Example @b 2
		%> @code
		%> >> import matfaust.Faust
		%> >> G = Faust.rand([2, 5], [10, 20], .5, 'dense')
		%>
		%> G =
		%>
		%> Faust size 19x16, density 1.37171, nnz_sum 417, 4 factor(s):
		%> - FACTOR 0 (real) DENSE, size 19x17, density 0.49226, nnz 159
		%> - FACTOR 1 (real) DENSE, size 17x10, density 0.517647, nnz 88
		%> - FACTOR 2 (real) DENSE, size 10x13, density 0.515385, nnz 67
		%> - FACTOR 3 (real) DENSE, size 13x16, density 0.495192, nnz 103
		%>
		%> @endcode
		%>
		%> <p>@b See @b also Faust.Faust.
		%==========================================================================================
		function F = rand(varargin)
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
				error('Faust.rand(): the number of arguments must be at least 2.')
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
				error('Faust.rand(): the first argument (num_factors) must be an integer or a vector of two integers.')
			end
			% set sizes of factors
			if(isscalar(dim_sizes) && mod(dim_sizes, 1) == 0)
				min_dim_size = dim_sizes;
				max_dim_size = dim_sizes;
			elseif(ismatrix(dim_sizes) && size(dim_sizes,1) == 1 && size(dim_sizes,2) == 2)
				min_dim_size = dim_sizes(1);
				max_dim_size = dim_sizes(2);
			else
				error('Faust.rand(): the second argument (dim_sizes) must be an integer or a vector of two integers.')
			end
			field = REAL;
			fac_type = MIXED;
			if(nargin >= 3)
				if(isnumeric(varargin{3}))
					density = varargin{3};
				else
					error('Faust.rand(): the third argument (density) must be a real number in [0;1].')
				end
				% set density type of factors
				if(nargin >= 4)
					if(ischar(varargin{4}))
						if(strcmp(varargin{4}, 'sparse'))
							fac_type = SPARSE;
						elseif(strcmp(varargin{4},'dense'))
							fac_type = DENSE;
						elseif(strcmp(varargin{4},'mixed'))
							fac_type = MIXED;
						end
					else
						error('Faust.rand(): the fourth argument (fac_type) must be among a character array among ''sparse'', ''dense'' or ''mixed''.')
					end
					%set the field of factors
					if(nargin >= 5)
						if(islogical(varargin{5}))
							if(varargin{5})
								field = REAL;
							else
								field = COMPLEX;
							end
						else
							error('Faust.rand(): the fifth argument (isreal) must be a boolean.')
						end
					else
						field = REAL;
					end
				else
					fac_type = MIXED;
				end
			else
				density = .1;
			end
			e = MException('FAUST:OOM', 'Out of Memory');
			if(field == COMPLEX)
				core_obj = mexFaustCplx('rand', fac_type, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density);
				is_real = false;
			else %if(field == REAL)
				core_obj = mexFaustReal('rand', fac_type, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density);
				is_real = true;
			end
			if(core_obj == 0)
				throw(e)
			end
			F = matfaust.Faust(core_obj, is_real);
		end

	end
	methods(Access = private, Static)
		function check_fact_mat(funcname, M)
			if(~ ismatrix(M) || isscalar(M))
				error([funcname,'() 1st argument (M) must be a matrix.'])
			end
			if(~ isnumeric(M))
				error([funcname, '() 1st argument (M) must be real or complex.'])
			end
			if(~ isreal(M))
				error([funcname, '() doesn''t yet support complex matrix factorization.'])
			end
		end
	end
end
