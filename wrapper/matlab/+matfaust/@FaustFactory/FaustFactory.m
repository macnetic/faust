%% class FaustFactory
%%

% ======================================================================
%> @brief     This factory class provides methods for generating a Faust especially by factorization of a dense matrix.
%>
%>
%>    This class gives access to the main factorization algorithms of
%>    FAµST. Those algorithms can factorize a dense matrix to a sparse product
%>    (i.e. a Faust object).
%>
%>    There are two algorithms for factorization.
%>
%>    - The first one is Palm4MSA :
%>    which stands for Proximal Alternating Linearized Minimization for
%>    Multi-layer Sparse Approximation. Note that Palm4MSA is not
%>    intended to be used directly. You should rather rely on the second algorithm.
%>
%>    - The second one is the Hierarchical Factorization algorithm:
%>    this is the central algorithm to factorize a dense matrix to a Faust.
%>    It makes iterative use of Palm4MSA to proceed with the factorization of a given
%>    dense matrix.
%>
%>    A more secondary functionality of this class is the Faust generation.
%>    Several methods are available:
%>
%>    - The pseudo-random generation of a Faust with FaustFactory.rand(),
%>    - the discrete Fourier transform with FaustFactory.dft(),
%>    - the Hadamard transform with FaustFactory.wht(),
%>    - and the identity Faust with FaustFactory.eye().
%>
% ======================================================================

classdef FaustFactory
	properties (SetAccess = public)

	end
	properties(SetAccess = private, Hidden = true, Constant)
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
		%> @retval [F, lambda] = fact_palm4msa(M, p) to optionally get lambda (scale).
		%>
		%> @b Example
		%>
		%> @code
		%>  import matfaust.*
		%>  import matfaust.factparams.*
		%>  M = rand(500, 32);
		%>  cons = ConstraintList('splin', 5, 500, 32, 'normcol', 1, 32, 32);
		%>  stop_crit = StoppingCriterion(200);
		%>  params = ParamsPalm4MSA(cons, stop_crit, 'is_update_way_R2L', false, 'init_lambda', 1.0);
		%>  F = FaustFactory.fact_palm4msa(M, params)
		%> @endcode
		%>
		%> F =
		%>
		%> Faust size 500x32, density 0.22025, nnz_sum 3524, 2 factor(s):
		%> - FACTOR 0 (real) SPARSE, size 500x32, density 0.15625, nnz 2500
		%> - FACTOR 1 (real) SPARSE, size 32x32, density 1, nnz 1024
		%>
		%>
		%==========================================================================================
		function  [F,lambda] = fact_palm4msa(M, p)
			import matfaust.Faust
			mex_constraints = cell(1, length(p.constraints));
			matfaust.FaustFactory.check_fact_mat('FaustFactory.fact_palm4msa', M)
			if(~ isa(p ,'matfaust.factparams.ParamsPalm4MSA'))
				error('p must be a ParamsPalm4MSA object.')
			end
			for i=1:length(p.constraints)
				cur_cell = cell(1, 4);
				cur_cell{1} = p.constraints{i}.name.conv2str();
				cur_cell{2} = p.constraints{i}.param;
				cur_cell{3} = p.constraints{i}.num_rows;
				cur_cell{4} = p.constraints{i}.num_cols;
				mex_constraints{i} = cur_cell;
			end
			if(~ p.is_mat_consistent(M))
				error('M''s number of columns must be consistent with the last residuum constraint defined in p. Likewise its number of rows must be consistent with the first factor constraint defined in p.')
			end
			% put mex_constraints in a cell array again because mex eats one level of array
			mex_params = struct('data', M, 'nfacts', p.num_facts, 'cons', {mex_constraints}, 'init_facts', {p.init_facts}, 'niter', p.stop_crit.num_its, 'sc_is_criterion_error', p.stop_crit.is_criterion_error, 'sc_error_treshold', p.stop_crit.error_treshold, 'sc_max_num_its', p.stop_crit.max_num_its, 'update_way', p.is_update_way_R2L);
			if(isreal(M))
				[lambda, core_obj] = mexPalm4MSAReal(mex_params);
			else
				[lambda, core_obj] = mexPalm4MSACplx(mex_params);
			end
			F = Faust(core_obj, isreal(M));
		end

		%==========================================================================================
		%> @brief Factorizes the matrix M with Hierarchical Factorization using the parameters set in p.
		%>
		%>
		%> @param M the dense matrix to factorize.
		%> @param p is a set of factorization parameters. It might be a fully defined instance of parameters (matfaust.factparams.ParamsHierarchicalFact) or a simplified expression which designates a pre-defined parametrization:
		%> - 'squaremat' to use pre-defined parameters typically used to factorize a Hadamard square matrix of order a power of two (see matfaust.demo.hadamard).
		%> - {'rectmat', j, k, s} to use pre-defined parameters used for instance in factorization of the MEG matrix which is a rectangular matrix of size m*n such that m < n (see matfaust.demo.bsl); j is the number of factors, k the sparsity of the main factor's columns, and s the sparsity of rows for all other factors except the residuum (that is the first factor here because the factorization is made toward the left -- is_side_fact_left == true, cf. matfaust.factparams.ParamsHierarchicalFact).
		%> </br>The residuum has a sparsity of P*rho^(num_facts-1). <br/> By default, rho == .8 and P = 1.4. It's possible to set custom values with for example p == { 'rectmat', j, k, s, 'rho', .4, 'P', .7}. <br/>The sparsity is here the number of non-zero elements.
		%> @note - The fully defined parameters (ParamsHierarchicalFact instance) used/generated by the function are available in the return result (so one can consult what precisely mean the simplified parameterizations and possibly adjust the attributes to factorize again).
		%> @note - This function has its shorthand matfaust.faust_fact(). For convenience you might use it like this:
		%> @code
		%> import matfaust.*
		%> F = faust_fact(M, p) % equiv. to FaustFactory.fact_hierarchical(M, p)
		%> @endcode
		%>
		%> @retval F The Faust object result of the factorization.
		%> @retval [F, lambda, p_obj] = fact_hierarchical(M, p) to optionally get lambda (scale) and the p_obj ParamsHierarchicalFact instance used to factorize.
		%>
		%> @b Example 1: Fully Defined Parameters for a Random Matrix Factorization
		%> @code
		%>  import matfaust.*
		%>  import matfaust.factparams.*
		%>  M = rand(500, 32);
		%>  fact_cons = ConstraintList('splin', 5, 500, 32, 'sp', 96, 32, 32, 'sp', 96, 32, 32);
		%>  res_cons = ConstraintList('normcol', 1, 32, 32, 'sp', 666, 32, 32, 'sp', 333, 32, 32);
		%>  stop_crit = StoppingCriterion(200);
		%>  stop_crit2 = StoppingCriterion(200);
		%>  params = ParamsHierarchicalFact(fact_cons, res_cons, stop_crit, stop_crit2, 'is_update_way_R2L', false, 'init_lambda', 1.0);
		%>  F = FaustFactory.fact_hierarchical(M, params)
		%>  @endcode
		%>  Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorization 1/3<br/>
		%>  Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorization 2/3<br/>
		%>  Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorization 3/3<br/>
		%>
		%>  F = 
		%>
		%>  Faust size 500x32, density 0.189063, nnz_sum 3025, 4 factor(s): 
		%>  - FACTOR 0 (real) SPARSE, size 500x32, density 0.15625, nnz 2500
		%>  - FACTOR 1 (real) SPARSE, size 32x32, density 0.09375, nnz 96
		%>  - FACTOR 2 (real) SPARSE, size 32x32, density 0.09375, nnz 96
		%>  - FACTOR 3 (real) SPARSE, size 32x32, density 0.325195, nnz 333
		%>
		%>  @b Example 2: Simplified Parameters for Hadamard Factorization
		%>@code
		%> import matfaust.*
		%> import matfaust.FaustFactory.*
		%> % generate a Hadamard Faust of size 32x32
		%> FH = wht(5);
		%> H = full(FH); % the full matrix version
		%> % factorize it
		%> FH2 = FaustFactory.fact_hierarchical(H, 'squaremat');
		%> % test the relative error
		%> norm(FH-FH2, 'fro')/norm(FH, 'fro') % the result is 1.1015e-16, the factorization is accurate
		%>@endcode
		%>FH =
		%>
		%>Faust size 32x32, density 0.3125, nnz_sum 320, 5 factor(s):
		%>- FACTOR 0 (real) SPARSE, size 32x32, density 0.0625, nnz 64
		%>- FACTOR 1 (real) SPARSE, size 32x32, density 0.0625, nnz 64
		%>- FACTOR 2 (real) SPARSE, size 32x32, density 0.0625, nnz 64
		%>- FACTOR 3 (real) SPARSE, size 32x32, density 0.0625, nnz 64
		%>- FACTOR 4 (real) SPARSE, size 32x32, density 0.0625, nnz 64
		%>
		%>FH2 =
		%>
		%>Faust size 32x32, density 0.3125, nnz_sum 320, 5 factor(s):
		%>- FACTOR 0 (real) SPARSE, size 32x32, density 0.0625, nnz 64
		%>- FACTOR 1 (real) SPARSE, size 32x32, density 0.0625, nnz 64
		%>- FACTOR 2 (real) SPARSE, size 32x32, density 0.0625, nnz 64
		%>- FACTOR 3 (real) SPARSE, size 32x32, density 0.0625, nnz 64
		%>- FACTOR 4 (real) SPARSE, size 32x32, density 0.0625, nnz 64
		%>
		%>
		%>@b Example 3: Simplified Parameters for a Rectangular Matrix Factorization (the BSL demo  MEG matrix)
		%>
		%> @code
		%> >> % in a matlab terminal
		%> >> import matfaust.*
		%> >> load('matrix_MEG.mat')
		%> >> MEG = matrix;
		%> >> num_facts = 9;
		%> >> k = 10;
		%> >> s = 8;
		%> >> MEG16 = FaustFactory.fact_hierarchical(MEG, {'rectmat', num_facts, k, s})
		%> @endcode
		%> MEG16 =
		%>
		%> Faust size 204x8193, density 0.0631655, nnz_sum 105573, 9 factor(s):
		%> - FACTOR 0 (real) SPARSE, size 204x204, density 0.293613, nnz 12219
		%> - FACTOR 1 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
		%> - FACTOR 2 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
		%> - FACTOR 3 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
		%> - FACTOR 4 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
		%> - FACTOR 5 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
		%> - FACTOR 6 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
		%> - FACTOR 7 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
		%> - FACTOR 8 (real) SPARSE, size 204x8193, density 0.0490196, nnz 81930
		%>
		%> @code
		%> >> % verify the constraint k == 10, on column 5
		%> >> fac9 = get_factor(MEG16,9);
		%> >> numel(nonzeros(fac9(:,5)))
		%> @endcode
		%> ans =
		%>
		%> 10
		%>
		%> @code
		%> >> % now verify the s constraint is respected on MEG16 factor 2
		%> >> numel(nonzeros(get_factor(MEG16, 2)))/size(MEG16,1)
		%> @endcode
		%>
		%>ans =
		%>
		%> 8
		%> <p> @b See @b also matfaust.faust_fact, factparams.ParamsHierarchicalFact, factparams.ParamsHierarchicalFactSquareMat, factparams.ParamsHierarchicalFactRectMat
		%==========================================================================================
		function varargout = fact_hierarchical(M, p)
			import matfaust.Faust
			import matfaust.factparams.*
			matfaust.FaustFactory.check_fact_mat('FaustFactory.fact_hierarchical', M)
			if(~ isa(p, 'ParamsHierarchicalFact') && ParamsFactFactory.is_a_valid_simplification(p))
				p = ParamsFactFactory.createParams(M, p);
			end
			mex_constraints = cell(2, p.num_facts-1);
			if(~ isa(p ,'ParamsHierarchicalFact'))
				error('p must be a ParamsHierarchicalFact object.')
			end
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
			if(~ p.is_mat_consistent(M))
				error('M''s number of columns must be consistent with the last residuum constraint defined in p. Likewise its number of rows must be consistent with the first factor constraint defined in p.')
			end
			% the setters for num_rows/cols verifies consistency with constraints
			mex_params = struct('nfacts', p.num_facts, 'cons', {mex_constraints}, 'niter1', p.stop_crits{1}.num_its,'niter2', p.stop_crits{2}.num_its, 'sc_is_criterion_error', p.stop_crits{1}.is_criterion_error, 'sc_error_treshold', p.stop_crits{1}.error_treshold, 'sc_max_num_its', p.stop_crits{1}.max_num_its, 'sc_is_criterion_error2', p.stop_crits{2}.is_criterion_error, 'sc_error_treshold2', p.stop_crits{2}.error_treshold, 'sc_max_num_its2', p.stop_crits{2}.max_num_its, 'nrow', p.data_num_rows, 'ncol', p.data_num_cols, 'fact_side', p.is_fact_side_left, 'update_way', p.is_update_way_R2L, 'verbose', p.is_verbose, 'init_lambda', p.init_lambda);
			if(isreal(M))
				[lambda, core_obj] = mexHierarchical_factReal(M, mex_params);
			else
				[lambda, core_obj] = mexHierarchical_factCplx(M, mex_params);
			end
			F = Faust(core_obj, isreal(M));
			varargout = {F, lambda, p};
		end

		function varargout = fgft_palm(U, Lap, p, varargin)
			import matfaust.Faust
			import matfaust.factparams.*
			% TODO: check U, Lap sizes, same field
			% TODO: refactor with fact_hierarchical
			if(length(varargin) == 1)
				init_D = varargin{1};
				if(~ ismatrix(init_D) || ~ isnumeric(init_D))
					error('fgft_palm arg. 4 must be a matrix')
				end
			elseif(length(varargin) > 1)
				error('fgft_palm, too many arguments.')
			else % nargin == 0
				init_D = ones(size(U,1));
				if(~ isreal(U))
					init_D = complex(init_D);
				end
			end
			matfaust.FaustFactory.check_fact_mat('FaustFactory.fgft_palm', U)
			if(~ isa(p, 'ParamsHierarchicalFact') && ParamsFactFactory.is_a_valid_simplification(p))
				p = ParamsFactFactory.createParams(U, p);
			end
			mex_constraints = cell(2, p.num_facts-1);
			if(~ isa(p ,'ParamsHierarchicalFact'))
				error('p must be a ParamsHierarchicalFact object.')
			end
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
			if(~ p.is_mat_consistent(U))
				error('U''s number of columns must be consistent with the last residuum constraint defined in p. Likewise its number of rows must be consistent with the first factor constraint defined in p.')
			end
			% the setters for num_rows/cols verifies consistency with constraints
			mex_params = struct('nfacts', p.num_facts, 'cons', {mex_constraints}, 'niter1', p.stop_crits{1}.num_its,'niter2', p.stop_crits{2}.num_its, 'sc_is_criterion_error', p.stop_crits{1}.is_criterion_error, 'sc_error_treshold', p.stop_crits{1}.error_treshold, 'sc_max_num_its', p.stop_crits{1}.max_num_its, 'sc_is_criterion_error2', p.stop_crits{2}.is_criterion_error, 'sc_error_treshold2', p.stop_crits{2}.error_treshold, 'sc_max_num_its2', p.stop_crits{2}.max_num_its, 'nrow', p.data_num_rows, 'ncol', p.data_num_cols, 'fact_side', p.is_fact_side_left, 'update_way', p.is_update_way_R2L, 'init_D', init_D, 'verbose', p.is_verbose, 'init_lambda', p.init_lambda);
			if(isreal(U))
				[lambda, core_obj] = mexHierarchical_factReal(U, mex_params, Lap);
			else
				[lambda, core_obj] = mexHierarchical_factCplx(U, mex_params, Lap);
			end
			F = Faust(core_obj, isreal(U));
			varargout = {F, lambda, p};
		end

		%==========================================================================================
		%> @brief Constructs a Faust implementing the Walsh-Hadamard Transform of order 2^n.
		%>
		%> The resulting Faust has n sparse factors of order 2^n, each one having 2 non-zero elements
		%> per row and per column.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b H = wht(n)
		%>
		%> @param n the power of two exponent for a Hadamard matrix of order 2^n and a factorization into n factors.
		%>
		%> @retval H the Faust implementing the Hadamard transform of dimension 2^n.
		%>
		%> @b Example
		%> @code
		%> % in a matlab terminal
		%> >> import matfaust.FaustFactory.*
		%> >> H = wht(10)
		%> @endcode
		%>
		%>
		%>H =
		%>
		%>Faust size 1024x1024, density 0.0195312, nnz_sum 20480, 10 factor(s):
		%>- FACTOR 0 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%>- FACTOR 1 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%>- FACTOR 2 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%>- FACTOR 3 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%>- FACTOR 4 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%>- FACTOR 5 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%>- FACTOR 6 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%>- FACTOR 7 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%>- FACTOR 8 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%>- FACTOR 9 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%>
		%==========================================================================================
		function H = wht(n)
			% check n (must be integer > 0)
			if(~ isreal(n) || n < 0 || abs(n-floor(n)) > 0)
				error('n must be an integer greater than zero')
			end
			if(n>31)
				error('Can''t handle a Hadamard Faust of order larger than 2^31')
			end
			core_obj = mexFaustReal('hadamard', n);
			is_real = true;
			e = MException('FAUST:OOM', 'Out of Memory');
			if(core_obj == 0)
				throw(e)
			end
			H = matfaust.Faust(core_obj, is_real);
		end

		%==========================================================================================
		%> @brief Constructs a Faust whose the full matrix is the Discrete Fourier Transform square matrix of order 2^n.
		%>
		%> The factorization algorithm used is Cooley-Tukey (FFT).
		%>
		%> The resulting Faust is complex and has (n+1) sparse factors whose the n first
		%> has 2 non-zero elements per row and per column. The last factor is
		%> a permutation matrix.
		%>
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b H = dft(n)
		%>
		%> @param n: the power of two exponent for a FFT of order 2^n and a
		%> factorization in n+1 factors.
		%>
		%>
		%> @retval F the Faust implementing the FFT transform of dimension 2^n.
		%>
		%> @b Example
		%> @code
		%> % in a matlab terminal
		%> >> import matfaust.FaustFactory.*
		%> >> F = dft(10)
		%> @endcode
		%>
		%>
		%> F =
		%>
		%> Faust size 1024x1024, density 0.0205078, nnz_sum 21504, 11 factor(s):
		%> - FACTOR 0 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%> - FACTOR 1 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%> - FACTOR 2 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%> - FACTOR 3 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%> - FACTOR 4 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%> - FACTOR 5 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%> - FACTOR 6 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%> - FACTOR 7 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%> - FACTOR 8 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%> - FACTOR 9 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
		%> - FACTOR 10 (complex) SPARSE, size 1024x1024, density 0.000976562, nnz 1024
		%==========================================================================================
		function H = dft(n)
			% check n (must be integer > 0)
			if(~ isreal(n) || n < 0 || abs(n-floor(n)) > 0)
				error('n must be an integer greater than zero')
			end
			if(n>31)
				error('Can''t handle a FFT Faust of order larger than 2^31')
			end
			core_obj = mexFaustCplx('fourier', n);
			is_real = false;
			e = MException('FAUST:OOM', 'Out of Memory');
			if(core_obj == 0)
				throw(e)
			end
			H = matfaust.Faust(core_obj, is_real);
		end

		%==========================================================================================
		%> @brief Faust identity matrix.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b FaustFactory.eye(m,n) or FaustFactory.eye([m,n]) forms a M-by-N Faust F = Faust(speye(M,N)).<br/>
		%> &nbsp;&nbsp;&nbsp; @b FaustFactory.eye(m) is a short for FaustFactory.eye(m,n).<br/>
		%> &nbsp;&nbsp;&nbsp; @b FaustFactory.eye(S, 'complex') or FaustFactory.eye(S, 'complex') or FaustFactory.eye(S, 'complex') with S the size, does the same as above but returns a complex Faust.</br>
		%>
		%> @b Example
		%> @code
		%> % in a matlab terminal
		%>>> import matfaust.FaustFactory
		%>>> FaustFactory.eye(4)
		%>
		%>ans =
		%>
		%>Faust size 4x4, density 0.25, nnz_sum 4, 1 factor(s):
		%>- FACTOR 0 (real) SPARSE, size 4x4, density 0.25, nnz 4
		%>
		%>>> full(FaustFactory.eye(4))
		%>
		%>ans =
		%>
		%>     1     0     0     0
		%>     0     1     0     0
		%>     0     0     1     0
		%>     0     0     0     1
		%>
		%>>> full(FaustFactory.eye(4,5))
		%>
		%>ans =
		%>
		%>     1     0     0     0     0
		%>     0     1     0     0     0
		%>     0     0     1     0     0
		%>     0     0     0     1     0
		%>
		%>>> full(FaustFactory.eye([5,4]))
		%>
		%>ans =
		%>
		%>     1     0     0     0
		%>     0     1     0     0
		%>     0     0     1     0
		%>     0     0     0     1
		%>     0     0     0     0
		%>
		%>>> full(FaustFactory.eye([5,4],'complex'))
		%>
		%>ans =
		%>
		%>   1.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
		%>   0.0000 + 0.0000i   1.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
		%>   0.0000 + 0.0000i   0.0000 + 0.0000i   1.0000 + 0.0000i   0.0000 + 0.0000i
		%>   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   1.0000 + 0.0000i
		%>   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
		%>
		%>>> full(FaustFactory.eye([4],'complex'))
		%>
		%>ans =
		%>
		%>   1.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
		%>   0.0000 + 0.0000i   1.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
		%>   0.0000 + 0.0000i   0.0000 + 0.0000i   1.0000 + 0.0000i   0.0000 + 0.0000i
		%>   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   1.0000 + 0.0000i
		%> @endcode
		%>
		%==========================================================================================
		function F = eye(varargin)
			if(nargin < 1)
				error('First argument is mandatory')
			end
			import matfaust.Faust
			if(ismatrix(varargin{1}))
				shape = varargin{1};
				ndim = size(shape,2);
				nrows = size(shape,1);
				if(ndim > 2)
					error('N-dimensional arrays are not supported.')
				elseif(nrows > 1)
					error('Size vector should be a row vector with real elements.')
				elseif(ndim == 2)
					m = shape(1);
					n = shape(2);
				elseif(ndim == 1)
					m = varargin{1};
					if(nargin > 1 && isnumeric(varargin{2}))
						n = varargin{2};
					else
						n = m;
					end
				else
					error('Size vector should be a row vector with real elements.')
				end
			else
				error('Size inputs must be numeric.')
			end
			la = varargin{nargin};
			if(nargin ~= 1 && ~ isnumeric(la) && (ischar(la) || ischar(cell2mat(la))))
				% F = Faust(sparse(1:m, 1:n, 1+eps(1)*j)); % hack to avoid passing through a full matrix
				if(strcmp(la,'complex'))
					F = Faust(eye(m,n,'like', sparse(1,1,1+i)));
				elseif(strcmp(la, 'real'))
					F = Faust(speye(m,n));
				else
					if(iscell(la))
						la = cell2mat(la)
					end
					error(['Unknown option: ' la])
				end
			else
				F = Faust(speye(m,n));
			end
		end


		%==========================================================================================
		%> @brief Generates a random Faust.
		%>
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b FaustFactory.rand(N,S) with N and S two integers, generates a Faust of N factors. All factors are square matrices of order S. The type of factors (dense or sparse) is a random choice.
		%>
		%> &nbsp;&nbsp;&nbsp; @b FaustFactory.rand([N1,N2],S) same as above except that here the number of factors is randomly chosen between N1 and N2 inclusively.
		%>
		%> &nbsp;&nbsp;&nbsp; @b FaustFactory.rand([N1,N2],[S1, S2]) or @b FaustFactory.rand(N, [S1, S2]) same as above except that here the factor matrices have random sizes; the number of rows and columns are both randomly chosen between S1 and S2 inclusively.
		%>
		%> &nbsp;&nbsp;&nbsp; @b FaustFactory.rand(N, S, D) or @b FaustFactory.rand([N1, N2], [S1, S2], D) same as above but specifying D the approximate density of each factor.
		%>
		%> &nbsp;&nbsp;&nbsp; @b FaustFactory.rand(@b N, @b S, @b {D, @b 'per_row'}) or @b FaustFactory.rand([@b N1, @b N2], [@b S1, @b S2], {@b D, @b 'per_row'}) same as above but specifying D, the density of each factor per row ('per_row') or per column ('per_col').
		%>
		%> &nbsp;&nbsp;&nbsp; @b @b FaustFactory\.rand(N, @b S, @b D, @b 'dense') or @b FaustFactory\.rand(@b [@b  N1, @b N2], [@b S1, @b S2], @b D, @b 'dense') same as above but generating only dense matrices as factors.
		%>
		%> &nbsp;&nbsp;&nbsp; @b FaustFactory\.rand(@b N, @b S, @b D, @b 'sparse') or @b FaustFactory\.rand([@b N1, @b N2], [@b S1, @b S2], @b D, @b 'sparse') same as above but generating only sparse matrices as factors.
		%>
		%> &nbsp;&nbsp;&nbsp; @b FaustFactory\.rand(@b N, @b S, @b D, @b 'sparse', @b 'complex'), @b FaustFactory\.rand([@b N1, @b N2], [@b S1, @b S2], @b D, @b 'sparse', @b false), FaustFactory\.rand(@b N, @b S, @b D, @b 'dense', @b 'complex') or @b FaustFactory\.rand([@b N1, @b N2], [@b S1, @b S2], @b D, @b 'dense', @b 'complex') same as above but generating a complex Faust, that is, matrices defined over the complex field.
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
		%> >> import matfaust.FaustFactory
		%> >> F = FaustFactory.rand(2, 10, .5, 'mixed', 'complex')
		%>
		%> F =
		%>
		%> Faust size 10x10, density 1, nnz_sum 100, 2 factor(s):
		%> - FACTOR 0 (complex) SPARSE, size 10x10, density 0.5, nnz 50
		%> - FACTOR 1 (complex) DENSE, size 10x10, density 0.5, nnz 50
		%> @endcode
		%> @b Example @b 2
		%> @code
		%> >> import matfaust.FaustFactory
		%> >> G = FaustFactory.rand([2, 5], [10, 20], .5, 'dense')
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
			import matfaust.FaustFactory
			if(nargin < 2)
				error('FaustFactory.rand(): the number of arguments must be at least 2.')
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
				error('FaustFactory.rand(): the argument 1 (num_factors) must be an integer or a vector of two integers.')
			end
			% set sizes of factors
			if(isscalar(dim_sizes) && mod(dim_sizes, 1) == 0)
				min_dim_size = dim_sizes;
				max_dim_size = dim_sizes;
			elseif(ismatrix(dim_sizes) && size(dim_sizes,1) == 1 && size(dim_sizes,2) == 2)
				min_dim_size = dim_sizes(1);
				max_dim_size = dim_sizes(2);
			else
				error('FaustFactory.rand(): the argument 2 (dim_sizes) must be an integer or a vector of two integers.')
			end
			field = FaustFactory.REAL;
			fac_type = FaustFactory.MIXED;
			per_row = true;
			density = -1; % default density: 5 elements per row or per column for each factor
			err_dens_not_num = 'FaustFactory.rand(): the argument 3 (density) must be a real number in [0;1] or a cell array of length 2 with density at first and ''per_row'' or ''per_col'' char array in second cell.';
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
								error('FaustFactory.rand(): when argument 3 (density) is a cell the first cell element must be a real number into [0;1] and the second a char array ''per_row'' or ''per_row''.')
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
							err_unknown_arg4or5 = ['FaustFactory.rand(): the argument ' int2str(i) ' (fac_type) must be among a character array among ''sparse'', ''dense'', ''mixed'', ''real'', or ''complex''.'];
							if(ischar(varargin{i}))
								if(strcmp(varargin{i}, 'sparse'))
									fac_type = FaustFactory.SPARSE;
								elseif(strcmp(varargin{i},'dense'))
									fac_type = FaustFactory.DENSE;
								elseif(strcmp(varargin{i},'mixed'))
									fac_type = FaustFactory.MIXED;
								elseif(strcmp(varargin{i}, 'real'))
									field = FaustFactory.REAL;
								elseif(strcmp(varargin{i},'complex'))
									field = FaustFactory.COMPLEX;
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
			if(field == FaustFactory.COMPLEX)
				core_obj = mexFaustCplx('rand', fac_type, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row);
				is_real = false;
			else %if(field == FaustFactory.REAL)
				core_obj = mexFaustReal('rand', fac_type, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row);
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
%			if(~ isreal(M))
%				error([funcname, '() doesn''t yet support complex matrix factorization.'])
%			end
		end
	end
end
