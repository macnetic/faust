%==========================================================================================
%> @brief Factorizes the matrix M with Hierarchical Factorization using the parameters set in p.
%>
%>
%> @param M the dense matrix to factorize.
%> @param p is a set of factorization parameters. It might be a fully defined instance of parameters (matfaust.factparams.ParamsHierarchical) or a simplified expression which designates a pre-defined parametrization:
%> - 'squaremat' to use pre-defined parameters typically used to factorize a Hadamard square matrix of order a power of two (see matfaust.demo.hadamard).
%> - {'rectmat', j, k, s} to use pre-defined parameters used for instance in factorization of the MEG matrix which is a rectangular matrix of size m*n such that m < n (see matfaust.demo.bsl); j is the number of factors, k the sparsity of the main factor's columns, and s the sparsity of rows for all other factors except the residuum (that is the first factor here because the factorization is made toward the left -- is_side_fact_left == true, cf. matfaust.factparams.ParamsHierarchical).
%> </br>The residuum has a sparsity of P*rho^(num_facts-1). <br/> By default, rho == .8 and P = 1.4. It's possible to set custom values with for example p == { 'rectmat', j, k, s, 'rho', .4, 'P', .7}. <br/>The sparsity is here the number of non-zero elements.
%> @param 'backend',int (optional) the backend (the C++ implementation) chosen. Must be 2016 (the default) or 2020 (which should be quicker for certain configurations - e.g. factorizing a Hadamard matrix).
%> @note - The fully defined parameters (ParamsHierarchical instance) used/generated by the function are available in the return result (so one can consult what precisely mean the simplified parameterizations and possibly adjust the attributes to factorize again).
%> @note - This function has its shorthand matfaust.faust_fact(). For convenience you might use it like this:
%> @code
%> import matfaust.*
%> F = faust_fact(M, p) % equiv. to matfaust.fact.hierarchical(M, p)
%> @endcode
%>
%> @retval F The Faust object result of the factorization.
%> @retval [F, lambda, p_obj] = hierarchical(M, p) to optionally get lambda (scale) and the p_obj ParamsHierarchical instance used to factorize.
%>
%> @b Example 1: Fully Defined Parameters for a Random Matrix Factorization
%> @code
%>  import matfaust.factparams.*
%>  import matfaust.fact.hierarchical
%>  M = rand(500, 32);
%>  fact_cons = ConstraintList('splin', 5, 500, 32, 'sp', 96, 32, 32, 'sp', 96, 32, 32);
%>  res_cons = ConstraintList('normcol', 1, 32, 32, 'sp', 666, 32, 32, 'sp', 333, 32, 32);
%> % or alternatively you can use projectors-functors
%> % import matfaust.proj.*
%> % fact_cons = {splin([500,32], 5), sp([32,32], 96), sp([32,32], 96)}
%> % res_cons = {normcol([32,32], 1), sp([32,32], 666), sp([32,32], 333)}
%>  stop_crit = StoppingCriterion(200);
%>  stop_crit2 = StoppingCriterion(200);
%>  params = ParamsHierarchical(fact_cons, res_cons, stop_crit, stop_crit2, 'is_update_way_R2L', false, 'init_lambda', 1.0);
%>  F = hierarchical(M, params, 'backend', 2016)
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
%> import matfaust.fact.hierarchical
%> % generate a Hadamard Faust of size 32x32
%> FH = wht(32);
%> H = full(FH); % the full matrix version
%> % factorize it
%> FH2 = hierarchical(H, 'squaremat');
%> % test the relative error
%> norm(FH-FH2, 'fro')/norm(FH, 'fro') % the result is about 1e-16, the factorization is accurate
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
%> >> import matfaust.fact.hierarchical
%> >> load('matrix_MEG.mat')
%> >> MEG = matrix;
%> >> num_facts = 9;
%> >> k = 10;
%> >> s = 8;
%> >> MEG16 = hierarchical(MEG, {'rectmat', num_facts, k, s})
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
%> >> fac9 = factors(MEG16,9);
%> >> numel(nonzeros(fac9(:,5)))
%> @endcode
%> ans =
%>
%> 10
%>
%> @code
%> >> % now verify the s constraint is respected on MEG16 factor 2
%> >> numel(nonzeros(factors(MEG16, 2)))/size(MEG16,1)
%> @endcode
%>
%>ans =
%>
%> 8
%> <p> @b See @b also matfaust.faust_fact, factparams.ParamsHierarchical, factparams.ParamsHierarchicalSquareMat, factparams.ParamsHierarchicalRectMat
%==========================================================================================
function varargout = hierarchical(M, p, varargin)
	import matfaust.Faust
	import matfaust.factparams.*
	import matfaust.fact.check_fact_mat
	check_fact_mat('matfaust.fact.hierarchical', M)
	if(~ isa(p, 'ParamsHierarchical') && ParamsFactFactory.is_a_valid_simplification(p))
		p = ParamsFactFactory.createParams(M, p);
	end
	if(~ isa(p ,'ParamsHierarchical'))
		error('p must be a ParamsHierarchical object.')
	end
	if(~ p.is_mat_consistent(M))
		error('M''s number of columns must be consistent with the last residuum constraint defined in p. Likewise its number of rows must be consistent with the first factor constraint defined in p.')
	end
	mex_params = p.to_mex_struct();
	backend = 2016;
	nargin = length(varargin);
	if(nargin > 0)
		backend = varargin{1};
		if(strcmp('backend', backend))
			if(nargin < 2)
				error('keyword argument ''backend'' must be followed by 2016 or 2020')
			else
				backend = varargin{2};
			end
		end
		if(~ (isscalar(backend) && floor(backend) == backend) || backend ~= 2016 && backend ~= 2020)
			backend
			error('backend must be a int equal to 2016 or 2020')
		end
	end
	if(backend == 2016)
		if(isreal(M))
			[lambda, core_obj] = mexHierarchical_factReal(M, mex_params);
		else
			[lambda, core_obj] = mexHierarchical_factCplx(M, mex_params);
		end
		F = Faust(core_obj, isreal(M));
	elseif(backend == 2020)
		if(isreal(M))
			[lambda, core_obj] = mexHierarchical2020Real(M, mex_params);
		else
			error('backend 2020 doesn''t handle yet the complex matrices')
		end
		F = Faust(core_obj, isreal(M));
	end
	varargout = {F, lambda, p};
end
