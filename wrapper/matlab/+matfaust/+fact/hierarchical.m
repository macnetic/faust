%==========================================================================================
%> @brief Factorizes the matrix M with Hierarchical Factorization using the parameters set in p.
%>
%>
%> @param M the dense matrix to factorize. The class(M) value can be double or single (only if isreal(M) is true). This class has a major impact on performance.
%> @param p is a set of factorization parameters. It might be a fully defined instance of parameters (matfaust.factparams.ParamsHierarchical) or a simplified expression which designates a pre-defined parametrization:
%> - 'squaremat' to use pre-defined parameters typically used to factorize a Hadamard square matrix of order a power of two (see matfaust.demo.hadamard).
%> - {'rectmat', j, k, s} to use pre-defined parameters used for instance in factorization of the MEG matrix which is a rectangular matrix of size m*n such that m < n (see matfaust.demo.bsl); j is the number of factors, k the sparsity of the main factor's columns, and s the sparsity of rows for all other factors except the residuum (that is the first factor here because the factorization is made toward the left -- is_side_fact_left == true, cf. matfaust.factparams.ParamsHierarchical).
%> </br>The residuum has a sparsity of P*rho^(num_facts-1). <br/> By default, rho == .8 and P = 1.4. It's possible to set custom values with for example p == { 'rectmat', j, k, s, 'rho', .4, 'P', .7}. <br/>The sparsity is here the number of non-zero elements.
%> @param 'backend',int (optional) the backend to use (the C++ implementation). Must be 2016 (the default) or 2020 (which should be faster for most of the factorizations).
%> @param 'gpu', bool (optional) set to true to execute the algorithm using the GPU implementation. This option is only available when backend==2020.
%>
%>
%> @note - If backend parameter is 2020 and regardless to the StoppingCriterion-s defined in p,
%> it is possible to stop any internal call to PALM4MSA manually at any iteration
%> by the key combination CTRL-C.
%> The last Faust computed in the PALM4MSA instance will be used to continue
%> the hierarchical factorization.
%> A typical use case is when the verbose mode is enabled and you see that the error
%> doesn't change anymore or only slightly, you might stop iterations by typing CTRL-C.
%>
%> @note - The fully defined parameters (ParamsHierarchical instance) used/generated by the function are available in the return result (so one can consult what precisely mean the simplified parameterizations and possibly adjust the attributes to factorize again).
%> @note - This function has its shorthand matfaust.faust_fact(). For convenience you might use it like this:
%> @code
%> import matfaust.*
%> F = faust_fact(M, p) % equiv. to matfaust.fact.hierarchical(M, p)
%> @endcode
%>
%> @retval F The Faust object result of the factorization.
%> @retval [F, lambda] = palm4msa(M, p) when optionally getting lambda (scale).
%>
%> @b Example 1: Fully Defined Parameters for a Random Matrix Factorization
%> @code
%> >> import matfaust.factparams.*
%> >> import matfaust.fact.hierarchical
%> >> M = rand(500, 32);
%> >> fact_cons = ConstraintList('splin', 5, 500, 32, 'sp', 96, 32, 32, 'sp', 96, 32, 32);
%> >> res_cons = ConstraintList('normcol', 1, 32, 32, 'sp', 666, 32, 32, 'sp', 333, 32, 32);
%> >> % or alternatively you can use projectors-functors
%> >> % import matfaust.proj.*
%> >> % fact_cons = {splin([500,32], 5), sp([32,32], 96), sp([32,32], 96)}
%> >> % res_cons = {normcol([32,32], 1), sp([32,32], 666), sp([32,32], 333)}
%> >> stop_crit = StoppingCriterion(200);
%> >> stop_crit2 = StoppingCriterion(200);
%> >> params = ParamsHierarchical(fact_cons, res_cons, stop_crit, stop_crit2, 'is_update_way_R2L', false, 'init_lambda', 1.0);
%> >> F = hierarchical(M, params, 'backend', 2016)
%> Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 1/3<br/>
%> Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 2/3<br/>
%> Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 3/3<br/>
%>
%>  F =
%>
%> Faust size 500x32, density 0.189063, nnz_sum 3025, 4 factor(s):
%> - FACTOR 0 (double) SPARSE, size 500x32, density 0.15625, nnz 2500
%> - FACTOR 1 (double) SPARSE, size 32x32, density 0.09375, nnz 96
%> - FACTOR 2 (double) SPARSE, size 32x32, density 0.09375, nnz 96
%> - FACTOR 3 (double) SPARSE, size 32x32, density 0.325195, nnz 333
%>
%> >>
%>  @endcode
%>
%>  @b Example 2: Simplified Parameters for Hadamard Factorization
%>@code
%> >> import matfaust.*
%> >> import matfaust.fact.hierarchical
%> >> % generate a Hadamard Faust of size 32x32
%> >> FH = wht(32);
%> >> H = full(FH); % the full matrix version
%> >> % factorize it
%> >> FH2 = hierarchical(H, 'hadamard')
%> Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 1/4
%> Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 2/4
%> Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 3/4
%> Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 4/4
%>
%> FH2 =
%>
%> Faust size 32x32, density 0.3125, nnz_sum 320, 5 factor(s):
%> - FACTOR 0 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 1 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 2 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 3 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 4 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%>
%> >> % test the relative error
%> >> norm(FH-FH2, 'fro')/norm(FH, 'fro') < 1e-15 % the result is about 1e-16, the factorization is accurate doctest: +ELLIPSIS
%>
%> ans =
%>
%>   ... logical ...
%>
%>    1
%>
%> >>
%>@endcode
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
%> >> MEG16 = hierarchical(MEG, {'rectmat', num_facts, k, s}) % too long for doctest, % doctest: +SKIP
%>
%> MEG16 =
%>
%> Faust size 204x8193, density 0.0631655, nnz_sum 105573, 9 factor(s):
%> - FACTOR 0 (double) SPARSE, size 204x204, density 0.293613, nnz 12219
%> - FACTOR 1 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
%> - FACTOR 2 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
%> - FACTOR 3 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
%> - FACTOR 4 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
%> - FACTOR 5 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
%> - FACTOR 6 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
%> - FACTOR 7 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
%> - FACTOR 8 (double) SPARSE, size 204x8193, density 0.0490196, nnz 81930
%>
%> >> % verify the constraint k == 10, on column 5
%> >> fac9 = factors(MEG16,9); % doctest: +SKIP
%> >> numel(nonzeros(fac9(:,5))) % doctest: +SKIP
%>
%> ans =
%>
%>      10
%>
%> >> % now verify the s constraint is respected on MEG16 factor 2
%> >> numel(nonzeros(factors(MEG16, 2)))/size(MEG16,1) % doctest: +SKIP
%>
%> ans =
%>
%>      8
%>
%> >>
%> @endcode
%>
%> <br/>
%> @b Example 4: Simplified Parameters for Discrete Fourier Transform Factorization
%> @code
%> >> import matfaust.*
%> >> import matfaust.fact.hierarchical
%> >> % generate a DFT Faust of size 32x32
%> >> FDFT = dft(32);
%> >> DFT = full(FDFT); % the full matrix version
%> >> % factorize it
%> >> FDFT2 = hierarchical(DFT, 'dft');
%>  Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 1/5
%>  Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 2/5
%>  Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 3/5
%>  Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 4/5
%>  Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 5/5
%>
%> >> % test the relative error
%> >> norm(FDFT-FDFT2, 'fro')/norm(FDFT, 'fro') < 1e-5 % the result is about 1e-6, the factorization is accurate doctest: +ELLIPSIS
%>
%> ans =
%>
%>   ... logical ...
%>
%>    1
%>
%> >> FDFT
%>
%> FDFT =
%>
%> Faust size 32x32, density 0.34375, nnz_sum 352, 6 factor(s):
%> - FACTOR 0 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 1 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 2 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 3 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 4 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 5 (complex) SPARSE, size 32x32, density 0.03125, nnz 32
%> >> FDFT2
%>
%> FDFT2 =
%>
%> Faust size 32x32, density 1.3125, nnz_sum 1344, 6 factor(s):
%> - FACTOR 0 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 1 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 2 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 3 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 4 (complex) SPARSE, size 32x32, density 0.0625, nnz 64
%> - FACTOR 5 (complex) SPARSE, size 32x32, density 1, nnz 1024
%> @endcode
%>
%> @b Example 5: Simplified Parameters for Hadamard Factorization without residual
%> constraints</b>
%>
%> This factorization parameterization is the same as the one shown in 2.
%> except that there is no constraints at all on residual factors. See
%> matfaust.factparams.ParamsHierarchicalNoResCons and
%> matfaust.factparams.ParamsHierarchicalWHTNoResCons for more details.
%> @code
%> >> import matfaust.wht
%> >> import matfaust.fact.hierarchical
%> >> % generate a Hadamard Faust of size 32x32
%> >> FH = wht(32);
%> >> H = full(FH); % the full matrix version
%> >> % factorize it
%> >> FH2 = hierarchical(H, 'hadamard_simple', 'backend', 2020);
%> Faust::hierarchical: 1/4<br/>
%> Faust::hierarchical: 2/4<br/>
%> Faust::hierarchical: 3/4<br/>
%> Faust::hierarchical: 4/4<br/>
%> >> % test the relative error
%> >> norm(H-full(FH2), 'fro')/norm(H, 'fro')
%>
%> ans =
%>
%>      0
%>
%> >> FH2
%>
%> FH2 =
%>
%>      Faust size 32x32, density 0.3125, nnz_sum 320, 5 factor(s):
%>      - FACTOR 0 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%>      - FACTOR 1 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%>      - FACTOR 2 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%>      - FACTOR 3 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%>      - FACTOR 4 (double) SPARSE, size 32x32, density 0.0625, nnz 64
%>
%> @endcode
%>
%> @b Example 6: Simplified Parameters for a Rectangular Matrix Factorization (the BSL demo MEG matrix) without residual constraints</b>
%>
%> The factorization parameterization shown here is the same as in 3.
%> except that there is no constraint at all on residual factors. See
%> matfaust.factparams.ParamsHierarchicalNoResCons and
%> matfaust.factparams.ParamsHierarchicalRectMatNoResCons for more details.
%> In the example below the MEG matrix is factorized according to the
%> parameterization shown in 3. (aka "MEG") and on the other hand with
%> the parameterization of interest here (aka "MEG_SIMPLE", with no
%> residual constraints), the approximate accuracy is quite the same so we
%> can conclude that on this case (as in 5.) removing residual constraints can
%> not only simplify the parameterization of hierarchical PALM4MSA but also
%> be as efficient.
%>
%> @code
%> >> import matfaust.fact.hierarchical
%> >> MEG = load('matrix_MEG.mat');
%> >> MEG = MEG.matrix.';
%> >> F1 = hierarchical(MEG, {'MEG', 5, 10, 8}, 'backend', 2020)
%> Faust::hierarchical: 1/4
%> Faust::hierarchical: 2/4
%> Faust::hierarchical: 3/4
%> Faust::hierarchical: 4/4
%>
%> F1 =
%>
%> Faust size 204x8193, density 0.0697972, nnz_sum 116657, 5 factor(s):
%> - FACTOR 0 (double) DENSE, size 204x204, density 0.716816, nnz 29831
%> - FACTOR 1 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
%> - FACTOR 2 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
%> - FACTOR 3 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
%> - FACTOR 4 (double) SPARSE, size 204x8193, density 0.0490196, nnz 81930
%>
%> >> F2 = hierarchical(MEG, {'MEG_SIMPLE', 5, 10, 8}, 'backend', 2020)
%> Faust::hierarchical: 1/4
%> Faust::hierarchical: 2/4
%> Faust::hierarchical: 3/4
%> Faust::hierarchical: 4/4
%>
%> F2 =
%>
%> Faust size 204x8193, density 0.0697972, nnz_sum 116657, 5 factor(s):
%> - FACTOR 0 (double) DENSE, size 204x204, density 0.716816, nnz 29831
%> - FACTOR 1 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
%> - FACTOR 2 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
%> - FACTOR 3 (double) SPARSE, size 204x204, density 0.0392157, nnz 1632
%> - FACTOR 4 (double) SPARSE, size 204x8193, density 0.0490196, nnz 81930
%>
%> >> % compare the errors:
%> >> norm(MEG-full(F2), 'fro') / norm(MEG, 'fro')
%>
%> ans =
%>
%>     0.1303
%>
%> >> norm(MEG-full(F1), 'fro') / norm(MEG, 'fro')
%>
%> ans =
%>
%>     0.1260
%> >>
%> @endcode
%>
%> <p> @b See @b also matfaust.faust_fact, factparams.ParamsHierarchical, factparams.ParamsHierarchicalWHT, factparams.ParamsHierarchicalRectMat,factparams.ParamsHierarchicalDFT
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
	mex_params = p.to_mex_struct(M);
	backend = 2016;
	nargin = length(varargin);
	gpu = false;
	is_float = strcmp(class(M), 'single');
	if(is_float)
		dtype = 'float';
	else
		dtype = 'double'; % also for complex double
	end
	if(nargin > 0)
		for i=1:nargin
			switch(varargin{i})
				case 'backend'
					if(nargin < i+1)
						error('keyword argument ''backend'' must be followed by 2016 or 2020')
					else
						backend = varargin{i+1};
					end
				case 'gpu'
					if(nargin == i || ~ islogical(varargin{i+1}))
						error('gpu keyword argument is not followed by a logical')
					else
						gpu = varargin{i+1};
					end
			end
		end
		if(~ (isscalar(backend) && floor(backend) == backend) || backend ~= 2016 && backend ~= 2020)
			backend
			error('backend must be a int equal to 2016 or 2020')
		end
		if(backend ~= 2020 && gpu == true)
			error('GPU implementation is only available for 2020 backend.')
		end
	end
	if(backend == 2016)
		if(isreal(M))
			if(is_float)
				[lambda, core_obj] = mexHierarchical_factRealFloat(M, mex_params);
			else
				[lambda, core_obj] = mexHierarchical_factReal(M, mex_params);
			end
		else
			[lambda, core_obj] = mexHierarchical_factCplx(M, mex_params);
		end
		F = Faust(core_obj, isreal(M), 'cpu', dtype, true); % 4th arg to copy factors
	elseif(backend == 2020)
		if(isreal(M))
			if(gpu)
				if(is_float)
					[lambda, core_obj] = mexHierarchical2020_gpu2RealFloat(M, mex_params);
				else
					[lambda, core_obj] = mexHierarchical2020_gpu2Real(M, mex_params);
				end
			else
				if(is_float)
					[lambda, core_obj] = mexHierarchical2020RealFloat(M, mex_params);
				else
					[lambda, core_obj] = mexHierarchical2020Real(M, mex_params);
				end
			end
		else
			if(gpu)
				[lambda, core_obj] = mexHierarchical2020_gpu2Cplx(M, mex_params);
			else
				[lambda, core_obj] = mexHierarchical2020Cplx(M, mex_params);
			end
		end
		F = Faust(core_obj, isreal(M), 'cpu', dtype, true); % 4th arg to copy factors
	end
	varargout = {F, lambda, p};
end
