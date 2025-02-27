classdef FaustFactoryTest < matlab.unittest.TestCase
	properties
		faust_paths
	end

	methods(TestClassSetup)
		function initializeRandSeed(this)
			rng('shuffle')
		end
	end

	methods (TestMethodSetup)
		function addFaustToPath(this)
			addpath(this.faust_paths{:})
			if exist('set_path') > 0
				set_path
			end
			import matfaust.Faust
		end
	end

	methods (Test)
		function test_palm4msa(this)
			disp('Test matfaust.fact.palm4msa()')
			import matfaust.*
			import matfaust.factparams.*
			num_facts = 2;
			is_update_way_R2L = false;
			init_lambda = 1.0;
			init_facts = cell(2,1);
			init_facts{1} = zeros(500,32);
			init_facts{2} = eye(32);
			%M = rand(500, 32)
			load([this.faust_paths{1},'../../../misc/data/mat/config_compared_palm2.mat']);
			% data matrix is loaded from file
			M = data;
			cons = cell(2,1);
			cons{1} = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5);
			cons{2} = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32, 32, 1.0);
			stop_crit = StoppingCriterion(200);
			params = ParamsPalm4MSA(cons, stop_crit, 'is_update_way_R2L', is_update_way_R2L, 'init_lambda', init_lambda, 'step_size', ParamsFact.DEFAULT_STEP_SIZE,...
				'constant_step_size', ParamsFact.DEFAULT_CONSTANT_STEP_SIZE);
			F = matfaust.fact.palm4msa(M, params);
			this.verifyEqual(size(F), size(M))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-M;
			disp('err: ')
			norm(E,'fro')/norm(M,'fro')
			% matrix to factorize and reference relative error come from
			% misc/test/src/C++/test_palm4MSA.cpp
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'), 0.2710, 'AbsTol', 0.0001)
		end

		function test_palm4msa2020(this)
			disp('Test matfaust.fact.palm4msa2020()')
			import matfaust.*
			import matfaust.factparams.*
			num_facts = 2;
			is_update_way_R2L = false;
			init_lambda = 1.0;
			init_facts = cell(2,1);
			init_facts{1} = zeros(500,32);
			init_facts{2} = eye(32);
			%M = rand(500, 32)
			load([this.faust_paths{1},'../../../misc/data/mat/config_compared_palm2.mat']);
			% data matrix is loaded from file
			M = data;
			cons = cell(2,1);
			cons{1} = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5);
			cons{2} = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32, 32, 1.0);
			stop_crit = StoppingCriterion(200);
			params = ParamsPalm4MSA(cons, stop_crit, 'is_update_way_R2L', is_update_way_R2L, 'init_lambda', init_lambda, 'step_size', ParamsFact.DEFAULT_STEP_SIZE,...
				'constant_step_size', ParamsFact.DEFAULT_CONSTANT_STEP_SIZE);
			F = matfaust.fact.palm4msa(M, params, 'backend', 2020);
			this.verifyEqual(size(F), size(M))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-M;
			disp('err: ')
			norm(E,'fro')/norm(M,'fro')
			% matrix to factorize and reference relative error come from
			% misc/test/src/C++/test_palm4MSA.cpp
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'), 0.2710, 'AbsTol', 0.0001)
		end

		function test_palm4msaCplxDeftInitFacts(this)
			disp('Test matfaust.fact.palm4msaCplxDeftInitFacts()')
			import matfaust.factparams.*
			num_facts = 2;
			is_update_way_R2L = false;
			init_lambda = 1.0;
			%M = rand(500, 32)
			load([this.faust_paths{1},'../../../misc/data/mat/config_compared_palm2.mat']);
			% data matrix is loaded from file
			M = data+j*data;
			cons = cell(2,1);
			cons{1} = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5);
			cons{2} = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32, 32, 1.0);
			stop_crit = StoppingCriterion(200);
			params = ParamsPalm4MSA(cons, stop_crit, 'is_update_way_R2L', is_update_way_R2L, 'init_lambda', init_lambda);
			F = matfaust.fact.palm4msa(M, params);
			this.verifyEqual(size(F), size(M))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-M;
			disp('err: ')
			norm(E,'fro')/norm(M,'fro')
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'), 0.272814, 'AbsTol', 0.0001)
		end

		function test_palm4msaCplx(this)
			disp('Test matfaust.fact.palm4msaCplx()')
			import matfaust.factparams.*
			num_facts = 2;
			is_update_way_R2L = false;
			init_lambda = 1.0;
			init_facts = cell(2,1);
			init_facts{1} = complex(zeros(500,32)); % type consistency with M is verified by palm4msa
			init_facts{2} = complex(eye(32));
			%M = rand(500, 32)
			load([this.faust_paths{1},'../../../misc/data/mat/config_compared_palm2.mat']);
			% data matrix is loaded from file
			M = data+j*data;
			cons = cell(2,1);
			cons{1} = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5);
			cons{2} = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32, 32, 1.0);
			stop_crit = StoppingCriterion(200);
			params = ParamsPalm4MSA(cons, stop_crit, 'is_update_way_R2L', is_update_way_R2L, 'init_lambda', init_lambda, 'init_facts', init_facts);
			F = matfaust.fact.palm4msa(M, params);
			this.verifyEqual(size(F), size(M))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-M;
			disp('err: ')
			norm(E,'fro')/norm(M,'fro')
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'), 0.272814, 'AbsTol', 0.0001)
		end

		function test_palm4msa2020Cplx(this)
			disp('Test matfaust.fact.palm4msa2020Cplx()')
			import matfaust.factparams.*
			num_facts = 2;
			is_update_way_R2L = false;
			init_lambda = 1.0;
			init_facts = cell(2,1);
			init_facts{1} = complex(zeros(500,32)); % type consistency with M is verified by palm4msa
			init_facts{2} = complex(eye(32));
			%M = rand(500, 32)
			load([this.faust_paths{1},'../../../misc/data/mat/config_compared_palm2.mat']);
			% data matrix is loaded from file
			M = data+j*data;
			cons = cell(2,1);
			cons{1} = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5);
			cons{2} = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32, 32, 1.0);
			stop_crit = StoppingCriterion(200);
			params = ParamsPalm4MSA(cons, stop_crit, 'is_update_way_R2L', is_update_way_R2L, 'init_lambda', init_lambda, 'init_facts', init_facts);
			F = matfaust.fact.palm4msa(M, params, 'backend', 2020);
			this.verifyEqual(size(F), size(M))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-M;
			disp('err: ')
			norm(E,'fro')/norm(M,'fro')
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'), 0.272814, 'AbsTol', 0.0001)
		end

		function test_palm4msa2020CplxDeftInitFacts(this)
			disp('Test matfaust.fact.palm4msa2020CplxDeftInitFacts()')
			import matfaust.factparams.*
			num_facts = 2;
			is_update_way_R2L = false;
			init_lambda = 1.0;
			load([this.faust_paths{1},'../../../misc/data/mat/config_compared_palm2.mat']);
			% data matrix is loaded from file
			M = data+j*data;
			cons = cell(2,1);
			cons{1} = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5);
			cons{2} = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32, 32, 1.0);
			stop_crit = StoppingCriterion(200);
			params = ParamsPalm4MSA(cons, stop_crit, 'is_update_way_R2L', is_update_way_R2L, 'init_lambda', init_lambda);
			F = matfaust.fact.palm4msa(M, params, 'backend', 2020)
			this.verifyEqual(size(F), size(M))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-M;
			disp('err: ')
			norm(E,'fro')/norm(M,'fro')
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'), 0.272814, 'AbsTol', 0.0001)
		end

		function test_hierarchical(this)
			disp('Test matfaust.fact.hierarchical()')
			import matfaust.factparams.*
			num_facts = 4;
			is_update_way_R2L = false;
			init_lambda = 1.0;
			%init_facts = cell(num_facts,1);
			%init_facts{1} = zeros(500,32);
			%for i=2:num_facts
			%init_facts{i} = zeros(32);
			%end
			%M = rand(500, 32)
			load([this.faust_paths{1},'../../../misc/data/mat/matrix_hierarchical_fact.mat'])
			% matrix is loaded from file
			M = matrix;
			fact_cons = cell(3,1);
			res_cons = cell(3, 1);
			fact_cons{1} = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5);
			fact_cons{2} = ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 96);
			fact_cons{3} = ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 96);
			res_cons{1} = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32, 32, 1);
			res_cons{2} =  ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 666);
			res_cons{3} =  ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 333);
			stop_crit = StoppingCriterion(200);
			stop_crit2 = StoppingCriterion(200);
			params = ParamsHierarchical(fact_cons, res_cons, stop_crit, stop_crit2);
			F = matfaust.fact.hierarchical(M, params);
			this.verifyEqual(size(F), size(M))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-M;
			disp('err: ')
			norm(E,'fro')/norm(M,'fro')
			% matrix to factorize and reference relative error come from
			% misc/test/src/C++/hierarchicalFactorization.cpp
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'), 0.26844, 'AbsTol', 0.001)
		end

		function test_hierarchical2020(this)
			disp('Test matfaust.fact.hierarchical()')
			import matfaust.factparams.*
			num_facts = 4;
			is_update_way_R2L = false;
			init_lambda = 1.0;
			%init_facts = cell(num_facts,1);
			%init_facts{1} = zeros(500,32);
			%for i=2:num_facts
			%init_facts{i} = zeros(32);
			%end
			%M = rand(500, 32)
			load([this.faust_paths{1},'../../../misc/data/mat/matrix_hierarchical_fact.mat'])
			% matrix is loaded from file
			M = matrix;
			fact_cons = cell(3,1);
			res_cons = cell(3, 1);
			fact_cons{1} = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5);
			fact_cons{2} = ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 96);
			fact_cons{3} = ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 96);
			res_cons{1} = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32, 32, 1);
			res_cons{2} =  ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 666);
			res_cons{3} =  ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 333);
			stop_crit = StoppingCriterion(200);
			stop_crit2 = StoppingCriterion(200);
			params = ParamsHierarchical(fact_cons, res_cons, stop_crit, stop_crit2);
			F = matfaust.fact.hierarchical(M, params, 'backend', 2020);
			this.verifyEqual(size(F), size(M))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-M;
			disp('err: ')
			norm(E,'fro')/norm(M,'fro')
			% matrix to factorize and reference relative error come from
			% misc/test/src/C++/hierarchicalFactorization.cpp
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'), 0.26844, 'AbsTol', 0.001)
		end

		function test_hierarchicalCplx(this)
			disp('Test matfaust.fact.hierarchicalCplx()')
			import matfaust.factparams.*
			num_facts = 4;
			is_update_way_R2L = false;
			init_lambda = 1.0;
			%init_facts = cell(num_facts,1);
			%init_facts{1} = zeros(500,32);
			%for i=2:num_facts
			%init_facts{i} = zeros(32);
			%end
			%M = rand(500, 32)
			load([this.faust_paths{1},'../../../misc/data/mat/matrix_hierarchical_fact.mat'])
			% matrix is loaded from file
			M = matrix+j*matrix;
			fact_cons = cell(3,1);
			res_cons = cell(3, 1);
			fact_cons{1} = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5);
			fact_cons{2} = ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 96);
			fact_cons{3} = ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 96);
			res_cons{1} = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32, 32, 1);
			res_cons{2} =  ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 666);
			res_cons{3} =  ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 333);
			stop_crit = StoppingCriterion(200);
			stop_crit2 = StoppingCriterion(200);
			params = ParamsHierarchical(fact_cons, res_cons, stop_crit, stop_crit2,...
				'init_lambda', init_lambda, 'is_update_way_R2L', is_update_way_R2L);
			F = matfaust.fact.hierarchical(M, params);
			this.verifyEqual(size(F), size(M))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-M;
			disp('err: ')
			norm(E,'fro')/norm(M,'fro')
			% matrix to factorize and reference relative error come from
			% misc/test/src/C++/hierarchicalFactorization.cpp
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'),0.2727, 'AbsTol', 0.0001)
		end

		function test_hierarchical2020Cplx(this)
			disp('Test matfaust.fact.hierarchical2020Cplx()')
			import matfaust.factparams.*
			num_facts = 4;
			is_update_way_R2L = false;
			init_lambda = 1.0;
			%init_facts = cell(num_facts,1);
			%init_facts{1} = zeros(500,32);
			%for i=2:num_facts
			%init_facts{i} = zeros(32);
			%end
			%M = rand(500, 32)
			load([this.faust_paths{1},'../../../misc/data/mat/matrix_hierarchical_fact.mat'])
			% matrix is loaded from file
			M = matrix+j*matrix;
			fact_cons = cell(3,1);
			res_cons = cell(3, 1);
			fact_cons{1} = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5);
			fact_cons{2} = ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 96);
			fact_cons{3} = ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 96);
			res_cons{1} = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32, 32, 1);
			res_cons{2} =  ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 666);
			res_cons{3} =  ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 333);
			stop_crit = StoppingCriterion(200);
			stop_crit2 = StoppingCriterion(200);
			params = ParamsHierarchical(fact_cons, res_cons, stop_crit, stop_crit2,...
				'init_lambda', init_lambda, 'is_update_way_R2L', is_update_way_R2L);
			F = matfaust.fact.hierarchical(M, params, 'backend', 2020);
			this.verifyEqual(size(F), size(M))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-M;
			disp('err: ')
			norm(E,'fro')/norm(M,'fro')
			% matrix to factorize and reference relative error come from
			% misc/test/src/C++/hierarchicalFactorization.cpp
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'),0.2727, 'AbsTol', 0.0001)
		end

		function test_fgft_givens(this)
			disp('Test matfaust.fact.fgft_givens()')
			load([this.faust_paths{1} '../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat'])
			% Lap and J available
			[F,D] = matfaust.fact.eigtj(Lap, 'nGivens', J, 'enable_large_Faust', true, 'nGivens_per_fac', 1);%, 0, 'verbosity', 1);
			this.verifyEqual(size(F), size(Lap))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = F*full(D)*F'-Lap;
			err = norm(E,'fro')/norm(Lap,'fro')
			% the error reference is from the C++ test,
			% misc/test/src/C++/EigTJ.cpp.in
			this.verifyEqual(err, 0.0326529, 'AbsTol', 0.00001)
			% verify it works the same using the eigtj() alias function
			[F2,D2] = matfaust.fact.eigtj(Lap, 'nGivens', J, 'enable_large_Faust', true, 'nGivens_per_fac', 1);%, 0, 'verbosity', 2);
			this.verifyEqual(full(F2),full(F))
			this.verifyEqual(D,D2)
		end

		function test_fgft_givens_sparse(this)
			disp('Test matfaust.fact.fgft_givens()')
			load([this.faust_paths{1} '../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat'])
			% Lap and J available
			[F,D] = matfaust.fact.eigtj(sparse(Lap), 'nGivens', J, 'enable_large_Faust', true, 'nGivens_per_fac', 1);%, 0, 'verbosity', 1);
			this.verifyEqual(size(F), size(Lap))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = F*full(D)*F'-Lap;
			err = norm(E,'fro')/norm(Lap,'fro');
			% the error reference is from the C++ test,
			% misc/test/src/C++/EigTJ.cpp.in
			this.verifyEqual(err, 0.0326529, 'AbsTol', 0.00001)
			% verify it works the same using the eigtj() alias function
			[F2,D2] = matfaust.fact.eigtj(Lap, 'nGivens', J, 'nGivens_per_fac', 1, 'enable_large_Faust', true);%, 0, 'verbosity', 2);
			this.verifyEqual(full(F2),full(F))
			this.verifyEqual(D,D2)
		end

		function test_fgft_givens_parallel(this)
			disp('Test matfaust.fact.fgft_givens() -- parallel')
			load([this.faust_paths{1} '../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat'])
			% Lap and J available
			t = size(Lap,1)/2;
			[F,D] = matfaust.fact.eigtj(Lap, 'nGivens', J, 'nGivens_per_fac', t, 'enable_large_Faust', true); %, 'verbosity', 2);
			this.verifyEqual(size(F), size(Lap))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = F*full(D)*F'-Lap;
			err = norm(E,'fro')/norm(Lap,'fro')
			% the error reference is from the C++ test,
			% misc/test/src/C++/EigTJParallel.cpp.in
			this.verifyEqual(err,0.0398154, 'AbsTol', 0.00001)
			% verify it works the same using the eigtj() alias function
			[F2,D2] = matfaust.fact.eigtj(Lap, 'nGivens', J, 'nGivens_per_fac', t, 'enable_large_Faust', true); %, 'verbosity', 2);
			this.verifyEqual(full(F2),full(F))
			this.verifyEqual(D,D2)
		end

		function test_fgft_givens_parallel_sparse(this)
			disp('Test matfaust.fact.fgft_givens_sparse() -- parallel')
			load([this.faust_paths{1} '../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat'])
			% Lap and J available
			t = size(Lap,1)/2;
			[F,D] = matfaust.fact.eigtj(sparse(Lap), 'nGivens', J, 'nGivens_per_fac', t, 'enable_large_Faust', true); %, 'verbosity', 2);
			this.verifyEqual(size(F), size(Lap))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = F*full(D)*F'-Lap;
			err = norm(E,'fro')/norm(Lap,'fro')
			% the error reference is from the C++ test,
			% misc/test/src/C++/EigTJParallel.cpp.in
			this.verifyEqual(err, 0.0398154, 'AbsTol', 0.00001)
			% verify it works the same using the eigtj() alias function
			[F2,D2] = matfaust.fact.eigtj(Lap, 'nGivens', J, 'nGivens_per_fac', t, 'enable_large_Faust', true); %, 'verbosity', 2);
			this.verifyEqual(full(F2),full(F))
			this.verifyEqual(D,D2)
		end

		function test_fgft_palm(this)
			disp('Test matfaust.fact.fgft_palm()')
			import matfaust.factparams.*
			num_facts = 4;
			is_update_way_R2L = false;
			init_lambda = 1.0;

			load([this.faust_paths{1},'../../../misc/data/mat/HierarchicalFactFFT_test_U_L_params.mat'])
			% U, Lap, init_D, params are loaded from file

			fact_cons = cell(3,1);
			res_cons = cell(3, 1);
			fact_cons{1} = ConstraintInt(ConstraintName(ConstraintName.SP), 128, 128, 12288);
			fact_cons{2} = ConstraintInt(ConstraintName(ConstraintName.SP), 128, 128, 6144);
			fact_cons{3} = ConstraintInt(ConstraintName(ConstraintName.SP), 128, 128, 3072);
			res_cons{1} = ConstraintInt(ConstraintName(ConstraintName.SP), 128, 128, 384);
			res_cons{2} =  ConstraintInt(ConstraintName(ConstraintName.SP), 128, 128, 384);
			res_cons{3} =  ConstraintInt(ConstraintName(ConstraintName.SP), 128, 128, 384);
			stop_crit = StoppingCriterion(params.niter1);
			stop_crit2 = StoppingCriterion(params.niter2);
			params.fact_side = 0; % forced
			params.verbose = false; % forced
			params.init_lambda = 128;
			params.step_size = 1e-16
			params.grad_calc_opt_mode = 2
			params = ParamsHierarchical(fact_cons, res_cons, stop_crit, stop_crit2, 'is_fact_side_left', params.fact_side == 1, 'is_update_way_R2L', params.update_way == 1, 'init_lambda', params.init_lambda, 'step_size', params.step_size, 'constant_step_size', false, 'is_verbose', params.verbose);
			diag_init_D = diag(init_D);
			[F,D,lambda] = matfaust.fact.fgft_palm(U, Lap, params, diag_init_D);
			this.verifyEqual(size(F), size(U))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-U;
			err = norm(E,'fro')/norm(U,'fro')
			% matrix to factorize and reference relative error come from
			% misc/test/src/C++/hierarchicalFactorizationFFT.cpp
			this.verifyEqual(err, 0.0555, 'AbsTol', 0.0001)
		end

		function test_svdtj(this)
			disp('Test matfaust.fact.svdtj')
			import matfaust.fact.svdtj
			M = rand(16, 32);
			[U5, S5, V5] = svdtj(M, 'tol', 1e-3, 'enable_large_Faust', false);
			this.verifyLessThanOrEqual(norm(U5 * S5 * V5' - M) / norm(M), 1e-3)
		end

		function testHadamard(this)
			disp('Test matfaust.wht()')
			import matfaust.*
			log2n = randi([1,6]);
			n = 2^log2n;
			H = wht(n, 'normed', false);
			fH = full(H);
			this.verifyEqual(nnz(fH), numel(fH));
			i = 1;
			while(i<n-1)
				j = i + 1;
				while(j<n)
					this.verifyEqual(fH(i,:)*fH(j,:)',0);
					j = j + 1;
				end
				i = i + 1;
			end
			this.assertEqual(full(wht(n)), full(normalize(H)), 'AbsTol', 10^-7)

			this.assertEqual(full(wht(n, 'normed', false)), full(H), 'AbsTol', 10^-7)
		end

		function testFourier(this)
			disp('Test matfaust.dft()')
			import matfaust.dft
			log2n = randi([1,10]);
			n = 2^log2n;
			F = dft(n, 'normed', false);
			fF = full(F);
			fftI = fft(eye(n));
			% this.verifyEqual(nnz(fH), numel(fH));
			this.verifyEqual(norm(fF-fftI), 0, 'AbsTol', 10^-11);
			this.assertEqual(full(dft(n)), full(normalize(F)), 'AbsTol', 10^-7)

			this.assertEqual(full(dft(n, 'normed', false)), full(F), 'AbsTol', 10^-7)
		end

		function testFourierDiagOpt(this)
			disp('Test matfaust.dftDiagOpt()')
			import matfaust.dft
			log2n = randi([1,10]);
			n = 2^log2n;
			normed_vals = {true, false};
			for i=1:2
				normed = normed_vals{i};
				F = dft(n, 'normed', normed, 'diag_opt', true);
				fF = full(F);
				frefF = full(dft(n, 'normed', normed, 'diag_opt', false));
				this.verifyEqual(norm(fF-frefF), 0, 'AbsTol', 10^-12);
			end
		end

		function testRandButterfly(this)
			disp('Test matfaust.rand_butterfly')
			import matfaust.wht
			import matfaust.rand_butterfly
			H = full(wht(32));
			classes = {'single', 'double'};
			fields = {'real', 'complex'};
			for i=1:2
				for j=1:2
					field = fields{j};
					if strcmp(field, 'real')
						class = classes{i};
					else % field == complex
						class = classes{2};
					end
					F = rand_butterfly(32, 'class', class, 'field', field);
					this.assertNotEqual(full(F), H);
					[ref_I, ref_J, ~] = find(H);
					[I, J, ~] = find(full(F));
					this.assertEqual(ref_I, I);
					this.assertEqual(ref_J, J);
				end
			end
		end

		function testCircAntiCirc(this)
			disp('Test matfaust.circ/anticirc')
			c = [1 2 3 4];
			C = [[1     4     3     2];
			[2     1     4     3]
			[3     2     1     4]
			[4     3     2     1]];
			this.assertEqual(C, real(full(matfaust.circ(c))))
			this.assertEqual(C, real(full(matfaust.circ(c, 'diag_opt', true))))
			this.assertEqual(matfaust.circ(c) * reshape(c, numel(c), 1), matfaust.circ(c, 'diag_opt', true) * reshape(c, numel(c), 1))
			A = [[2     3     4     1]
			[3     4     1     2]
			[4     1     2     3]
			[1     2     3     4]];
			this.assertEqual(A, real(full(matfaust.anticirc(c))))
			this.assertEqual(A, real(full(matfaust.anticirc(c, 'diag_opt', true))))
			this.assertEqual(matfaust.anticirc(c) * reshape(c, numel(c), 1), matfaust.anticirc(c, 'diag_opt', true) * reshape(c, numel(c), 1))
			% not a power of two
			c = [1 2 3 4 5];
			C = [[1     5     4     3     2]
			[2     1     5     4     3]
			[3     2     1     5     4]
			[4     3     2     1     5]
			[5     4     3     2     1]];
			this.assertEqual(C, real(full(matfaust.circ(c))), 'AbsTol', 1e-8)
			this.assertEqual(C, real(full(matfaust.circ(c, 'diag_opt', true))), 'AbsTol', 1e-8)
			this.assertEqual(matfaust.circ(c) * reshape(c, numel(c), 1), matfaust.circ(c, 'diag_opt', true) * reshape(c, numel(c), 1))
			A = [[2     3     4     5    1]
			[3     4     5     1    2]
			[4     5     1     2    3]
			[5     1     2     3    4]
			[1     2     3     4    5]];
			this.assertEqual(A, real(full(matfaust.anticirc(c))), 'AbsTol', 1e-8)
			this.assertEqual(A, real(full(matfaust.anticirc(c, 'diag_opt', true))), 'AbsTol', 1e-8)
			this.assertEqual(matfaust.anticirc(c) * reshape(c, numel(c), 1), matfaust.anticirc(c, 'diag_opt', true) * reshape(c, numel(c), 1))
		end

		function testToeplitz(this)
			disp('Test matfaust.toeplitz')
			c = [1 2 3 4];
			r = [1 7 8];
			T = [[1., 7., 8.]
			[2., 1., 7.]
			[3., 2., 1.]
			[4., 3., 2.]];
			this.verifyEqual(T, real(full(matfaust.toeplitz(c, r))), 'AbsTol', 1e-6)
			this.verifyEqual(T, real(full(matfaust.toeplitz(c, r, 'diag_opt', true))), 'AbsTol', 1e-6)
			this.verifyEqual(toeplitz(c, r), real(full(matfaust.toeplitz(c, r))), 'AbsTol', 1e-6)
			this.verifyEqual(toeplitz(c, r), real(full(matfaust.toeplitz(c, r, 'diag_opt', true))), 'AbsTol', 1e-6)
			this.verifyEqual(toeplitz(c), real(full(matfaust.toeplitz(c))), 'AbsTol', 1e-6)
			this.verifyEqual(toeplitz(c), real(full(matfaust.toeplitz(c, 'diag_opt', true))), 'AbsTol', 1e-6)
			this.verifyEqual(toeplitz(c) * reshape(c, numel(c), 1) , matfaust.toeplitz(c, 'diag_opt', true) * reshape(c, numel(c), 1) , 'AbsTol', 1e-6)
			this.verifyEqual(toeplitz(c) * reshape(c, numel(c), 1) , matfaust.toeplitz(c, 'diag_opt', false) * reshape(c, numel(c), 1) , 'AbsTol', 1e-6)
		end

		function test_sp(this)
			disp('Test matfaust.proj.sp')
			import matfaust.proj.sp;
			min_n = 5;
			min_m = 5;
			m = randi([min_m 128], 1);
			n = randi([min_n 128], 1);
			M = rand(m,n);
			k = randi([1 m*n], 1);
			p = sp([m n], k, 'normalized', true);
			Mp = p(M);
			for i=1:n
				this.verifyLessThanOrEqual(length(nonzeros(Mp(:,i))), k)
			end
			this.verifyEqual(norm(Mp, 'fro'), 1, 'AbsTol', 10^-3)
		end

		function test_splincol(this)
			disp('Test matfaust.proj.splincol')
			import matfaust.proj.splincol
			min_n = 5;
			min_m = 5;
			m = randi([min_m 128], 1);
			n = randi([min_n 128], 1);
			k = randi([1 m], 1);
			M = rand(m,n);
			p = splincol([m, n], k, 'normalized', true);
			Mp = p(M);
			% TODO: test it is union of splin and spcol
			this.verifyEqual(norm(Mp, 'fro'), 1, 'AbsTol', 10^-3)
		end

		function test_spcol(this)
			disp('Test matfaust.proj.spcol')
			import matfaust.proj.spcol
			min_n = 5;
			min_m = 5;
			m = randi([min_m 128], 1);
			n = randi([min_n 128], 1);
			M = rand(m,n);
			k = randi([1 m], 1);
			p = spcol([m, n], k, 'normalized', true);
			Mp = p(M);
			for i=1:n
				this.verifyLessThanOrEqual(length(nonzeros(Mp(:,i))), k)
			end
			this.verifyEqual(norm(Mp, 'fro'), 1, 'AbsTol', 10^-3)
		end

		function test_normcol(this)
			disp('Test matfaust.proj.normcol')
			import matfaust.proj.normcol
			min_n = 5;
			min_m = 5;
			m = randi([min_m 128], 1);
			n = randi([min_n 128], 1);
			M = rand(m,n);
			s = rand(1,1)*50;
			p = normcol([m, n], 's', s);
			Mp = p(M);
			for i=1:n
				this.verifyLessThanOrEqual(norm(Mp(:,i)), s+1e-6)
			end
		end

		function test_normlin(this)
			disp('Test matfaust.proj.normlin')
			import matfaust.proj.normlin
			min_n = 5;
			min_m = 5;
			m = randi([min_m 128], 1);
			n = randi([min_n 128], 1);
			M = rand(m,n);
			s = rand()*50;
			p = normlin([m, n], 's', s);
			Mp = p(M);
			for i=1:m
				this.verifyLessThanOrEqual(norm(Mp(i,:)), s+1e-6)
			end
		end

		function test_const(this)
			disp('Test matfaust.proj.const')
			import matfaust.proj.const
			min_n = 5;
			min_m = 5;
			m = randi([min_m 128], 1);
			n = randi([min_n 128], 1);
			M = rand(m,n);
			C = rand(m,n);
			s = rand()*50;
			p = const(C);
			Mp = p(M);
			this.verifyEqual(norm(C-Mp)/norm(C), 0)
		end

		function test_supp(this)
			disp('Test matfaust.proj.supp')
			import matfaust.proj.supp
			min_n = 5;
			min_m = 5;
			m = randi([min_m 128], 1);
			n = randi([min_n 128], 1);
			M = rand(m,n);
			S = rand(m,n) > .5;
			p = supp(S, 'normalized', false);
			Mp = p(M);
			this.assertTrue(all(all(Mp > 0 == S)))
			this.assertTrue(all(all(Mp(Mp > 0) == M(Mp > 0))))
		end

		function test_skperm(this)
			disp('Test matfaust.proj.skperm')
			import matfaust.proj.skperm
			M = [-0.04440802, -0.17569296, -0.02557815, -0.15559154;
			-0.0083095,  -3.38725936, -0.78484126, -0.4883618;
			-1.48942563, -1.71787215, -0.84000212, -3.71752454;
			-0.88957883, -0.19107863, -5.92900636, -6.51064175];
			p = skperm(size(M), 2, 'normalized', false);
			pM_ref = [-0.04440802,0.,-0.02557815,0.;
			-0.0083095,-3.38725936,0.,0.;
			0.,-1.71787215,0.,-3.71752454;
			0.,0.,-5.92900636,-6.51064175];
			pM = p(M);
			this.verifyEqual(norm(pM-pM_ref)/norm(pM_ref), 0)
		end

		function test_sptriu(this)
			disp('Test matfaust.proj.sptriu')
			import matfaust.proj.sptriu
			M = rand(5, 5);
			k = 2;
			p = sptriu(size(M), k);
			M_ = p(M);
			this.verifyEqual(numel(nonzeros(M_)), k)
			this.verifyEqual(norm(tril(M_, -1)), 0)
		end

		function test_sptril(this)
			disp('Test matfaust.proj.sptril')
			import matfaust.proj.sptril
			M = rand(5, 5);
			k = 2;
			p = sptril(size(M), k);
			M_ = p(M);
			this.verifyEqual(numel(nonzeros(M_)), k)
			this.verifyEqual(norm(triu(M_, 1)), 0)
		end

		function test_spsymm(this)
			disp('Test matfaust.proj.spsymm')
			import matfaust.proj.spsymm
			M = rand(5, 5);
			k = 2;
			p = spsymm(size(M), k);
			M_ = p(M);
			this.verifyEqual(M_.', M_)
			this.verifyTrue(numel(nonzeros(M_)) - k <= 1)
		end

		function test_butterfly(this)
			disp('Test matfaust.fact.butterfly')
			import matfaust.wht
			import matfaust.dft
			import matfaust.fact.butterfly
			D = full(dft(64));
			H = full(wht(64));
			types = {'left', 'right', 'bbtree'};
			for t=1:2
				if(t == 2)
					M = D;
				else
					M = H;
				end
				for i=1:3
					type = types{i};
					args = {M, 'type', type};
					if ~ isreal(M)
						args = [args, {'perm', 'bitrev'}];
					end
					F = butterfly(args{:});
					this.verifyEqual(norm(full(F)-M)/norm(M), 0, 'AbsTol', 1e-6);
				end
			end
		end

		function test_palm4msa_mhtp(this)
			disp('Test matfaust.fact.palm4msa_mhtp')
			import matfaust.fact.palm4msa_mhtp
			import matfaust.factparams.*
			import matfaust.proj.*
			M = rand(500,32);
			projs = { splin([500,32], 5), normcol([32,32], 1.0)};
			stop_crit = StoppingCriterion(200);
			param = ParamsPalm4MSA(projs, stop_crit);
			mhtp_param = MHTPParams('num_its', 60, 'palm4msa_period', 10);
			G = palm4msa_mhtp(M, param, mhtp_param);
			err = norm(full(G)-M)/norm(M)
			this.verifyLessThanOrEqual(err, 0.20)
		end

		function test_hierarchical_mhtp(this)
			disp('Test matfaust.fact.hierarchical_mhtp')
			import matfaust.fact.hierarchical_mhtp
			import matfaust.factparams.ParamsHierarchical
			import matfaust.factparams.StoppingCriterion
			import matfaust.factparams.MHTPParams
			import matfaust.proj.*
			M = rand(500,32);
			fact_projs = { splin([500,32], 5), sp([32,32], 96), sp([32, 32], 96)};
			res_projs = { normcol([32,32], 1), sp([32,32], 666), sp([32, 32], 333)};
			stop_crit1 = StoppingCriterion(200);
			stop_crit2 = StoppingCriterion(200);
			% 50 iterations of MHTP will run every 100 iterations of PALM4MSA (each time PALM4MSA is called by the hierarchical algorithm)
			mhtp_param = MHTPParams('num_its', 150, 'palm4msa_period', 100);
			param = ParamsHierarchical(fact_projs, res_projs, stop_crit1, stop_crit2);
			F = hierarchical_mhtp(M, param, mhtp_param);
			err = norm(full(F)-M)/norm(M)
			this.verifyLessThanOrEqual(err, 0.20)
		end

		function test_hierarchical_dft(this)
			disp('Test matfaust.fact.hierarchical')
			import matfaust.fact.hierarchical
			import matfaust.dft
			DFT = full(dft(32));
			F = hierarchical(DFT, 'dft', 'backend', 2020)
			err = norm(full(F)-DFT)/norm(DFT)
			this.verifyLessThanOrEqual(err, 1e-3)
			DFT = full(dft(32));
			F = hierarchical(DFT, 'dft', 'backend', 2016)
			err = norm(full(F)-DFT)/norm(DFT)
			this.verifyLessThanOrEqual(err, 1e-3)
		end

		function test_dct(this)
			disp('Test matfaust.dct')
			n = 128;
			% x = rand(n, 1);
			x = ones(n, 1);
			y_ref = zeros(n, 1);
			for k=1:n
				for i=1:n
					y_ref(k) = y_ref(k) + 2 * x(i) * cos(pi * (k-1) * (2 * (i-1) + 1) / 2 / n);
				end
			end
			DCT = matfaust.dct(n, 'normed', false);
			y_test = real(DCT*x);
			this.verifyEqual(y_ref, y_test, 'AbsTol', 1e-6);
			DCT_normed = matfaust.dct(n);
			DCT_normalized = normalize(DCT);
			this.verifyEqual(DCT_normed*x, DCT_normalized*x, 'AbsTol', 1e-6);
		end

		function test_dst(this)
			disp('Test matfaust.dst')
			n = 128;
			% x = rand(n, 1);
			x = ones(n, 1);
			y_ref = zeros(n, 1);
			for k=1:n
				for i=1:n
					y_ref(k) = y_ref(k) + x(i) * sin(pi / n * ((i-1) + 1 / 2) * k);
				end
			end
			y_ref = y_ref * 2;
			DST = matfaust.dst(n, 'normed', false);
			y_test = real(DST*x);
			this.verifyEqual(y_ref, y_test, 'AbsTol', 1e-6);
			DST_normed = matfaust.dst(n);
			DST_normalized = normalize(DST);
			this.verifyEqual(DST_normed*x, DST_normalized*x, 'AbsTol', 1e-6);
		end

	end

	methods
		function faust_test = FaustFactoryTest(varargin)
			faust_test.faust_paths = varargin
			run(faust_test)
		end
	end

end

