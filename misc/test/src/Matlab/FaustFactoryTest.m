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
			set_path
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
			F = matfaust.fact.palm4msa(M, params)
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
			F = matfaust.fact.palm4msa(M, params, 'backend', 2020)
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
			disp('Test matfaust.fact.palm4msaCplx()')
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
			F = matfaust.fact.palm4msa(M, params)
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
			F = matfaust.fact.palm4msa(M, params)
			this.verifyEqual(size(F), size(M))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-M;
			disp('err: ')
			norm(E,'fro')/norm(M,'fro')
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'), 0.272814, 'AbsTol', 0.0001)
		end

		function test_palm4msa2020Cplx(this)
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
			F = matfaust.fact.palm4msa(M, params, 'backend', 2020)
			this.verifyEqual(size(F), size(M))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-M;
			disp('err: ')
			norm(E,'fro')/norm(M,'fro')
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'), 0.272814, 'AbsTol', 0.0001)
		end

		function test_palm4msa2020CplxDeftInitFacts(this)
			disp('Test matfaust.fact.palm4msaCplx()')
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
			F = matfaust.fact.hierarchical(M, params)
			this.verifyEqual(size(F), size(M))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-M;
			disp('err: ')
			norm(E,'fro')/norm(M,'fro')
			% matrix to factorize and reference relative error come from
			% misc/test/src/C++/hierarchicalFactorization.cpp
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'), 0.26844, 'AbsTol', 0.00001)
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
			F = matfaust.fact.hierarchical(M, params, 'backend', 2020)
			this.verifyEqual(size(F), size(M))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-M;
			disp('err: ')
			norm(E,'fro')/norm(M,'fro')
			% matrix to factorize and reference relative error come from
			% misc/test/src/C++/hierarchicalFactorization.cpp
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'), 0.26844, 'AbsTol', 0.00001)
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
			F = matfaust.fact.hierarchical(M, params)
			this.verifyEqual(size(F), size(M))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-M;
			disp('err: ')
			norm(E,'fro')/norm(M,'fro')
			% matrix to factorize and reference relative error come from
			% misc/test/src/C++/hierarchicalFactorization.cpp
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'),0.27353, 'AbsTol', 0.0001)
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
			% misc/test/src/C++/GivensFGFT.cpp.in
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
			err = norm(E,'fro')/norm(Lap,'fro')
			% the error reference is from the C++ test,
			% misc/test/src/C++/GivensFGFT.cpp.in
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
			% misc/test/src/C++/GivensFGFTParallel.cpp.in
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
			% misc/test/src/C++/GivensFGFTParallel.cpp.in
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
			params.fact_side = 0 % forced
			params.verbose = 0 % forced
			params.init_lambda = 128;
			params = ParamsHierarchical(fact_cons, res_cons, stop_crit, stop_crit2, 'is_fact_side_left', params.fact_side == 1, 'is_update_way_R2L', params.update_way == 1, 'init_lambda', params.init_lambda, 'step_size', params.stepsize, 'constant_step_size', false, 'is_verbose', params.verbose ~= 1);
			diag_init_D = diag(init_D)
			[F,D,lambda] = matfaust.fact.fgft_palm(U, Lap, params, diag_init_D)
			this.verifyEqual(size(F), size(U))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-U;
			err = norm(E,'fro')/norm(U,'fro')
			% matrix to factorize and reference relative error come from
			% misc/test/src/C++/hierarchicalFactorizationFFT.cpp
			this.verifyEqual(err, 0.05550, 'AbsTol', 0.00001)
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
			this.verifyEqual(norm(fF-fftI), 0, 'AbsTol', 10^-12);
			this.assertEqual(full(dft(n)), full(normalize(F)), 'AbsTol', 10^-7)

			this.assertEqual(full(dft(n, 'normed', false)), full(F), 'AbsTol', 10^-7)
		end

		function test_sp(this)
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
			import matfaust.proj.skperm
			M = [-0.04440802, -0.17569296, -0.02557815, -0.15559154;
			-0.0083095,  -3.38725936, -0.78484126, -0.4883618;
			-1.48942563, -1.71787215, -0.84000212, -3.71752454;
			-0.88957883, -0.19107863, -5.92900636, -6.51064175]
			p = skperm(size(M), 2)
			pM_ref = [-0.04440802,0.,-0.02557815,0.;
			-0.0083095,-3.38725936,0.,0.;
			0.,-1.71787215,0.,-3.71752454;
			0.,0.,-5.92900636,-6.51064175]
			pM = p(M)
			this.verifyEqual(norm(pM-pM_ref)/norm(pM_ref), 0)
		end

		function test_butterfly(this)
			import matfaust.wht
			import matfaust.dft
			import matfaust.fact.butterfly
			D = full(dft(64));
			H = full(wht(64));
			dirs = {'left', 'right'}
			for t=1:2
				if(t == 2)
					M = D;
				else
					M = H;
				end
				for i=1:2
					dir = dirs{i};
					F = butterfly(M, 'dir', dir);
					this.verifyEqual(norm(full(F)-M)/norm(M), 0, 'AbsTol', 1e-6);
				end
			end
		end

		function test_palm4msa_mhtp(this)
			import matfaust.fact.palm4msa_mhtp
			import matfaust.factparams.*
			import matfaust.proj.*
			M = rand(500,32);
			projs = { splin([500,32], 5), normcol([32,32], 1.0)};
			stop_crit = StoppingCriterion(200);
			param = ParamsPalm4MSA(projs, stop_crit);
			mhtp_param = MHTPParams('num_its', 60, 'palm4msa_period', 10);
			G = palm4msa_mhtp(M, param, mhtp_param)
			err = norm(full(G)-M)/norm(M)
			this.verifyLessThanOrEqual(err, 0.20)
		end

		function test_hierarchical_mhtp(this)
			import matfaust.fact.hierarchical_mhtp
			import matfaust.factparams.ParamsHierarchical
			import matfaust.factparams.StoppingCriterion
			import matfaust.factparams.MHTPParams
			import matfaust.proj.*
			M = rand(500,32);
			fact_projs = { splin([500,32], 5), sp([32,32], 96), sp([32, 32], 96)};
			res_projs = { normcol([32,32], 1), sp([32,32], 666), sp([32, 32], 333)};
			stop_crit1 = StoppingCriterion(200)
			stop_crit2 = StoppingCriterion(200)
			% 50 iterations of MHTP will run every 100 iterations of PALM4MSA (each time PALM4MSA is called by the hierarchical algorithm)
			mhtp_param = MHTPParams('num_its', 150, 'palm4msa_period', 100)
			param = ParamsHierarchical(fact_projs, res_projs, stop_crit1, stop_crit2)
			F = hierarchical_mhtp(M, param, mhtp_param)
			err = norm(full(F)-M)/norm(M)
			this.verifyLessThanOrEqual(err, 0.20)
		end

		function test_hierarchical_dft(this)
			disp('Test hierarchical dft')
			import matfaust.fact.hierarchical
			import matfaust.dft
			DFT = full(dft(32))
			F = hierarchical(DFT, 'dft', 'backend', 2020)
			err = norm(full(F)-DFT)/norm(DFT)
			this.verifyLessThanOrEqual(err, 1e-3)
			DFT = full(dft(32))
			F = hierarchical(DFT, 'dft', 'backend', 2016)
			err = norm(full(F)-DFT)/norm(DFT)
			this.verifyLessThanOrEqual(err, 1e-3)
		end

	end

	methods
		function faust_test = FaustFactoryTest(varargin)
			faust_test.faust_paths = varargin
			run(faust_test)
		end
	end

end

