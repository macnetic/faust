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

		function test_palm4msaCplx(this)
			disp('Test matfaust.fact.palm4msaCplx()')
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
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'), 0.29177, 'AbsTol', 0.0001)
		end

		function test_hierarchical(this)
			disp('Test matfaust.fact.hierarchical()')
			import matfaust.*
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
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'), 0.26851, 'AbsTol', 0.00001)
		end

		function test_hierarchicalCplx(this)
			disp('Test matfaust.fact.hierarchicalCplx()')
			import matfaust.*
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
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'), 1.0063, 'AbsTol', 0.0001)
		end

		function test_fgft_givens(this)
			disp('Test matfaust.fact.fgft_givens()')
			import matfaust.*
			load([this.faust_paths{1} '../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat'])
			% Lap and J available
			[F,D] = matfaust.fact.fgft_givens(Lap, 'nGivens', J);%, 0, 'verbosity', 1);
			this.verifyEqual(size(F), size(Lap))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = F*full(D)*F'-Lap;
			err = norm(E,'fro')/norm(Lap,'fro')
			% the error reference is from the C++ test,
			% misc/test/src/C++/GivensFGFT.cpp.in
			this.verifyEqual(err, 0.0326529, 'AbsTol', 0.00001)
			% verify it works the same using the eigtj() alias function
			[F2,D2] = matfaust.fact.eigtj(Lap, 'nGivens', J);%, 0, 'verbosity', 2);
			this.verifyEqual(full(F2),full(F))
			this.verifyEqual(D,D2)
		end

		function test_fgft_givens_sparse(this)
			disp('Test matfaust.fact.fgft_givens()')
			import matfaust.*
			load([this.faust_paths{1} '../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat'])
			% Lap and J available
			[F,D] = matfaust.fact.fgft_givens(sparse(Lap), 'nGivens', J);%, 0, 'verbosity', 1);
			this.verifyEqual(size(F), size(Lap))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = F*full(D)*F'-Lap;
			err = norm(E,'fro')/norm(Lap,'fro')
			% the error reference is from the C++ test,
			% misc/test/src/C++/GivensFGFT.cpp.in
			this.verifyEqual(err, 0.0326529, 'AbsTol', 0.00001)
			% verify it works the same using the eigtj() alias function
			[F2,D2] = matfaust.fact.eigtj(Lap, 'nGivens', J);%, 0, 'verbosity', 2);
			this.verifyEqual(full(F2),full(F))
			this.verifyEqual(D,D2)
		end

		function test_fgft_givens_parallel(this)
			disp('Test matfaust.fact.fgft_givens_sparse() -- parallel')
			import matfaust.*
			load([this.faust_paths{1} '../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat'])
			% Lap and J available
			t = size(Lap,1)/2;
			[F,D] = matfaust.fact.fgft_givens(Lap, 'nGivens', J, 'nGivens_per_fac', t); %, 'verbosity', 2);
			this.verifyEqual(size(F), size(Lap))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = F*full(D)*F'-Lap;
			err = norm(E,'fro')/norm(Lap,'fro')
			% the error reference is from the C++ test,
			% misc/test/src/C++/GivensFGFTParallel.cpp.in
			this.verifyEqual(err,0.0398154, 'AbsTol', 0.00001)
			% verify it works the same using the eigtj() alias function
			[F2,D2] = matfaust.fact.eigtj(Lap, 'nGivens', J, 'nGivens_per_fac', t); %, 'verbosity', 2);
			this.verifyEqual(full(F2),full(F))
			this.verifyEqual(D,D2)
		end

		function test_fgft_givens_parallel_sparse(this)
			disp('Test matfaust.fact.fgft_givens_sparse() -- parallel')
			import matfaust.*
			load([this.faust_paths{1} '../../../misc/data/mat/test_GivensDiag_Lap_U_J.mat'])
			% Lap and J available
			t = size(Lap,1)/2;
			[F,D] = matfaust.fact.fgft_givens(sparse(Lap), 'nGivens', J, 'nGivens_per_fac', t); %, 'verbosity', 2);
			this.verifyEqual(size(F), size(Lap))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = F*full(D)*F'-Lap;
			err = norm(E,'fro')/norm(Lap,'fro')
			% the error reference is from the C++ test,
			% misc/test/src/C++/GivensFGFTParallel.cpp.in
			this.verifyEqual(err, 0.0398154, 'AbsTol', 0.00001)
			% verify it works the same using the eigtj() alias function
			[F2,D2] = matfaust.fact.eigtj(Lap, 'nGivens', J, 'nGivens_per_fac', t); %, 'verbosity', 2);
			this.verifyEqual(full(F2),full(F))
			this.verifyEqual(D,D2)
		end

		function test_fgft_palm(this)
			disp('Test matfaust.fact.fgft_palm()')
			import matfaust.*
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
			this.verifyEqual(err, 0.05539, 'AbsTol', 0.00001)
		end

		function testHadamard(this)
			disp('Test matfaust.wht()')
			import matfaust.*
			log2n = randi([1,6]);
			n = 2^log2n;
			H = wht(n, false);
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
		end

		function testFourier(this)
			disp('Test matfaust.dft()')
			import matfaust.dft
			log2n = randi([1,10]);
			n = 2^log2n;
			F = dft(n, false);
			fF = full(F);
			fftI = fft(eye(n));
			% this.verifyEqual(nnz(fH), numel(fH));
			this.verifyEqual(norm(fF-fftI), 0, 'AbsTol', 10^-12);
			this.assertEqual(full(dft(n)), full(normalize(F)), 'AbsTol', 10^-7)
		end
	end

	methods
		function faust_test = FaustFactoryTest(varargin)
			faust_test.faust_paths = varargin
			run(faust_test)
		end
	end

end

