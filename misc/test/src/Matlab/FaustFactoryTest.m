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
		function test_fact_palm4msa(this)
			disp('Test FaustFactory.fact_palm4msa()')
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
			F = FaustFactory.fact_palm4msa(M, params)
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

		function test_fact_palm4msaCplx(this)
			disp('Test FaustFactory.fact_palm4msaCplx()')
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
			F = FaustFactory.fact_palm4msa(M, params)
			this.verifyEqual(size(F), size(M))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-M;
			disp('err: ')
			norm(E,'fro')/norm(M,'fro')
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'), 0.29177, 'AbsTol', 0.0001)
		end

		function test_fact_hierarchical(this)
			disp('Test FaustFactory.fact_hierarchical()')
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
			params = ParamsHierarchicalFact(fact_cons, res_cons, stop_crit, stop_crit2);
			F = FaustFactory.fact_hierarchical(M, params)
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

		function test_fact_hierarchicalCplx(this)
			disp('Test FaustFactory.fact_hierarchicalCplx()')
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
			params = ParamsHierarchicalFact(fact_cons, res_cons, stop_crit, stop_crit2,...
					'init_lambda', init_lambda, 'is_update_way_R2L', is_update_way_R2L);
			F = FaustFactory.fact_hierarchical(M, params)
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

		function test_fgft_palm(this)
			disp('Test FaustFactory.fgft_palm()')
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
			params.verbose = 1 % forced
			params.init_lambda = 128;
			params = ParamsHierarchicalFact(fact_cons, res_cons, stop_crit, stop_crit2, 'is_fact_side_left', params.fact_side == 1, 'is_update_way_R2L', params.update_way == 1, 'init_lambda', params.init_lambda, 'step_size', params.stepsize, 'constant_step_size', false, 'is_verbose', true);
			diag_init_D = diag(init_D)
			[F,lambda] = FaustFactory.fgft_palm(U, Lap, params, diag_init_D)
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
			disp('Test FaustFactory.wht()')
			import matfaust.*
			n = randi([1,6]);
			H = FaustFactory.wht(n);
			fH = full(H);
			this.verifyEqual(nnz(fH), numel(fH));
			i = 1;
			while(i<2^n-1)
				j = i + 1;
				while(j<2^n)
					this.verifyEqual(fH(i,:)*fH(j,:)',0);
					j = j + 1;
				end
				i = i + 1;
			end
		end

		function testFourier(this)
			disp('Test FaustFactory.dft()')
			import matfaust.*
			n = randi([1,10]);
			F = FaustFactory.dft(n);
			fF = full(F);
			fftI = fft(eye(2^n));
			% this.verifyEqual(nnz(fH), numel(fH));
			this.verifyEqual(norm(fF-fftI), 0, 'AbsTol',10^-12);
		end
	end

	methods
		function faust_test = FaustFactoryTest(varargin)
			faust_test.faust_paths = varargin
			run(faust_test)
		end
	end

end

