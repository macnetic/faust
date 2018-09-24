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
			import matfaust.*;
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
			params = ParamsPalm4MSA(num_facts, is_update_way_R2L, init_lambda, cons, stop_crit, init_facts, ParamsFact.DEFAULT_STEP_SIZE, ParamsFact.DEFAULT_CONSTANT_STEP_SIZE);
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
			import matfaust.*;
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
			params = ParamsPalm4MSA(num_facts, is_update_way_R2L, init_lambda, cons, stop_crit, init_facts, ParamsFact.DEFAULT_STEP_SIZE, ParamsFact.DEFAULT_CONSTANT_STEP_SIZE);
			F = FaustFactory.fact_palm4msa(M, params)
			this.verifyEqual(size(F), size(M))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-M;
			disp('err: ')
			norm(E,'fro')/norm(M,'fro')
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'), 0.9094, 'AbsTol', 0.0001)
		end

		function test_fact_hierarchical(this)
			disp('Test FaustFactory.fact_hierarchical()')
			import matfaust.*;
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
			params = ParamsHierarchicalFact(num_facts, is_update_way_R2L, init_lambda,...
				fact_cons, res_cons, size(M,1), size(M,2), {stop_crit, stop_crit2});
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
			import matfaust.*;
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
			params = ParamsHierarchicalFact(num_facts, is_update_way_R2L, init_lambda,...
				fact_cons, res_cons, size(M,1), size(M,2), {stop_crit, stop_crit2});
			F = FaustFactory.fact_hierarchical(M, params)
			this.verifyEqual(size(F), size(M))
			%disp('norm F: ')
			%norm(F, 'fro')
			E = full(F)-M;
			disp('err: ')
			norm(E,'fro')/norm(M,'fro')
			% matrix to factorize and reference relative error come from
			% misc/test/src/C++/hierarchicalFactorization.cpp
			this.verifyEqual(norm(E,'fro')/norm(M,'fro'), 0.99275, 'AbsTol', 0.00001)
		end	end

	methods
		function faust_test = FaustFactoryTest(varargin)
			faust_test.faust_paths = varargin
			run(faust_test)
		end
	end

end

