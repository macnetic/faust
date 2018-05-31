classdef FaustTest < matlab.unittest.TestCase
    % Tests the Faust class.

	properties
		test_faust
		faust_paths
	end

	properties (Constant = true)
		MAX_NUM_FACTORS = 64
		MAX_DIM_SIZE = 1024
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
		end

		function instantiateTestFaust(this)
			num_factors = randi(this.MAX_NUM_FACTORS)
			density = 0.1
			d2 = randi(this.MAX_DIM_SIZE);
			factors = cell(1,num_factors);
			for i = 1:num_factors
				d1 = d2;
				d2 = randi(this.MAX_DIM_SIZE);
				factors{i} = sprand(d1, d2, density);
			end
			this.test_faust = Faust(factors,1.0,0)
		end
    end

	methods(TestMethodTeardown)
		function deleteTestFaust(this)
			delete(this.test_faust)
		end
	end

	methods (Test)
        function testSave(this)
			rand_suffix = int2str(randi(10000))
			filepath = [ tempdir filesep 'test_faust' rand_suffix '.mat']
			save(this.test_faust, filepath)
			F = Faust(filepath)
			for i = 1:get_nb_factor(F)
				this.verifyEqual(get_fact(this.test_faust,i), get_fact(F,i))
			end
			delete(filepath)
        end

		function testSave2(this)
			rand_suffix = int2str(randi(10000))
			filepath_ref = [ tempdir filesep 'ref_faust' rand_suffix '.mat']
			filepath_test = [ tempdir filesep 'test_faust' rand_suffix '.mat']
			nb_fact=get_nb_factor(this.test_faust);

			faust_factors=cell(1,nb_fact);

			for i=1:nb_fact
				faust_factors{i}=get_fact(this.test_faust,i);
			end
			save(filepath_ref,'faust_factors');
			save(this.test_faust,filepath_test)
			ref_F = Faust(filepath_ref)
			test_F = Faust(filepath_test)
			this.verifyEqual(get_nb_factor(ref_F), get_nb_factor(test_F))
			for i = 1:get_nb_factor(ref_F)
				this.verifyEqual(get_fact(ref_F,i), get_fact(test_F,i))
			end

			delete(filepath_ref)
			delete(filepath_test)
		end
    end

	methods
		function faust_test = FaustTest(varargin)
			faust_test.faust_paths = varargin
			run(faust_test)
		end
	end
end


