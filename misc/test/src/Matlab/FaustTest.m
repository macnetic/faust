classdef FaustTest < matlab.unittest.TestCase
    % Tests the Faust class.

	properties
		test_faust
		faust_paths
        factors
        num_factors
    end

	properties (Constant = true)
		MAX_NUM_FACTORS = 32
		MAX_DIM_SIZE = 512
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
            this.factors = factors
            this.num_factors = num_factors
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

        function testSize(this)
           disp('testSize()')
           [num_rows, num_cols] = size(this.test_faust)
           num_facs = this.num_factors
           size_first_fac = size(this.factors{1})
           size_last_fac = size(this.factors{num_facs})
           this.verifyEqual(num_rows, size_first_fac(1))
           this.verifyEqual(num_cols, size_last_fac(2))
        end

        function testGetFactorAndConstructor(this)
            disp('testGetFactorAndConstructor()')
            for i = 1:this.num_factors
               this.verifyEqual(get_fact(this.test_faust, i), full(this.factors{i}))
            end
        end

        function testGetNumFactors(this)
            disp('testGetNumFactors()')
            this.verifyEqual(get_nb_factor(this.test_faust), this.num_factors)
        end

        function testNorm2(this)
            disp('testNorm2()')
            ref_F = this.mulFactors()
            this.verifyEqual(norm(ref_F,2), norm(this.test_faust, 2), 'RelTol',0.05)
        end

        function testnnz(this)
            disp('testnnz()')
            this.verifyEqual(nnz(this.test_faust), nnzero_count(this.factors))
        end

        function testDensity(this)
            disp('testDensity()')
            ref_nlines = size(this.factors{1},1)
            ref_ncols = size(this.factors{this.num_factors},2)
            this.verifyEqual(nnzero_count(this.factors)/ref_nlines/ref_ncols, density(this.test_faust), 'AbsTol', 0.01)
        end

        function testrcg(this)
            disp('testrcg()')
            ref_nlines = size(this.factors{1},1)
            ref_ncols = size(this.factors{this.num_factors},2)
            this.verifyEqual(ref_nlines*ref_ncols/nnzero_count(this.factors), RCG(this.test_faust), 'RelTol', 0.01)
        end

        function testend(this)
            disp('testend()')
            prod_F = this.mulFactors()
            for i=1:size(this.factors{1},1)
                this.verifyEqual(prod_F(i,end), this.test_faust(i,end),'RelTol', 0.01)
            end
            for j=1:size(this.factors{this.num_factors},2)
                this.verifyEqual(prod_F(end,j), this.test_faust(end,j), 'RelTol', 0.01)
            end
        end

        function testsubsref(this)
            disp('testsubsref()')
            % test whole faust
            ref_F =  this.mulFactors();
            test_F = this.test_faust(1:end,1:end);
            this.verifyEqual(ref_F(1:end,1:end), test_F, 'RelTol', 0.01)
            % test a random row
            row_i = randi(1,size(this.test_faust,1));
            test_F = this.test_faust(row_i,:);
            this.verifyEqual(test_F, ref_F(row_i,:), 'RelTol', 0.01)
            % test a random col
            col_i = randi(1,size(this.test_faust,2));
            test_F = this.test_faust(:,col_i);
            this.verifyEqual(test_F, ref_F(:,col_i), 'RelTol', 0.01)
        end

        function testFull(this)
            disp('Test Faust.full()')
            ref_faust = this.mulFactors();
            test_faust = full(this.test_faust);
            this.verifyEqual(test_faust, ref_faust, 'RelTol', 10^-3)
        end

        function testTranspose(this)
            disp('Test Faust.transpose()')
            % test full matrix to avoid longer time for comparisons (element by
            % element)
            transp_faust = full(transpose(this.test_faust)); 
            test_faust = full(this.test_faust);
            this.verifyEqual(test_faust,transpose(transp_faust), 'RelTol', 10^-3)
        end

        function testCtranspose(this)
            disp('Test Faust.ctranspose()')
            %TODO (when the function will be implemeted)
        end

        function testConj(this)
            disp('Test Faust.conj()')
            %TODO (when the funciton will be implemented)
        end
    end

	methods
		function faust_test = FaustTest(varargin)
			faust_test.faust_paths = varargin
			run(faust_test)
		end

        function prod_F = mulFactors(this)
            first_fac_size = size(this.factors{1})
            prod_F = eye(first_fac_size(1));
            for i=1:this.num_factors
                prod_F = prod_F*this.factors{i};
            end
        end

    end
end


