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
			import matfaust.Faust
		end

		function instantiateTestFaust(this)
			import matfaust.Faust
			num_factors = randi(this.MAX_NUM_FACTORS)
			density = 0.1
			d2 = randi(this.MAX_DIM_SIZE);
			factors = cell(1,num_factors);
			for i = 1:num_factors
				d1 = d2;
				d2 = randi(this.MAX_DIM_SIZE);
				factors{i} = sprand(d1, d2, density);
			end
			this.test_faust = Faust(factors,1.0)
            this.factors = factors;
            this.num_factors = num_factors
		end
    end

	methods(TestMethodTeardown)
		function deleteTestFaust(this)
			delete(this.test_faust)
		end
	end

	methods(Static)
		function sumnnz = nnzero_count(factors)
			sumnnz = 0
			disp('size:');size(factors,2)
			for i = 1:size(factors,2)
				sumnnz = sumnnz + nnz(factors{i})
			end
		end
	end

	methods (Test)
        function testSave(this)
			import matfaust.Faust
			rand_suffix = int2str(randi(10000))
			filepath = [ tempdir filesep 'test_faust' rand_suffix '.mat']
			save(this.test_faust, filepath)
			F = Faust(filepath)
			for i = 1:get_num_factors(F)
				this.verifyEqual(full(get_factor(this.test_faust,i)), full(get_factor(F,i)))
			end
			delete(filepath)
        end

		function testSave2(this)
			import matfaust.Faust
			rand_suffix = int2str(randi(10000))
			filepath_ref = [ tempdir filesep 'ref_faust' rand_suffix '.mat']
			filepath_test = [ tempdir filesep 'test_faust' rand_suffix '.mat']
			nb_fact=get_num_factors(this.test_faust);

			faust_factors=cell(1,nb_fact);

			for i=1:nb_fact
				faust_factors{i}=get_factor(this.test_faust,i);
			end
			save(filepath_ref,'faust_factors');
			save(this.test_faust,filepath_test)
			ref_F = Faust(filepath_ref)
			test_F = Faust(filepath_test)
			this.verifyEqual(get_num_factors(ref_F), get_num_factors(test_F))
			for i = 1:get_num_factors(ref_F)
				this.verifyEqual(full(get_factor(ref_F,i)), full(get_factor(test_F,i)))
			end

			delete(filepath_ref)
			delete(filepath_test)
		end

        function testSize(this)
           disp('testSize()')
           [num_rows, num_cols] = size(this.test_faust);
           num_facs = this.num_factors;
           size_first_fac = size(this.factors{1});
           size_last_fac = size(this.factors{num_facs});
           this.verifyEqual(num_rows, size_first_fac(1))
           this.verifyEqual(num_cols, size_last_fac(2))
        end

		function testNumel(this)
			disp('testNumel()')
			test_numel = numel(this.test_faust);
			ref_numel = numel(this.mulFactors());
			this.verifyEqual(test_numel, ref_numel)
		end

        function testGetFactorAndConstructor(this)
            disp('testGetFactorAndConstructor()')
            for i = 1:this.num_factors
               this.verifyEqual(full(get_factor(this.test_faust, i)), full(this.factors{i}))
            end
			% test get_factor on transpose Faust
			tF = this.test_faust.'
			for i = 1:this.num_factors
               this.verifyEqual(full(transpose(get_factor(tF, this.num_factors-i+1))), full(this.factors{i}))
            end
			% test get_factor on conj Faust
			cF = conj(this.test_faust)
			for i = 1:this.num_factors
               this.verifyEqual(full(conj(get_factor(cF, i))), full(this.factors{i}))
            end
			% test get_factor on transpose Faust
			tcF = this.test_faust.'
			for i = 1:this.num_factors
               this.verifyEqual(full(conj(transpose(get_factor(tcF, this.num_factors-i+1)))), full(this.factors{i}))
		   end
		end

        function testGetNumFactors(this)
            disp('testGetNumFactors()')
            this.verifyEqual(get_num_factors(this.test_faust), this.num_factors)
        end

        function testNorm2(this)
            disp('testNorm2()')
            ref_F = this.mulFactors();
            this.verifyEqual(norm(ref_F,2), norm(this.test_faust, 2), 'RelTol',0.05)
        end

        function testNorm1(this)
            disp('testNorm1()')
            ref_F = this.mulFactors();
            this.verifyEqual(norm(ref_F,1), norm(this.test_faust, 1), 'RelTol',0.05)
        end

        function testNormFro(this)
            disp('testNormFro()')
            ref_F = this.mulFactors();
            this.verifyEqual(norm(ref_F, 'fro'), norm(this.test_faust, 'fro'), 'RelTol',0.05)
        end

        function testNormInf(this)
            disp('testNormInf()')
            ref_F = this.mulFactors();
            this.verifyEqual(norm(ref_F, inf), norm(this.test_faust, inf), 'RelTol',0.05)
		end

        function testnnz(this)
            disp('testnnz()')
			this.verifyEqual(nnz_sum(this.test_faust), FaustTest.nnzero_count(this.factors))
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
            this.verifyEqual(ref_nlines*ref_ncols/FaustTest.nnzero_count(this.factors), rcg(this.test_faust), 'RelTol', 0.01)
        end

        function testend(this)
            disp('testend()')
            prod_F = this.mulFactors();
            for i=1:size(this.factors{1},1)
				end_F = full(this.test_faust(i,end));
				this.verifyEqual(prod_F(i,end), end_F(1,1),'RelTol', 0.01)
            end
            for j=1:size(this.factors{this.num_factors},2)
				end_F = full(this.test_faust(end,j));
				this.verifyEqual(prod_F(end,j), end_F(1,1), 'RelTol', 0.01)
            end
        end

        function testsubsref(this)
            disp('testsubsref()')
            % test whole faust
            ref_F =  this.mulFactors();
            test_F = full(this.test_faust(1:end,1:end));
            this.verifyEqual(ref_F(1:end,1:end), test_F, 'RelTol', 0.01)
            % test a random row
            row_i = randi([1,size(this.test_faust,1)]);
            test_F = full(this.test_faust(row_i,:));
            this.verifyEqual(test_F, ref_F(row_i,:), 'RelTol', 0.01)
            % test a random col
            col_i = randi([1,size(this.test_faust,2)]);
            test_F = full(this.test_faust(:,col_i));
            this.verifyEqual(test_F, ref_F(:,col_i), 'RelTol', 0.01)
			% test double-slice (on rows and cols at the same time)
			tested_fausts = {this.test_faust, ctranspose(this.test_faust)}
			for i=1,length(tested_fausts)
				row_i1 = randi([1,size(tested_fausts{i},1)]);
				row_i2 = randi([row_i1,size(tested_fausts{i},1)]);
				col_i1 = randi([1,size(tested_fausts{i},2)]);
				col_i2 = randi([col_i1,size(tested_fausts{i},2)]);
				% slice test faust, affect the resulting Faust to G and compare the full matrices
				full_test_F = full(tested_fausts{i});
				G = tested_fausts{i}(row_i1:row_i2, col_i1:col_i2);
				full_G = full(G);
				this.verifyEqual(full_test_F(row_i1:row_i2, col_i1:col_i2),full_G,'RelTol', 0.0001)
				% now slice G to get H
				row_i1 = randi([1,size(G,1)]);
				row_i2 = randi([row_i1,size(G,1)]);
				col_i1 = randi([1,size(G,2)]);
				col_i2 = randi([col_i1,size(G,2)]);
				H = G(row_i1:row_i2, col_i1:col_i2);
				full_H = full(H)
				% ensure consistence of slicing on full matrices
				this.verifyEqual(full_G(row_i1:row_i2, col_i1:col_i2), full_H, 'RelTol', 0.0001)
				% test second slice on transpose of sliced Faust G
				tG = transpose(G)
				full_tG = full(tG)
				I = tG(col_i1:col_i2, row_i1:row_i2)
				full_I = full(I)
				this.verifyEqual(full_I, full_tG(col_i1:col_i2, row_i1:row_i2), 'RelTol', 0.0001)
				%============= test indexing by arbitrary vectors
				F = tested_fausts{i};
				num_rows = randi([1,size(F,1)*2]);
				num_cols = randi([1,size(F,2)*2]);
				row_ids = zeros(1,num_rows);
				col_ids = zeros(1,num_cols);
				for j=1:num_rows
					row_ids(1,j) = randi([1,size(F,1)]);
				end
				for j=1:num_cols
					col_ids(1,j) = randi([1,size(F,2)]);
				end
				row_ids = row_ids(1,1:num_rows);
				col_ids = col_ids(1,1:num_cols);
				fF = full(F);
				% test sub-indexing on rows, all columns;
				ref_sfF = fF(row_ids,:);
				test_sF = full(F(row_ids,:));
				this.verifyEqual(ref_sfF, test_sF, 'RelTol', 0.0001);
				% test sub-indexing on cols, all rows;
				ref_sfF = fF(:,col_ids);
				test_sF = full(F(:,col_ids));
				this.verifyEqual(ref_sfF, test_sF, 'RelTol', 0.0001);
				% test sub-indexing on rows and cols;
				ref_sfF = fF(row_ids,col_ids);
				test_sF = full(F(row_ids,col_ids));
				this.verifyEqual(ref_sfF, test_sF, 'RelTol', 0.0001);
				%============= test slice indexing with negative steps;
				row_step = randi([row_i1-row_i2-1, -1])
				col_step = randi([col_i1-col_i2-1, -1])
				ref_sfF = fF(row_i2:row_step:row_i1,col_i2:col_step:col_i1);
				test_sF = full(F(row_i2:row_step:row_i1,col_i2:col_step:col_i1));
				this.verifyEqual(ref_sfF, test_sF, 'RelTol', 0.0001);
			end
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
            %TODO
        end

        function testConj(this)
            disp('Test Faust.conj()')
            %TODO
        end


		function testMul(this)
			disp('Test Faust.mtimes()')
			F = this.test_faust;
			rmat = rand(size(F,2),size(F,2));
			ref_full_faust = this.mulFactors();
			ref_mat = ref_full_faust*rmat;
			test_mat = F*rmat;
			this.verifyEqual(test_mat,ref_mat, 'RelTol', 10^-3)
			% do the same for a complex matrix
			cmat = rand(size(F,2),size(F,2)) + i*rand(size(F,2),size(F,2));
			ref_mat = ref_full_faust*cmat;
			test_mat = F*cmat;
			this.verifyEqual(test_mat,ref_mat, 'RelTol', 10^-3)
			% mul by real scalar
			r = rand();
			test_rF = full(F*r);
			test_commu_rF = full(r*F)
			ref_rF = ref_full_faust*r;
			this.verifyEqual(test_rF,ref_rF, 'RelTol', 10^-3);
			this.verifyEqual(test_commu_rF,ref_rF, 'RelTol', 10^-3);
			this.verifyNotEqual(test_rF, ref_full_faust*(r+1))
			this.verifyNotEqual(test_commu_rF, ref_full_faust*(r+1))
			this.assertLessThan(norm(full(F'*r)-full(F)'*r)/norm(full(F)'*r), eps(1.))
			this.assertLessThan(norm(full(F.'*r)-full(F).'*r)/norm(full(F).'*r), eps(1.))
			% TODO: mul by complex scalar (when impl.)
		end

		function testcat(this)
			import matfaust.*
			disp('Test cat')
			FAUST=0
			SPARSE=1
			FULL=2
			for typeG=0:2
				for dimcat=1:2
					other_dim = mod(dimcat,2)+1;
					F = this.test_faust;
					%=============== test vert (or horz) cat
					G = FaustFactory.rand(randi(FaustTest.MAX_NUM_FACTORS), size(F,other_dim));
					G_num_factors = get_num_factors(G);
					% set a Faust with a random number of rows (or cols) from G
					H_facs = cell(1,G_num_factors+1);
					if dimcat == 1
						H_facs{1} = rand(randi(FaustTest.MAX_DIM_SIZE-1)+1,size(F,other_dim));
						for i=1:G_num_factors
							H_facs{i+1} = get_factor(G, i);
						end
					else
						H_facs{G_num_factors+1} = rand(size(F,other_dim),randi(FaustTest.MAX_DIM_SIZE-1)+1);
						for i=1:G_num_factors
							H_facs{i} = get_factor(G, i);
						end
					end
					H = Faust(H_facs);
					if(typeG == SPARSE)
						H_ = sparse(full(H));
					elseif(typeG == FAUST)
						H_ = H;
					else % typeG == FULL
						H_ = full(H);
					end
					if dimcat == 1
						C = vertcat(F,H_);
						D = [F;H_];
					else
						C = horzcat(F,H_);
						D = [F H_];
					end
					this.verifyEqual(full(C), full(D), 'AbsTol', 1e-11)
					this.verifyEqual(full(C), cat(dimcat, full(F),full(H)),'AbsTol', 1e-11)
					C = cat(dimcat,F,H)
					this.verifyEqual(full(C), cat(dimcat, full(F),full(H)),'AbsTol', 1e-11)
				end
			end
		end

        function testDelete(this)
            disp('Test Faust.delete()')
            tFaust = transpose(this.test_faust);
            delete(this.test_faust);
            this.verifyError(@() size(this.test_faust),'MATLAB:class:InvalidHandle')
            this.verifyEqual(size(tFaust), [size(this.factors{this.num_factors},2), size(this.factors{1},1)])
        end

    end

	methods
		function faust_test = FaustTest(varargin)
			faust_test.faust_paths = varargin
			run(faust_test)
		end

        function prod_F = mulFactors(this)
            first_fac_size = size(this.factors{1});
            prod_F = eye(first_fac_size(1));
            for i=1:this.num_factors
                prod_F = prod_F*this.factors{i};
            end
        end

    end
end


