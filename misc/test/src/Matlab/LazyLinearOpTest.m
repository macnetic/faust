classdef LazyLinearOpTest < matlab.unittest.TestCase

	properties
		lop;
        lopA;
        lop2;
        lop2A;
        lop3;
        lop3A;
        faust_paths;
        from_func;
    end

	properties (Constant = true)
        %
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
			import matfaust.rand
            import matfaust.lazylinop.aslazylinearoperator
            F1 = rand(10, 15);
			this.lop = aslazylinearoperator(F1);
            this.from_func = false;
            this.lopA = full(F1);
            F2 = rand(10, 15);
			this.lop2 = aslazylinearoperator(F2);
            this.lop2A = full(F2);
            F3 = rand(15, 10);
			this.lop3 = aslazylinearoperator(F3);
            this.lop3A = full(F3);
		end
    end

	methods(TestMethodTeardown)

	end

	methods(Static)

	end

	methods (Test)

        function testSize(this)
            this.verifyEqual(size(this.lop), size(full(this.lop)))
        end

        function testTransp(this)
            lopT = this.lop.';
            this.verifyEqual(full(lopT), this.lopA.', 'AbsTol', 1e-6)
            this.verifyEqual(size(lopT, 1), size(this.lop, 2))
            this.verifyEqual(size(lopT, 2), size(this.lop, 1))
        end

        function testConj(this)
            lopC = conj(this.lop);
            this.verifyEqual(full(lopC), conj(this.lopA), 'AbsTol', 1e-6)
        end

        function testAdjoint(this)
            lopH = this.lop';
            this.verifyEqual(full(lopH), this.lopA', 'AbsTol', 1e-6)
            this.verifyEqual(size(lopH, 1), size(this.lop, 2))
            this.verifyEqual(size(lopH, 2), size(this.lop, 1))
        end

        function testAdd(this)
            ladd = this.lop + this.lop2
            this.verifyEqual(full(ladd), this.lopA + this.lop2A, 'AbsTol', 1e-6)
            M = rand(size(this.lop))
            ladd2 = this.lop + M
            this.verifyTrue(isa(ladd2, 'matfaust.lazylinop.LazyLinearOp'))
            this.verifyEqual(full(ladd2), this.lopA + M, 'AbsTol', 1e-6)
        end

        function testSub(this)
            ladd = this.lop - this.lop2
            this.verifyEqual(full(ladd), this.lopA - this.lop2A, 'AbsTol', 1e-6)
            M = rand(size(this.lop))
            ladd2 = this.lop - M
            this.verifyTrue(isa(ladd2, 'matfaust.lazylinop.LazyLinearOp'))
            this.verifyEqual(full(ladd2), this.lopA - M, 'AbsTol', 1e-6)
        end

        function testmtimes(this)
			import matfaust.lazylinop.*
            lmul = this.lop * this.lop3;
            this.verifyTrue(isLazyLinearOp(lmul))
            this.verifyEqual(full(lmul), this.lopA * this.lop3A, 'AbsTol', 1e-6)

            M = rand(size(this.lop, 2), 15)
            lmul2 = this.lop * M
            this.verifyFalse(isa(lmul2, 'matfaust.lazylinop.LazyLinearOp'))
            this.verifyTrue(ismatrix(lmul2))
            this.verifyEqual(full(lmul2), this.lopA * M, 'AbsTol', 1e-6)

            if ~ this.from_func
                lmul2 = this.lop * sparse(M)
                this.verifyFalse(isLazyLinearOp(lmul2))
                this.verifyTrue(ismatrix(lmul2))
                this.verifyEqual(full(lmul2), this.lopA * M, 'AbsTol', 1e-6)
            end

            lmul3 = this.lop * M(:, 1)
            this.verifyFalse(isLazyLinearOp(lmul3))
            this.verifyTrue(ismatrix(lmul3))
            this.verifyEqual(full(lmul3), this.lopA * M(:,1), 'AbsTol', 1e-6)

            lmulS = aslazylinearoperator(sparse(M)) * sparse(M.')
            this.verifyFalse(isLazyLinearOp(lmulS))
            this.verifyTrue(ismatrix(lmulS))
            this.verifyEqual(full(lmulS), M * M.', 'AbsTol', 1e-6)

            lmul4 = this.lop * 5
            this.verifyTrue(isLazyLinearOp(lmul4))
            this.verifyEqual(size(lmul4), size(this.lop))
            this.verifyEqual(full(lmul4), this.lopA * 5, 'AbsTol', 1e-6)


        end

        function testpagemtimes(this)
            M = rand(size(this.lop, 2), 15, 2, 2)
            lmul = pagemtimes(this.lop, M);
            mul_ref = pagemtimes(this.lopA , M);
            for i=1:2
                for j=1:2
                    this.verifyEqual(lmul(:,:, i, j), mul_ref(:, :, i, j), 'AbsTol', 1e-6)
                end
            end
        end

        function testCat(this)
            import matfaust.lazylinop.LazyLinearOp
            lop = this.lop;
            lop2 = this.lop2;

            % vertcat
            lcat = [lop ; lop2]; % using directly this doesn't work
            this.verifyTrue(LazyLinearOp.isLazyLinearOp(lcat))
            this.verifyEqual(full(lcat), [this.lopA ; this.lop2A], 'AbsTol', 1e-6)

            % horzcat
            lcat = [lop , lop2]; % using directly this doesn't work
            this.verifyTrue(LazyLinearOp.isLazyLinearOp(lcat))
            this.verifyEqual(full(lcat), [this.lopA, this.lop2A], 'AbsTol', 1e-6)

            % auto cat
            lcat = [lop , lop]; % using directly this doesn't work
            this.verifyTrue(LazyLinearOp.isLazyLinearOp(lcat))
            this.verifyEqual(full(lcat), [this.lopA, this.lopA], 'AbsTol', 1e-6)
        end

        function testChainOps(this)
            lchain = this.lop + this.lop2;
            lchain = lchain * this.lop3;
            lchain = 2 * lchain;
            lop3 = this.lop3;
            lop3A = this.lop3A;
            lchain = [lchain ; lop3]
            v = rand(size(lchain, 2), 1)
            this.verifyEqual(lchain * v, [2 * ((this.lopA + this.lop2A) * this.lop3A) ; lop3A] * v, 'AbsTol', 1e-6)
        end

        function testSubsref(this)
            n1 = floor(size(this.lop, 1) / 2);
            n2 = floor(size(this.lop, 2) / 2);
            lslice = this.lop(3:n1, 3:n2);
            lsliceA = this.lopA(3:n1, 3:n2);
            this.verifyEqual(full(lslice), lsliceA);
        end

        function testReal(this)
            import matfaust.lazylinop.aslazylinearoperator
            cF = matfaust.rand(10, 15, 'field', 'complex')
            lcF = aslazylinearoperator(cF)
            this.verifyEqual(full(real(lcF)), full(real(cF)), 'AbsTol', 1e-6)
        end

        function testImag(this)
            import matfaust.lazylinop.aslazylinearoperator
            cF = matfaust.rand(10, 15, 'field', 'complex')
            lcF = aslazylinearoperator(cF)
            this.verifyEqual(full(imag(lcF)), full(imag(cF)), 'AbsTol', 1e-6)
        end

        function testAsLazyLinearOperator(this)
            import matfaust.lazylinop.aslazylinearoperator
            import matfaust.lazylinop.isLazyLinearOp
            cF = matfaust.rand(10, 15, 'field', 'complex');
            lcF = aslazylinearoperator(cF);
            this.verifyTrue(isLazyLinearOp(lcF))
            this.verifyEqual(size(cF), size(lcF))
        end

    end


	methods
		function lop_test = LazyLinearOpTest(varargin)
			lop_test.faust_paths = varargin
			%ENOTE char arrays concat doesn't work without space or comma separator between arrays
			try % try to call only a single unit test
				run(lop_test, varargin{end})
			catch
				run(lop_test)
			end
		end


    end
end
