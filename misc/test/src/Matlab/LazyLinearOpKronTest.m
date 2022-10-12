classdef LazyLinearOpKronTest < LazyLinearOpTest

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
            import matfaust.lazylinop.asLazyLinearOp
            import matfaust.lazylinop.kron
            F1 = rand(10, 15);
            F2 = rand(10, 15);
			lop_A = asLazyLinearOp(F1);
			lop_B = asLazyLinearOp(F2);
            this.lop = kron(lop_A, lop_B)
            this.lopA = full(this.lop);
			this.lop2 = asLazyLinearOp(rand(size(this.lop, 1), size(this.lop, 2)));
            this.lop2A = full(this.lop2);
			this.lop3 = asLazyLinearOp(rand(size(this.lop, 2), 10));
            this.lop3A = full(this.lop3);
		end
    end

	methods(TestMethodTeardown)

	end

	methods(Static)

	end


	methods
%		function lop_test = LazyLinearOpKronTest(varargin)
%			lop_test.faust_paths = varargin
%			%ENOTE char arrays concat doesn't work without space or comma separator between arrays
%			try % try to call only a single unit test
%				run(lop_test, varargin{end})
%			catch
%				run(lop_test)
%			end
%		end
    end
end
