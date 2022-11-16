classdef LazyLinearOpFFTTest < LazyLinearOpTest

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
            import matfaust.lazylinop.*
			this.lop = LazyLinearOp.create_from_funcs(@(x) fft(x), @(x) 8 * ifft(x), [8, 8]);
            this.from_func = true;
            this.lopA = full(this.lop);
            F2 = rand(8, 8);
			this.lop2 = aslazylinearoperator(F2);
            this.lop2A = full(F2);
            F3 = rand(8, 8);
			this.lop3 = aslazylinearoperator(F3);
            this.lop3A = full(F3);
		end
    end

	methods(TestMethodTeardown)

	end

	methods(Static)

	end


end
