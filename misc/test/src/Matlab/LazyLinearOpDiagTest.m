classdef LazyLinearOpDiagTest < LazyLinearOpTest

	properties (Constant = true)
        %
	end

	methods(TestClassSetup)
	end

    methods (TestMethodSetup)

        function addFaustToPath(this)
			addpath(this.faust_paths{:})
        	set_path
		end

		function instantiateTestFaust(this)
            import matfaust.lazylinop.aslazylinearoperator
            import matfaust.lazylinop.diag
            v = rand(10, 1);
            w = rand(10, 1);
			this.lop = diag(v, -2);
			this.lop2 = diag(v, 2);
            this.lopA = full(this.lop);
            this.lop2A = full(this.lop2);
			this.lop3 = diag(w, 2);
            this.lop3A = full(this.lop3);
		end
    end

	methods(TestMethodTeardown)

	end

end
