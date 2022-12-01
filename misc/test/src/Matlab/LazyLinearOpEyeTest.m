classdef LazyLinearOpEyeTest < LazyLinearOpTest

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
            import matfaust.lazylinop.eye
			this.lop = eye(10, 15, -2);
			this.lop2 = eye(10, 15);
            this.lopA = full(this.lop);
            this.lop2A = full(this.lop2);
			this.lop3 = eye(size(this.lop, 2), 10, 2);
            this.lop3A = full(this.lop3);
		end
    end

	methods(TestMethodTeardown)

	end

end
