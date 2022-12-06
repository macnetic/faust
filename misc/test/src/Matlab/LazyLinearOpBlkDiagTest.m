classdef LazyLinearOpBlkDiagTest < LazyLinearOpTest

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
            import matfaust.rand
            A = rand(12, 14);
            B = rand(13, 15);
            C = rand(14, 16);
			this.lop = matfaust.lazylinop.blkdiag(A, B, C)
            this.lopA = blkdiag(full(A), full(B), full(C));
			this.lop2 = matfaust.lazylinop.blkdiag(A+A, B+B, C+C);
            this.lop2A = blkdiag(2*full(A), 2*full(B), 2*full(C));
			this.lop3 = matfaust.lazylinop.blkdiag(C.', B.', A.');
            this.lop3A = blkdiag(full(C.'), full(B.'), full(A.'));
		end
    end

	methods(TestMethodTeardown)

	end

end
