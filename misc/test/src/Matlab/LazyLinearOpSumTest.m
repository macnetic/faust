classdef LazyLinearOpSumTest < LazyLinearOpTest

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
            B = rand(12, 14);
            C = rand(12, 14);
            function S = add(varargin)
                S = 0;
                for i=1:length(varargin)
                    S = S + varargin{i};
                end
            end
			this.lop = matfaust.lazylinop.sum(A, B, C)
            this.lopA = add(full(A), full(B), full(C));
			this.lop2 = matfaust.lazylinop.sum(A+A, B+B, C+C);
            this.lop2A = add(2*full(A), 2*full(B), 2*full(C));
			this.lop3 = matfaust.lazylinop.sum(C.', B.', A.');
            this.lop3A = add(full(C.'), full(B.'), full(A.'));
		end
    end

	methods(TestMethodTeardown)

	end

end
