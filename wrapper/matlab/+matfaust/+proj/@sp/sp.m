%==================================================
%> @brief Functor for the SP projector.
%>
%> A, the image matrix, is such that \f$ \| A \|_0 = k,  \| A\|_F = 1 (if normalized is true).\f$.
%==================================================
classdef sp < matfaust.proj.proj_gen
	properties
	end
	methods
		%===============================================
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b sp(shape, k): returns a SP projector (functor), shape defines the size of the input matrix (e.g. [1, 10]), k defines the sparsity of the output matrix (k nnz coefficients).<br/>
		%> &nbsp;&nbsp;&nbsp; @b sp(@b shape,@b k,@b 'normalized', bool,@b 'pos', bool): the optional parameters are set. By default both normalized and pos are false.
		%>
		%> @param shape: vector of size 2, to define the size of the input matrix.
		%> @param k: the sparsity parameter.
		%> @param 'normalized', true: normalizes the projection image according to its Frobenius norm.
		%> @param 'normalized', false: (the default) no normalization.
		%> @param 'pos', true: skips the negative values (replaced by zero) of the input matrix.
		%> @param 'pos', false: (the default) negative values are not skipped.
		%> @retval the sp projector.
		%===============================================
		function proj = sp(shape, k, varargin)
			import matfaust.factparams.ConstraintInt
			% default values
			proj.constraint = ConstraintInt('sp', shape(1), shape(2), k, varargin{:});
		end
	end
end
