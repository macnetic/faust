%==================================================
%> @brief Functor for the SKPERM projector.
%>
%> This projector returns the matrix with k nonzeros of the input matrix per column and per row that maximizes the Frobenius norm.
%==================================================
classdef skperm < matfaust.proj.proj_gen
	properties
	end
	methods
		%===============================================
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b skperm(shape, k): returns a SKPERM projector (functor), shape defines the size of the input matrix (e.g. [1, 10]), k defines the sparsity of the output matrix (k nnz coefficients per row or per column).<br/>
		%> &nbsp;&nbsp;&nbsp; @b splincol(@b shape,@b k,@b 'normalized', bool,@b 'pos', bool): the optional parameters are set. By default both normalized and pos are false.
		%>
		%> @param shape: vector of size 2, to define the size of the input matrix.
		%> @param k: the sparsity parameter.
		%> @param 'normalized', true: (the default) normalizes the projection image according to its Frobenius norm.
		%> @param 'normalized', false: no normalization.
		%> @param 'pos', true: skips the negative values (replaced by zero) of the input matrix.
		%> @param 'pos', false: (the default) negative values are not skipped.
		%> @retval proj the skperm projector.
		%===============================================
		function proj = skperm(shape, k, varargin)
			import matfaust.factparams.ConstraintInt
			% default values
			proj.constraint = ConstraintInt('skperm', shape(1), shape(2), k, varargin{:});
		end
	end
end
