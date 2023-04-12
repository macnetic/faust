%==================================================
%> @brief Functor for the SYMM_SP projector.
%>
%> A, the image matrix, is such that A is symmetric and \f$ \| A \|_0 = k,  \| A\|_F = 1 \f$ (if normalized == True).
%==================================================
classdef symm_sp < matfaust.proj.proj_gen
	properties
	end
	methods
		%===============================================
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b symm_sp(@b shape, @b k): returns a SYMM_SP projector (functor), shape defines the size of the input matrix (e.g. [1, 10]), k defines the sparsity of the output matrix (k nnz coefficients).<br/>
		%> &nbsp;&nbsp;&nbsp; @b symm_sp(@b shape,@b k,@b 'normalized', bool,@b 'pos', bool): the optional parameters are set. By default both normalized and pos are false.
		%>
		%> @param shape: vector of size 2, to define the size of the input matrix.
		%> @param k: the sparsity parameter (the number of nonzeros of the projection image.
		%> @param 'normalized', true: normalizes the projection image according to its Frobenius norm.
		%> @param 'normalized', false: (the default) no normalization.
		%> @param 'pos', true: skips the negative values (replaced by zero) of the input matrix.
		%> @param 'pos', false: (the default) negative values are not skipped.
        %>
		%> @retval proj the sp projector.
        %>
        %> <br/>
        %>
        %> @b Example
        %> @code
		%> >> import matfaust.proj.symm_sp
		%> >> rng(42) % for reproducibility
		%> >> M = rand(5, 5);
		%> >> p = symm_sp(size(M), 2)
		%> >> p(M)
		%>
		%> ans =
		%>
		%>          0         0         0         0         0
		%>          0         0    0.7071         0         0
		%>          0    0.7071         0         0         0
		%>          0         0         0         0         0
		%>          0         0         0         0         0
        %> @endcode
		%===============================================
		function proj = symm_sp(shape, k, varargin)
			import matfaust.factparams.ConstraintInt
			% default values
			proj.constraint = ConstraintInt('symm_sp', shape(1), shape(2), k, varargin{:});
		end
	end
end
