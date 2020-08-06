%==================================================
%> @brief Functor for the SPCOL projector.
%>
%> A, the image matrix, is defined by \f$ \forall j \in \{1,\ldots,shape(2)\} \f$ the j-th column \f$  A_{*,j} \f$ is such that \f$ \|A_{*,j}\|_0 = k,  \| A\|_F = 1 \f$ (if normalized is true).
%==================================================
classdef spcol < matfaust.proj.proj_gen
	properties
	end
	methods
		%===============================================
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b spcol(shape, k): returns a SPCOL projector (functor), shape defines the size of the input matrix (e.g. [1, 10]), k defines the sparsity of the output matrix (k nnz coefficients per column).<br/>
		%> &nbsp;&nbsp;&nbsp; @b spcol(@b shape,@b k,@b 'normalized', bool,@b 'pos', bool): the optional parameters are set. By default both normalized and pos are false.
		%>
		%> @param shape: vector of size 2, to define the size of the input matrix.
		%> @param k: the sparsity parameter.
		%> @param 'normalized', true: normalizes the projection image according to its Frobenius norm.
		%> @param 'normalized', false: (the default) no normalization.
		%> @param 'pos', true: skips the negative values (replaced by zero) of the input matrix.
		%> @param 'pos', false: (the default) negative values are not skipped.
		%> @retval the spcol projector.
		%===============================================
		function proj = spcol(shape, k, varargin)
			import matfaust.factparams.ConstraintInt
			proj.constraint = ConstraintInt('spcol', shape(1), shape(2), k, varargin{:});
		end
	end
end
