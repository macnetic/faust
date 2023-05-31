% ==================================================
%> @brief Functor for the SPTRIL projector.
%>
%> A, the image matrix, is such that the upper triangular part is 0 and \f$ \| A \|_0 = k,  \| A\|_F = 1 \f$ (if normalized == True).
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b sptril(@b shape, @b k): returns a SPTRIL projector (functor), shape defines the size of the input matrix (e.g. [1, 10]), k defines the sparsity of the output matrix (k nnz coefficients).<br/>
%> &nbsp;&nbsp;&nbsp; @b sptril(@b shape,@b k,@b 'normalized', bool,@b 'pos', bool): the optional parameters are set. By default both normalized and pos are false.
%>
%> @param shape: vector of size 2, to define the size of the input matrix.
%> @param k: the sparsity parameter (the number of nonzeros of the projection image.
%> @param 'normalized', true: normalizes the projection image according to its Frobenius norm.
%> @param 'normalized', false: (the default) no normalization.
%> @param 'pos', true: skips the negative values (replaced by zero) of the input matrix.
%> @param 'pos', false: (the default) negative values are not skipped.
%>
%> @retval proj the sptril projector.
%>
%> <br/>
%>
%> @b Example
%> @code
%>
%> >> import matfaust.proj.sptril
%> >> rng(42) % reproducibility
%> >> M = rand(5, 5);
%> >> p = sptril(size(M), 2);
%> >> p(M)
%>
%> ans =
%>
%>          0         0         0         0         0
%>     0.7392         0         0         0         0
%>          0    0.6735         0         0         0
%>          0         0         0         0         0
%>          0         0         0         0         0
%>
%> >>
%> @endcode
%>
%> @b see @b also matfaust.proj.sptriu
% ==================================================
classdef sptril < matfaust.proj.proj_gen
	properties
	end
	methods

		function proj = sptril(shape, k, varargin)
			import matfaust.factparams.ConstraintInt
			% default values
			proj.constraint = ConstraintInt('sptril', shape(1), shape(2), k, varargin{:});
		end
	end
end
