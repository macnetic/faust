% ==================================================
%> @brief Functor for the SP projector.
%>
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
%> @retval sp projector.
%>
%> A, the image matrix, is such that \f$ \| A \|_0 = k,  \| A\|_F = 1\f$ (if normalized is true).
%>
%> @b Example
%> @code
%> >> import matfaust.proj.sp
%> >> rng(42)
%> >> M = rand(5, 5)
%>
%> M =
%>
%>     0.3745    0.1560    0.0206    0.1834    0.6119
%>     0.9507    0.0581    0.9699    0.3042    0.1395
%>     0.7320    0.8662    0.8324    0.5248    0.2921
%>     0.5987    0.6011    0.2123    0.4319    0.3664
%>     0.1560    0.7081    0.1818    0.2912    0.4561
%>
%> >> p = sp(size(M), 3);
%> >> p(M)
%>
%> ans =
%>
%>          0         0         0         0         0
%>     0.5902         0    0.6021         0         0
%>          0    0.5377         0         0         0
%>          0         0         0         0         0
%>          0         0         0         0         0
%>
%> >>
%> @endcode
% ==================================================
classdef sp < matfaust.proj.proj_gen
	properties
	end
	methods

		function proj = sp(shape, k, varargin)
			import matfaust.factparams.ConstraintInt
			% default values
			proj.constraint = ConstraintInt('sp', shape(1), shape(2), k, varargin{:});
		end
	end
end
