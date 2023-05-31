% ==================================================
%> @brief Functor for the SPCOL projector.
%>
%> A, the image matrix, is defined by \f$ \forall j \in \{1,\ldots,shape(2)\} \f$ the j-th column \f$  A_{*,j} \f$ is such that \f$ \|A_{*,j}\|_0 = k,  \| A\|_F = 1 \f$ (if normalized is true).
%>
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
%> @retval spcol projector.
%>
%> @b Example
%> @code
%> >> import matfaust.proj.spcol
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
%> >> p = spcol(size(M), 3, 'normalized', false);
%> >> p(M)
%>
%> ans =
%>
%>          0         0         0         0    0.6119
%>     0.9507         0    0.9699    0.3042         0
%>     0.7320    0.8662    0.8324    0.5248         0
%>     0.5987    0.6011    0.2123    0.4319    0.3664
%>          0    0.7081         0         0    0.4561
%>
%> >>
%> @endcode
% ==================================================
classdef spcol < matfaust.proj.proj_gen
	properties
	end
	methods

		function proj = spcol(shape, k, varargin)
			import matfaust.factparams.ConstraintInt
			proj.constraint = ConstraintInt('spcol', shape(1), shape(2), k, varargin{:});
		end
	end
end
