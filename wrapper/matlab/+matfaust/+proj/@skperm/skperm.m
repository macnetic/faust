% ==================================================
%> @brief Functor for the SKPERM projector.
%>
%> This projector returns the matrix with k nonzeros of the input matrix per column and per row that maximizes the Frobenius norm.
%>
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
%>
%>
%> @b Example
%> @code
%> >> import matfaust.proj.*
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
%> >> p = skperm(size(M), 3, 'normalized', false);
%> >> p(M)
%>
%> ans =
%>
%>          0         0    0.0206    0.1834    0.6119
%>     0.9507         0    0.9699    0.3042         0
%>     0.7320    0.8662    0.8324         0         0
%>     0.5987    0.6011         0         0    0.3664
%>          0    0.7081         0    0.2912    0.4561
%>
%> >>
%> @endcode
%>
% ==================================================
classdef skperm < matfaust.proj.proj_gen
	properties
	end
	methods

		function proj = skperm(shape, k, varargin)
			import matfaust.factparams.ConstraintInt
			% default values
			proj.constraint = ConstraintInt('skperm', shape(1), shape(2), k, varargin{:});
		end
	end
end
