% ==================================================
%> @brief Functor for the SPLINCOL projector.
%>
%> It's the union of SPLIN and SPCOL projectors.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b splincol(shape, k): returns a SPLINCOL projector (functor), shape defines the size of the input matrix (e.g. [1, 10]), k defines the sparsity of the output matrix (k nnz coefficients per row or per column).<br/>
%> &nbsp;&nbsp;&nbsp; @b splincol(@b shape,@b k,@b 'normalized', bool,@b 'pos', bool): the optional parameters are set. By default both normalized and pos are false.
%>
%> @param shape: vector of size 2, to define the size of the input matrix.
%> @param k: the sparsity parameter.
%> @param 'normalized', true: normalizes the projection image according to its Frobenius norm.
%> @param 'normalized', false: (the default) no normalization.
%> @param 'pos', true: skips the negative values (replaced by zero) of the input matrix.
%> @param 'pos', false: (the default) negative values are not skipped.
%> @retval splincol projector.
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
%> >> p1 = splin(size(M), 3, 'normalized', false);
%> >> p2 = spcol(size(M), 3, 'normalized', false);
%> >> p = splincol(size(M), 3, 'normalized', false);
%> >> p1M = p1(M)
%>
%> p1M =
%>
%>     0.3745         0         0    0.1834    0.6119
%>     0.9507         0    0.9699    0.3042         0
%>     0.7320    0.8662    0.8324         0         0
%>     0.5987    0.6011         0    0.4319         0
%>          0    0.7081         0    0.2912    0.4561
%>
%> >> p2M = p2(M)
%>
%> p2M =
%>
%>          0         0         0         0    0.6119
%>     0.9507         0    0.9699    0.3042         0
%>     0.7320    0.8662    0.8324    0.5248         0
%>     0.5987    0.6011    0.2123    0.4319    0.3664
%>          0    0.7081         0         0    0.4561
%>
%> >> p(M)
%
%> ans =
%>
%>     0.3745         0         0    0.1834    0.6119
%>     0.9507         0    0.9699    0.3042         0
%>     0.7320    0.8662    0.8324    0.5248         0
%>     0.5987    0.6011    0.2123    0.4319    0.3664
%>          0    0.7081         0    0.2912    0.4561
%>
%> >> p1M(p1M == p2M) = 0;
%> >> p1M + p2M
%>
%> ans =
%>
%>     0.3745         0         0    0.1834    0.6119
%>     0.9507         0    0.9699    0.3042         0
%>     0.7320    0.8662    0.8324    0.5248         0
%>     0.5987    0.6011    0.2123    0.4319    0.3664
%>          0    0.7081         0    0.2912    0.4561
%>
%> >> all(all(p1M + p2M == p(M)))
%>
%> ans =
%>
%>  logical
%>
%>     1
%>
%> >>
%> @endcode
% ==================================================
classdef splincol < matfaust.proj.proj_gen
	properties
	end
	methods

		function proj = splincol(shape, k, varargin)
			import matfaust.factparams.ConstraintInt
			% default values
			proj.constraint = ConstraintInt('splincol', shape(1), shape(2), k, varargin{:});
		end
	end
end
