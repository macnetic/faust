% ==================================================
%> @brief Functor for the HANKEL projector.
%>
%> Each antidiagonal of the projector output matrix is filled with a constant value which is the mean of the elements in the same antidiagonal of the input matrix.
%>
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b hankel(shape): returns a HANKEL projector (functor), shape defines the size of the input matrix (e.g. [1, 10]).<br/>
%> &nbsp;&nbsp;&nbsp; @b hankel(@b shape,@b 'normalized', bool,@b 'pos', bool): the optional parameters are set. By default both normalized and pos are false.
%>
%> @param 'normalized', true: normalizes the projection image according to its Frobenius norm.
%> @param 'normalized', false: (the default) no normalization.
%> @param 'pos', true: skips the negative values (replaced by zero) of the input matrix.
%> @param 'pos', false: (the default) negative values are not skipped.
%> @retval hankel projector.
%>
%>
%> @b Example
%> @code
%> >> import matfaust.proj.hankel
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
%> >> p = hankel(size(M), 3);
%> >> p(M)
%>
%> ans =
%>
%>     0.1632    0.2411    0.1177    0.2852    0.2184
%>     0.2411    0.1177    0.2852    0.2184    0.1726
%>     0.1177    0.2852    0.2184    0.1726    0.1316
%>     0.2852    0.2184    0.1726    0.1316    0.1433
%>     0.2184    0.1726    0.1316    0.1433    0.1987
%>
%> >>
%> @endcode
% ==================================================
classdef hankel < matfaust.proj.proj_gen
	properties
	end
	methods

		function proj = hankel(shape, varargin)
			import matfaust.factparams.ConstraintMat
			proj.constraint = ConstraintMat('hankel', shape, varargin{:});
		end
	end
end
