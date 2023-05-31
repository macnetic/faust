% ==================================================
%> @brief Functor for the anti-circulant projector.
%>
%> The output matrix of the projector is an anti-circulant matrix.
%>
%> Each constant used to fill a pair of diagonals of the output matrix is the mean of all the input matrix entries that it replaces in the output matrix.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b anticirc(shape): returns an anticirculant projector (functor), shape defines the size of the input matrix (e.g. [1, 10]).<br/>
%> &nbsp;&nbsp;&nbsp; @b anticirc(@b shape,@b 'normalized', bool,@b 'pos', bool): the optional parameters are set. By default both normalized and pos are false.
%>
%> @param 'normalized', true: normalizes the projection image according to its Frobenius norm.
%> @param 'normalized', false: (the default) no normalization.
%> @param 'pos', true: skips the negative values (replaced by zero) of the input matrix.
%> @param 'pos', false: (the default) negative values are not skipped.
%>
%> @retval the anticirc projector.
%>
%>
%> @b Example
%> @code
%> >> import matfaust.proj.anticirc
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
%> >> p = anticirc(size(M), 3);
%> >> p(M)
%>
%> ans =
%>
%>     0.1726    0.1773    0.1293    0.2708    0.2207
%>     0.1773    0.1293    0.2708    0.2207    0.1726
%>     0.1293    0.2708    0.2207    0.1726    0.1773
%>     0.2708    0.2207    0.1726    0.1773    0.1293
%>     0.2207    0.1726    0.1773    0.1293    0.2708
%>
%> >>
%> @endcode
% ==================================================
classdef anticirc < matfaust.proj.proj_gen
	properties
	end
	methods

		function proj = anticirc(shape, varargin)
			import matfaust.factparams.ConstraintMat
			proj.constraint = ConstraintMat('anticirc', shape, varargin{:});
		end
	end
end
