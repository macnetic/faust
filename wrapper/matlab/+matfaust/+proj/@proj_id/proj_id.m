% ==================================================
%> @brief Functor for the ID projector.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b proj_id([M, N]): returns a ID projector (functor). This identity projector returns the input matrix whatever it is. Optionally a normalization or a positivity filter can be applied.<br/>
%> &nbsp;&nbsp;&nbsp; @b proj_id(@b [M, N],@b 'normalized', bool,@b 'pos', bool): the optional parameters are set. By default both normalized and pos are false.
%>
%> @param mat_size: the size of the input matrix (eg. [M, N]).
%> @param 'normalized', true: normalizes the projection image according to its Frobenius norm.
%> @param 'normalized', false: (the default) no normalization.
%> @param 'pos', true: skips the negative values (replaced by zero) of the input matrix.
%> @param 'pos', false: (the default) negative values are not skipped.
%> @retval the proj_id projector.
%>
%>
%> @b Example
%> @code
%> >> import matfaust.proj.proj_id
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
%> >> p = proj_id(size(M), 'normalized', false);
%> >> p(M)
%>
%> ans =
%>
%>     0.3745    0.1560    0.0206    0.1834    0.6119
%>     0.9507    0.0581    0.9699    0.3042    0.1395
%>     0.7320    0.8662    0.8324    0.5248    0.2921
%>     0.5987    0.6011    0.2123    0.4319    0.3664
%>     0.1560    0.7081    0.1818    0.2912    0.4561
%>
%> >>
%> @endcode
% ==================================================
classdef proj_id < matfaust.proj.proj_gen
	properties
	end
	methods
		function proj = proj_id(mat_size, varargin)
			import matfaust.factparams.ConstraintMat
			proj.constraint = ConstraintMat('id', mat_size, varargin{:});
		end
	end
end
