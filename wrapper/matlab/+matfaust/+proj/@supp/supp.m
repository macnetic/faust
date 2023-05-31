% ==================================================
%> @brief Functor for the the SUPP projector.
%>
%> A, the image matrix, is such that nonzeros(A) == nonzeros(S)
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b supp(S): returns a SUPP projector (functor), S defines the support matrix (entries are 0 or 1) for which the nonzeros entries define the values to keep from the input matrix into the ouput matrix of the projector.<br/>
%> &nbsp;&nbsp;&nbsp; @b supp(@b S,@b 'normalized', bool,@b 'pos', bool): the optional parameters are set. By default both normalized and pos are false.
%>
%> @param 'normalized', true: normalizes the projection image according to its Frobenius norm.
%> @param 'normalized', false: (the default) no normalization.
%> @param 'pos', true: skips the negative values (replaced by zero) of the input matrix.
%> @param 'pos', false: (the default) negative values are not skipped.
%> @retval proj the supp projector.
%>
%> @b Example
%> @code
%> >> import matfaust.proj.supp
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
%> >> S = zeros(5, 5);
%> >> S(M > .5) = 1; % the support of values > .5 in M
%> >> p = supp(S, 'normalized', false);
%> >> p(M)
%>
%> ans =
%>
%>          0         0         0         0    0.6119
%>     0.9507         0    0.9699         0         0
%>     0.7320    0.8662    0.8324    0.5248         0
%>     0.5987    0.6011         0         0         0
%>          0    0.7081         0         0         0
%> >>
%> @endcode
% ==================================================
classdef supp < matfaust.proj.proj_gen
	properties
	end
	methods

		function proj = supp(S, varargin)
			import matfaust.factparams.ConstraintMat
			proj.constraint = ConstraintMat('supp', S, varargin{:});
		end
	end
end
