% ==================================================
%> @brief Functor for the TOEPLITZ projector.
%>
%> Each diagonal of the projector output matrix is filled with a constant value which is the mean of the elements in the same diagonal of the input matrix.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b toeplitz(shape): returns a TOEPLITZ projector (functor), shape defines the size of the input matrix (e.g. [1, 10]).<br/>
%> &nbsp;&nbsp;&nbsp; @b toeplitz(@b shape,@b 'normalized', bool,@b 'pos', bool): the optional parameters are set. By default both normalized and pos are false.
%>
%> @param 'normalized', true: normalizes the projection image according to its Frobenius norm.
%> @param 'normalized', false: (the default) no normalization.
%> @param 'pos', true: skips the negative values (replaced by zero) of the input matrix.
%> @param 'pos', false: (the default) negative values are not skipped.
%> @retval proj the toeplitz projector.
%>
%>
%> @b Example
%> @code
%> >> import matfaust.proj.toeplitz
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
%> >> p = toeplitz(size(M), 3);
%> >> p(M)
%>
%> ans =
%>
%>     0.1838    0.2153    0.0878    0.0689    0.2612
%>     0.2477    0.1838    0.2153    0.0878    0.0689
%>     0.2156    0.2477    0.1838    0.2153    0.0878
%>     0.2789    0.2156    0.2477    0.1838    0.2153
%>     0.0666    0.2789    0.2156    0.2477    0.1838
%>
%%> >>
%> @endcode
% ==================================================
classdef toeplitz < matfaust.proj.proj_gen
	properties
	end
	methods

		function proj = toeplitz(S, varargin)
			import matfaust.factparams.ConstraintMat
			proj.constraint = ConstraintMat('toeplitz', S, varargin{:});
		end
	end
end
