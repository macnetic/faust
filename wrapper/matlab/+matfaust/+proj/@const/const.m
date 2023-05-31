% ==================================================
%> @brief Functor for the CONST projector.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b const(C): returns a CONST projector (functor), C defines the constant matrix to which the ouptut matrix of projector must be equal (before normalization).<br/>
%> &nbsp;&nbsp;&nbsp; @b const(@b C,@b 'normalized', bool,@b 'pos', bool): the optional parameters are set. By default both normalized and pos are false.
%>
%> @param 'normalized', true: normalizes the projection image according to its Frobenius norm.
%> @param 'normalized', false: (the default) no normalization.
%> @param 'pos', true: skips the negative values (replaced by zero) of the input matrix.
%> @param 'pos', false: (the default) negative values are not skipped.
%> @retval p the const projector. <br/>
%>
%> @b Example
%> @code
%> >> import matfaust.proj.const
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
%> >> C = rand(5, 5)
%>
%> C =
%>
%>     0.7852    0.6075    0.8084    0.1220    0.6625
%>     0.1997    0.1705    0.3046    0.4952    0.3117
%>     0.5142    0.0651    0.0977    0.0344    0.5201
%>     0.5924    0.9489    0.6842    0.9093    0.5467
%>     0.0465    0.9656    0.4402    0.2588    0.1849
%>
%> >> p = const(C);
%> >> p(M)
%>
%> ans =
%>
%>     0.7852    0.6075    0.8084    0.1220    0.6625
%>     0.1997    0.1705    0.3046    0.4952    0.3117
%>     0.5142    0.0651    0.0977    0.0344    0.5201
%>     0.5924    0.9489    0.6842    0.9093    0.5467
%>     0.0465    0.9656    0.4402    0.2588    0.1849
%>
%> >>
%> @endcode
%>
% ==================================================
classdef const < matfaust.proj.proj_gen
	properties
	end
	methods
		function proj = const(C, varargin)
			import matfaust.factparams.ConstraintMat
			shape = size(C);
			proj.constraint = ConstraintMat('const', C, varargin{:});
		end
	end
end
