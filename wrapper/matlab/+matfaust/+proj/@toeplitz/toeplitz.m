%==================================================
%> @brief Functor for the TOEPLITZ projector.
%>
%> Each diagonal of the projector output matrix is filled with a constant value which is the mean of the elements in the same diagonal of the input matrix.
%==================================================
classdef toeplitz < matfaust.proj.proj_gen
	properties
	end
	methods

		%===============================================
		%> @brief The toeplitz projector.
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
		%===============================================
		function proj = toeplitz(S, varargin)
			import matfaust.factparams.ConstraintMat
			proj.constraint = ConstraintMat('toeplitz', S, varargin{:});
		end
	end
end
