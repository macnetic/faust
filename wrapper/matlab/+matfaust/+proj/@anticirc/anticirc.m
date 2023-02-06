%==================================================
%> @brief Functor for the anti-circulant projector.
%>
%> The output matrix of the projector is an anti-circulant matrix.
%>
%> Each constant used to fill a pair of diagonals of the output matrix is the mean of all the input matrix entries that it replaces in the output matrix.
%==================================================
classdef anticirc < matfaust.proj.proj_gen
	properties
	end
	methods

		%===============================================
		%> @brief The anti-circulant projector.
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
		%===============================================
		function proj = anticirc(shape, varargin)
			import matfaust.factparams.ConstraintMat
			proj.constraint = ConstraintMat('anticirc', shape, varargin{:});
		end
	end
end
