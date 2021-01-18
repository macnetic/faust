%==================================================
%> @brief Functor for the HANKEL projector.
%>
%> Each antidiagonal of the projector output matrix is filled with a constant value which is the mean of the elements in the same antidiagonal of the input matrix.
%==================================================
classdef hankel < matfaust.proj.proj_gen
	properties
	end
	methods

		%===============================================
		%> @brief the Hankel projector.
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
		%> @retval the hankel projector.
		%===============================================
		function proj = hankel(shape, varargin)
			import matfaust.factparams.ConstraintMat
			proj.constraint = ConstraintMat('hankel', shape, varargin{:});
		end
	end
end
