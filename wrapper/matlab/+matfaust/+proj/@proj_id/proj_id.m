%==================================================
%> Functor for the ID projector.
%==================================================
classdef proj_id < matfaust.proj.proj_gen
	properties
	end
	methods
		%===============================================
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
		%===============================================
		function proj = proj_id(mat_size, varargin)
			import matfaust.factparams.ConstraintMat
			proj.constraint = ConstraintMat('id', mat_size, varargin{:});
		end
	end
end
