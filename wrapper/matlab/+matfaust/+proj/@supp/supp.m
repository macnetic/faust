%==================================================
%> @brief Functor for the the SUPP projector.
%>
%> A, the image matrix, is such that nonzeros(A) == nonzeros(S)
%==================================================
classdef supp < matfaust.proj.proj_gen
	properties
	end
	methods

		%===============================================
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b supp(S): returns a SUPP projector (functor), S defines the support matrix (entries are 0 or 1) for which the nonzeros entries define the values to keep from the input matrix into the ouput matrix of the projector.<br/>
		%> &nbsp;&nbsp;&nbsp; @b supp(@b S,@b 'normalized', bool,@b 'pos', bool): the optional parameters are set. By default both normalized and pos are false.
		%>
		%> @param 'normalized', true: normalizes the projection image according to its Frobenius norm.
		%> @param 'normalized', false: (the default) no normalization.
		%> @param 'pos', true: skips the negative values (replaced by zero) of the input matrix.
		%> @param 'pos', false: (the default) negative values are not skipped.
		%> @retval the supp projector.
		%===============================================
		function proj = supp(S, varargin)
			import matfaust.factparams.ConstraintMat
			proj.constraint = ConstraintMat('supp', S, varargin{:});
		end
	end
end
