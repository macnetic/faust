%==================================================
%> Functor for the CONST projector.
%==================================================
classdef const < matfaust.proj.proj_gen
	properties
	end
	methods
		%===============================================
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b const(C): returns a CONST projector (functor), C defines the constant matrix to which the ouptut matrix of projector must be equal (before normalization).<br/>
		%> &nbsp;&nbsp;&nbsp; @b const(@b C,@b 'normalized', bool,@b 'pos', bool): the optional parameters are set. By default both normalized and pos are false.
		%>
		%> @param 'normalized', true: normalizes the projection image according to its Frobenius norm.
		%> @param 'normalized', false: (the default) no normalization.
		%> @param 'pos', true: skips the negative values (replaced by zero) of the input matrix.
		%> @param 'pos', false: (the default) negative values are not skipped.
		%===============================================
		function proj = const(C, varargin)
			import matfaust.factparams.ConstraintMat
			shape = size(C);
			proj.constraint = ConstraintMat('const', C, varargin{:});
		end
	end
end
