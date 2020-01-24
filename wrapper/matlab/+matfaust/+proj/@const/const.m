%==================================================
%> Functor for the CONST projector.
%==================================================
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
