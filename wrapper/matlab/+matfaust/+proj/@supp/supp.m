classdef supp < matfaust.factparams.proj_gen
	properties
	end
	methods
		function proj = supp(S, varargin)
			import matfaust.factparams.ConstraintMat
			proj.constraint = ConstraintMat('supp', S, varargin{:});
		end
	end
end
