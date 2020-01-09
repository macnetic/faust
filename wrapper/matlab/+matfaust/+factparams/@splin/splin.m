classdef splin < matfaust.factparams.proj_gen
	properties
	end
	methods
		function proj = splin(shape, k, varargin)
			import matfaust.factparams.ConstraintInt
			proj.constraint = ConstraintInt('splin', shape(1), shape(2), k, varargin{:});
		end
	end
end
