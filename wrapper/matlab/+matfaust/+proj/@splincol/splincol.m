classdef splincol < matfaust.factparams.proj_gen
	properties
	end
	methods
		function proj = splincol(shape, k, varargin)
			import matfaust.factparams.ConstraintInt
			% default values
			proj.constraint = ConstraintInt('splincol', shape(1), shape(2), k, varargin{:});
		end
	end
end
