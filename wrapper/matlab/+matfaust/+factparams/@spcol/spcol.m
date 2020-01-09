classdef spcol < matfaust.factparams.proj_gen
	properties
	end
	methods
		function proj = spcol(shape, k, varargin)
			import matfaust.factparams.ConstraintInt
			proj.constraint = ConstraintInt('spcol', shape(1), shape(2), k, varargin{:});
		end
	end
end
