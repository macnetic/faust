classdef normcol < matfaust.factparams.proj_gen
	properties
	end
	methods
		function proj = normcol(shape, normval, varargin)
			import matfaust.factparams.ConstraintReal
			proj.constraint = ConstraintReal('normcol', shape(1), shape(2), normval, varargin{:});
		end
	end
end
