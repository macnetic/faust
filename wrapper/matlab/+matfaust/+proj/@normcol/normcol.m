%==================================================
%> Functor that implements the NORMCOL projector. A, the image matrix, is defined by \f$ \forall j \in \{1,\ldots,shape(2)\} \| A_{*,j} \|_2 = s \| \f$.
%==================================================
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
