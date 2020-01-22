%==================================================
%> Functor that implements the NORMLIN projector. A, the image matrix, is defined by \f$ \forall i \in \{1,\ldots,shape(1)\} \| A_{i, *} \|_2 = s \| \f$.
%==================================================
classdef normlin < matfaust.factparams.proj_gen
	properties
	end
	methods
		function proj = normcol(shape, normval, varargin)
			import matfaust.factparams.ConstraintReal
			proj.constraint = ConstraintReal('normlin', shape(1), shape(2), normval, varargin{:});
		end
	end
end
