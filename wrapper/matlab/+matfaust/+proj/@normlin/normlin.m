%==================================================
%> @brief Functor for the NORMLIN projector.
%>
%> A, the image matrix, is defined by \f$ \forall i \in \{1,\ldots,shape(1)\} \f$ the i-th row \f$  A_{i,*} \f$ is such that \f$ \| A_{i, *} \|_2 = s \f$.
%==================================================
classdef normlin < matfaust.proj.proj_gen
	properties
	end
	methods
		function proj = normlin(shape, normval, varargin)
			import matfaust.factparams.ConstraintReal
			proj.constraint = ConstraintReal('normlin', shape(1), shape(2), normval, varargin{:});
		end
	end
end
