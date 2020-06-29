%==================================================
%> @brief Functor for the SPCOL projector.
%>
%> A, the image matrix, is defined by \f$ \forall j \in \{1,\ldots,shape(2)\} \f$ the j-th column \f$  A_{*,j} \f$ is such that \f$ \|A_{*,j}\|_0 = k,  \| A\|_F = 1 (if normalized is true)\f$.
%==================================================
classdef spcol < matfaust.proj.proj_gen
	properties
	end
	methods
		function proj = spcol(shape, k, varargin)
			import matfaust.factparams.ConstraintInt
			proj.constraint = ConstraintInt('spcol', shape(1), shape(2), k, varargin{:});
		end
	end
end
