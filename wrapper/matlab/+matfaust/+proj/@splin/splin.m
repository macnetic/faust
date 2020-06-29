%==========================================
%> @brief Functor for the SPLIN projector.
%>
%> A, the image matrix, is defined by \f$ \forall i \in \{1,\ldots,shape(1)\} \f$ the i-th row \f$  A_{i,*} \f$ is such that \f$ \|A_{i,*}\|_0 = k,  \| A\|_F = 1 (if normalized is true)\f$.
%==========================================
classdef splin < matfaust.proj.proj_gen
	properties
	end
	methods
		function proj = splin(shape, k, varargin)
			import matfaust.factparams.ConstraintInt
			proj.constraint = ConstraintInt('splin', shape(1), shape(2), k, varargin{:});
		end
	end
end
