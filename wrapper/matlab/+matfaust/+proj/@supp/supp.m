%==================================================
%> @brief Functor that implements the SUPP projector.
%>
%> A, the image matrix, is such that nonzeros(A) == nonzeros(S)
%==================================================
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
