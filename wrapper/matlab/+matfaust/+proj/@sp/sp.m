%==================================================
%> @brief Functor that implements the SP projector. A, the projected matrix, is such that \f$ \| A \|_0 = k,  \| A\|_F = 1\f$.
%==================================================
classdef sp < matfaust.factparams.proj_gen
	properties
	end
	methods
		function proj = sp(shape, k, varargin)
			import matfaust.factparams.ConstraintInt
			% default values
			proj.constraint = ConstraintInt('sp', shape(1), shape(2), k, varargin{:});
		end
	end
end
