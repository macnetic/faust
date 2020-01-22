%==================================================
%> @brief Functor that implements the SPLINCOL projector.
%>
%> It's the union of SPLIN and SPCOL projectors.
%==================================================
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
