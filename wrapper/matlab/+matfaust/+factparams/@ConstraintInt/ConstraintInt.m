% ======================================================================
%> This class represents an integer constraint on a matrix.
% ======================================================================
classdef ConstraintInt < matfaust.factparams.ConstraintGeneric
	properties (SetAccess = public)

	end
	methods
		function constraint = ConstraintInt(name, num_rows, num_cols, param, varargin)
			% check param is a int
			if(~ isreal(param) || ~ isscalar(param))
				error('ConstraintInt must receive an integer as param argument.')
			end

			constraint = constraint@matfaust.factparams.ConstraintGeneric(name, num_rows, num_cols, floor(param), varargin{:});
			if(~ isa(constraint.name, 'matfaust.factparams.ConstraintName') || ~ constraint.name.is_int_constraint())
				error(['ConstraintInt first argument must be a ConstraintName with a int type name ', ...
					'(name.is_int_constraint() must return True).'])
			end
		end
	end
end
