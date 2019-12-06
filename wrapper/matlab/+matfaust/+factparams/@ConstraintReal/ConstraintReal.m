% =========================================================
%> This class represents a real constraint on a matrix.
% =========================================================
classdef ConstraintReal < matfaust.factparams.ConstraintGeneric
	properties (SetAccess = public)

	end
	methods
		function constraint = ConstraintReal(name, num_rows, num_cols, param, varargin)
			% check param is a real scalar
			if(~ isreal(param) || ~ isscalar(param))
				error('ConstraintReal must receive a real as param argument.')
			end
			constraint = constraint@matfaust.factparams.ConstraintGeneric(name, num_rows, num_cols, param, varargin{:});
			if(constraint.default_normalized && constraint.name.name == matfaust.factparams.ConstraintName.CONST)
				constraint.normalized = false;
			end
			if(~ isa(constraint.name, 'matfaust.factparams.ConstraintName') || ~ constraint.name.is_real_constraint())
				error('ConstraintReal first argument must be a ConstraintName with a real type name.')
			end
		end
	end
end
