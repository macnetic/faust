classdef ConstraintReal < matfaust.ConstraintGeneric
	properties (SetAccess = public)

	end
	methods
		function constraint = ConstraintReal(name, num_rows, num_cols, param)
			% check param is a real scalar
			if(~ isreal(param) || ~ isscalar(param))
				error('ConstraintReal must receive a real as param argument.')
			end
			constraint = constraint@matfaust.ConstraintGeneric(name, num_rows, num_cols, param);
			if(~ name.is_real_constraint())
				error('ConstraintReal first argument must be a ConstraintName with a real type name.')
			end
		end
	end
end
