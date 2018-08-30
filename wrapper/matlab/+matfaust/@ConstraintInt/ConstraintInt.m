classdef ConstraintInt < matfaust.ConstraintGeneric
	properties (SetAccess = public)

	end
	methods
		function constraint = ConstraintInt(name, num_rows, num_cols, param)
			% check param is a int
			if(~ isreal(param) || ~ isinteger(int64(param)) || all(size(param) ~= [1, 1]))
				error('ConstraintInt must receive an integer as param argument.')
			end
			constraint = constraint@matfaust.ConstraintGeneric(name, num_rows, num_cols, floor(param));
			if(~ name.is_int_constraint())
				error('ConstraintInt first argument must be a ConstraintName with a int type name.')
			end

		end
	end
end
