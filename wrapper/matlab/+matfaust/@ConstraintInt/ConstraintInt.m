%% class ConstraintInt
%%

classdef ConstraintInt < matfaust.ConstraintGeneric
	properties (SetAccess = public)

	end
	methods
		function constraint = ConstraintInt(name, num_rows, num_cols, param)
			% check param is a int
			if(~ isreal(param) || ~ isscalar(param))
				error('ConstraintInt must receive an integer as param argument.')
			end
			if(~ isa(name, 'matfaust.ConstraintName') || ~ name.is_int_constraint())
				error(['ConstraintInt first argument must be a ConstraintName with a int type name ', ...
					'(name.is_int_constraint() must return True).'])
			end
			constraint = constraint@matfaust.ConstraintGeneric(name, num_rows, num_cols, floor(param));
		end
	end
end
