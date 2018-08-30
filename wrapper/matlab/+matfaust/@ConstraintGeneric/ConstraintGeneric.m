classdef ConstraintGeneric
	properties (SetAccess = public)
		name % obj of type ConstraintName
		num_rows
		num_cols
		param % type determined by child classes
	end
	methods
		function constraint = ConstraintGeneric(name, num_rows, num_cols, param)
			if(~ isa(name, 'matfaust.ConstraintName'))
				error('ConstraintGeneric first argument must be a matfaust.ConstraintName object.')
			end
			constraint.name = name;
			if(~ isreal(num_rows) || ~ isinteger(int64(num_rows)) || all(size(num_rows) ~= [1, 1]))
				error('ConstraintGeneric 2nd argument must be an integer.')
			end
			if(~ isreal(num_cols) || ~ isinteger(int64(num_cols)) || all(size(num_cols) ~= [1, 1]))
				error('ConstraintGeneric 3rd argument must be an integer.')
			end
			constraint.num_rows = floor(num_rows);
			constraint.num_cols = floor(num_cols);
			% child classes verify the param type
			constraint.param = param;
		end
	end
end
