%% class ConstraintGeneric
%%

classdef ConstraintGeneric
	properties (SetAccess = public)
		name % obj of type ConstraintName
		num_rows
		num_cols
		param % type determined by child classes
	end
	methods
		function constraint = ConstraintGeneric(name, num_rows, num_cols, param)
			%ENOTE: ischar(name{:})
			if(ischar(name) || iscell(name) && ischar(name{:}))
				name = matfaust.factparams.ConstraintName(name);
			end
			constraint.name = name;
			if(~ isreal(num_rows) || ~ isscalar(num_rows))
				error('ConstraintGeneric 2nd argument must be an integer.')
			end
			if(~ isreal(num_cols) || ~ isscalar(num_cols))
				error('ConstraintGeneric 3rd argument must be an integer.')
			end
			constraint.num_rows = floor(num_rows);
			constraint.num_cols = floor(num_cols);
			% child classes verify the param type
			constraint.param = param;
		end
	end
end
