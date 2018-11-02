%% class ConstraintMat
%%

classdef ConstraintMat < matfaust.ConstraintGeneric
	properties (SetAccess = public)

	end
	methods
		function constraint = ConstraintMat(name, param)
			% check param is a mat
			if(~ ismatrix(param) || ~ isnumeric(param))
				error('ConstraintMat must receive a matrix as param argument.')
			end
			constraint = constraint@matfaust.ConstraintGeneric(name, size(param, 1), size(param, 2), param);
			if(~ isa(name, 'matfaust.ConstraintName') || ~ name.is_mat_constraint())
				error('ConstraintMat first argument must be a ConstraintName with a matrix type name.')
			end
		end
	end
end
