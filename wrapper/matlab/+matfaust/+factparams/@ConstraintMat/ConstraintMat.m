% ======================================================================
%> This class represents a matrix-based constraint on a matrix.
% ======================================================================
classdef ConstraintMat < matfaust.factparams.ConstraintGeneric
	properties (SetAccess = public)

	end
	methods
		function constraint = ConstraintMat(name, param)
			% check param is a mat
			if(~ ismatrix(param) || ~ isnumeric(param))
				error('ConstraintMat must receive a matrix as param argument.')
			end
			constraint = constraint@matfaust.factparams.ConstraintGeneric(name, size(param, 1), size(param, 2), param);
			if(~ isa(constraint.name, 'matfaust.factparams.ConstraintName') || ~ constraint.name.is_mat_constraint())
				error('ConstraintMat first argument must be a ConstraintName with a matrix type name.')
			end
		end
	end
end
