% ======================================================================
%> This class represents a matrix-based constraint on a matrix.
% ======================================================================
classdef ConstraintMat < matfaust.factparams.ConstraintGeneric
	properties (SetAccess = public)

	end
	methods
		function constraint = ConstraintMat(name, param, varargin)
			% check param is a mat
			if(~ ismatrix(param) || ~ isnumeric(param) && ~ islogical(param))
				error('ConstraintMat must receive a matrix as param argument.')
			end
			if(islogical(param))
				param = real(param);
			end
			if(isa(name, 'matfaust.factparams.ConstraintName') && name.name == matfaust.factparams.ConstraintName.BLKDIAG || isstr(name) && matfaust.factparams.ConstraintName.str2name_int(name) == matfaust.factparams.ConstraintName.BLKDIAG)
				nrows = param(end,1);
				ncols = param(end,2);
			else
				nrows = size(param, 1);
				ncols = size(param, 2);
			end
			constraint = constraint@matfaust.factparams.ConstraintGeneric(name, nrows, ncols, param, varargin{:});
			if(~ isa(constraint.name, 'matfaust.factparams.ConstraintName') || ~ constraint.name.is_mat_constraint())
				error('ConstraintMat first argument must be a ConstraintName with a matrix type name.')
			end
			if(constraint.default_normalized && constraint.name.name == matfaust.factparams.ConstraintName.CONST)
				% for CONST proj the default is to not normalize
				constraint.normalized = false;
				% for SUPP it is true (which is handled by ConstraintGeneric parent
			end
		end
	end
end
