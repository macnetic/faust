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
				nrows = param(end, 1);
				ncols = param(end, 2);
			elseif(isa(name, 'matfaust.factparams.ConstraintName') && name.name == matfaust.factparams.ConstraintName.ID || isstr(name) && matfaust.factparams.ConstraintName.str2name_int(name) == matfaust.factparams.ConstraintName.ID)
				nrows = param(1);
				ncols = param(2);
			else
				nrows = size(param, 1);
				ncols = size(param, 2);
			end
			constraint = constraint@matfaust.factparams.ConstraintGeneric(name, nrows, ncols, param, varargin{:});
			if(~ isa(constraint.name, 'matfaust.factparams.ConstraintName') || ~ constraint.name.is_mat_constraint())
				error('ConstraintMat first argument must be a ConstraintName with a matrix type name.')
			end

		end

		function [normalized, pos] = get_default_parameters(self)

			import matfaust.factparams.ConstraintName
			name = self.name.name;
			switch(name)
				%TODO: add ID
				case ConstraintName.ID
					normalized = false;
					pos = false;
				case {ConstraintName.TOEPLITZ,
					ConstraintName.CIRC,
					ConstraintName.HANKEL,
					ConstraintName.SUPP,
					ConstraintName.BLKDIAG}
					normalized = true;
					pos = false;
				case ConstraintName.CONST
					normalized = false;
					pos = false;
				otherwise
					error('Invalid ConstraintName')
			end
		end

	end
end
