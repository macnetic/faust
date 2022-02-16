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
			if(~ isa(constraint.name, 'matfaust.factparams.ConstraintName') || ~ constraint.name.is_real_constraint())
				error('ConstraintReal first argument must be a ConstraintName with a real type name.')
			end
		end

		function [normalized, pos] = get_default_parameters(self)
			% all Real constraints have the same default parameters
			normalized = false;
			pos = false;
		end

	end

end
