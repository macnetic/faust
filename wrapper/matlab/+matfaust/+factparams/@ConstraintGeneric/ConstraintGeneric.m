% ======================================================================
%> This is the parent class for representing a factor constraint in FAuST factorization algorithms.
% ======================================================================
classdef (Abstract) ConstraintGeneric
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

		function pM = project(this, M, varargin)
			% optional key-value arguments
			argc = length(varargin);
			normalized = true;
			if(argc > 0)
				for i=1:argc
					switch(varargin{i})
						case 'normalized'
							if(argc == i || ~ islogical(normalized))
								error('normalized keyword arg. is not followed by a boolean')
							end
							normalized = varargin{i+1};
						otherwise
							if(isstr(varargin{i}))
								error([ varargin{i} ' unrecognized argument'])
							end
					end
				end
			end
			if(isreal(M))
				pM = mexFaustReal('prox', M, this.name.name, this.param, normalized);
			else
				pM = mexFaustCplx('prox', M, this.name.name, this.param, normalized);
			end
		end
	end
end
