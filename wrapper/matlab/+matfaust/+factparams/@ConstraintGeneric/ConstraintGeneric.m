% ======================================================================
%> This is the parent class for representing a factor constraint in FAuST factorization algorithms.
% ======================================================================
classdef (Abstract) ConstraintGeneric
	properties (SetAccess = public)
		name % obj of type ConstraintName
		num_rows
		num_cols
		param % type determined by child classes
		normalized % flag to normalize or not the projection image afterward
		default_normalized % flag to determine if the normalized flag has its default value or user customed one
		pos % flag to set negative matrix elements to zero through the projection
	end
	methods
		function constraint = ConstraintGeneric(name, num_rows, num_cols, param, varargin)
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
			% optional key-value arguments
			argc = length(varargin);
			constraint.normalized = true;
			constraint.default_normalized = true;
			constraint.pos = false;
			if(argc > 0)
				for i=1:argc
					switch(varargin{i})
						case 'normalized'
							if(argc == i || ~ islogical(varargin{i+1}))
								error('normalized keyword arg. is not followed by a boolean')
							end
							constraint.normalized = varargin{i+1};
							constraint.default_normalized = false;
						case 'pos'
							if(argc == i || ~ islogical(varargin{i+1}))
								error('pos keyword arg. is not followed by a boolean')
							end
							constraint.pos = varargin{i+1};
						otherwise
							if(isstr(varargin{i}))
								error([ varargin{i} ' unrecognized argument'])
							end
					end
				end
			end
		end

		function pM = project(this, M)
			if(isreal(M))
				pM = mexFaustReal('prox', M, this.name.name, this.param, this.normalized, this.pos);
			else
				pM = mexFaustCplx('prox', M, this.name.name, this.param, this.normalized, this.pos);
			end
		end
	end
end
