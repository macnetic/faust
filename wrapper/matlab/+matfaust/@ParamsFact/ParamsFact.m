classdef (Abstract) ParamsFact
	properties (SetAccess = public)
		num_facts
		is_update_way_R2L
		init_lambda
		init_facts
		step_size
		constant_step_size
		constraints
		is_verbose
	end
	properties (Constant)
		DEFAULT_STEP_SIZE = 10^-16
		DEFAULT_VERBOSITY = false
		DEFAULT_CONSTANT_STEP_SIZE = false
	end
	methods
		function p = ParamsFact(varargin)
			MIN_NARGIN = 5;
			if(nargin < MIN_NARGIN)
				error(['matfaust.ParamsFact() must receive at least ', int2str(MIN_NARGIN),' arguments.'])
			end
			num_facts = varargin{1};
			is_update_way_R2L = varargin{2};
			init_lambda = varargin{3};
			init_facts = varargin{4};
			constraints = varargin{5};
			% set default values
			step_size = matfaust.ParamsFact.DEFAULT_STEP_SIZE;
			is_verbose = matfaust.ParamsFact.DEFAULT_VERBOSITY;
			constant_step_size = matfaust.ParamsFact.DEFAULT_CONSTANT_STEP_SIZE;
			if(~ isscalar(num_facts) || ~ isinteger(int64(num_facts)))
				error('matfaust.ParamsFact 1st argument (num_facts) must be an integer.')
			end
			if(~ islogical(is_update_way_R2L))
				error('matfaust.ParamsFact 2nd argument (is_update_way_R2L) must be logical.')
			end
			if(~ isscalar(init_lambda))
				error('matfaust.ParamsFact 3rd argument (init_lambda) must be a scalar.')
			end
			if(~ iscell(init_facts)) % TODO: check init_facts length
				error('matfaust.ParamsFact 4th argument (init_facts) must be a cell array.')
			end
			for i = 1:length(init_facts)
				if(~ ismatrix(init_facts{i}) || ~ isnumeric(init_facts{i}))
					error('matfaust.ParamsFact 4th argument (init_facts) must contain matrices.')
					%TODO: check matrix dimensions
				end
			end
			if(~ iscell(constraints))
				error('matfaust.ParamsFact 5th argument (constraints) must be a cell array.')
			end
			for i = 1:length(constraints) %TODO: check constraints length in sub-class
				if(~ isa(constraints{i}, 'matfaust.ConstraintGeneric'))
					error('matfaust.ParamsFact 5th argument (constraints) must contain matfaust.ConstraintGeneric objects.')
				end
			end
			if(nargin > MIN_NARGIN)
				step_size = varargin{MIN_NARGIN+1};
				if(nargin > MIN_NARGIN+1)
					constant_step_size = varargin{MIN_NARGIN+2};
					if(nargin > MIN_NARGIN+2)
						is_verbose = varargin{MIN_NARGIN+3};
					end
				end
			end
			if(~ isscalar(step_size))
				error(['matfaust.ParamsHierarchicalFact ', int2str(MIN_NARGIN+1), 'th argument (step_size) must be an integer.'])
			end
			if(~ islogical(constant_step_size))
				error(['matfaust.ParamsFact ', int2str(MIN_NARGIN+2), 'th argument (constant_step_size) must be logical.'])
			end
			if(~ islogical(is_verbose))
				error(['matfaust.ParamsFact ',int2str(MIN_NARGIN+3),' argument (is_verbose) must be logical.'])
			end
			p.num_facts = num_facts;
			p.is_update_way_R2L = is_update_way_R2L;
			p.init_lambda = init_lambda;
			p.init_facts = init_facts;
			p.constraints = constraints;
			p.step_size = step_size;
			p.is_verbose = is_verbose;
			p.constant_step_size = constant_step_size;
		end
	end
end
