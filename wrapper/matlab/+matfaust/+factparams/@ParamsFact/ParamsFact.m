%% class ParamsFact
%%

classdef (Abstract) ParamsFact
	properties (SetAccess = public)
		num_facts
		is_update_way_R2L
		init_lambda
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
			MIN_NARGIN = 4;
			if(nargin < MIN_NARGIN)
				error(['matfaust.factparams.ParamsFact() must receive at least ', int2str(MIN_NARGIN),' arguments.'])
			end
			num_facts = varargin{1};
			is_update_way_R2L = varargin{2};
			init_lambda = varargin{3};
			constraints = varargin{4};
			% set default values
			step_size = matfaust.factparams.ParamsFact.DEFAULT_STEP_SIZE;
			is_verbose = matfaust.factparams.ParamsFact.DEFAULT_VERBOSITY;
			constant_step_size = matfaust.factparams.ParamsFact.DEFAULT_CONSTANT_STEP_SIZE;
			if(~ isscalar(num_facts) || ~ isreal(num_facts))
				error('matfaust.factparams.ParamsFact 1st argument (num_facts) must be an integer.')
			else
				num_facts = floor(num_facts);
			end
			if(~ islogical(is_update_way_R2L))
				error('matfaust.factparams.ParamsFact 2nd argument (is_update_way_R2L) must be logical.')
			end
			if(~ isscalar(init_lambda))
				error('matfaust.factparams.ParamsFact 3rd argument (init_lambda) must be a scalar.')
			end

			if(~ iscell(constraints))
				error('matfaust.factparams.ParamsFact 4th argument (constraints) must be a cell array.')
			end
			for i = 1:length(constraints) %TODO: check constraints length in sub-class
				if(~ isa(constraints{i}, 'matfaust.factparams.ConstraintGeneric'))
					error('matfaust.factparams.ParamsFact 5th argument (constraints) must contain matfaust.factparams.ConstraintGeneric objects.')
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
				step_size
				error(['matfaust.factparams.ParamsHierarchicalFact ', int2str(MIN_NARGIN+1), 'th argument (step_size) must be a real.'])
			end
			if(~ islogical(constant_step_size))
				error(['matfaust.factparams.ParamsFact ', int2str(MIN_NARGIN+2), 'th argument (constant_step_size) must be logical.'])
			end
			if(~ islogical(is_verbose))
				error(['matfaust.factparams.ParamsFact ',int2str(MIN_NARGIN+3),' argument (is_verbose) must be logical.'])
			end
			p.num_facts = num_facts;
			p.is_update_way_R2L = is_update_way_R2L;
			p.init_lambda = init_lambda;
			p.constraints = constraints;
			p.step_size = step_size;
			p.is_verbose = is_verbose;
			p.constant_step_size = constant_step_size;
		end
	end
end
