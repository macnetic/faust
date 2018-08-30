classdef ParamsPalm4MSA < matfaust.ParamsFact
	properties (SetAccess = public)
		stop_crit
	end
	methods
		function p = ParamsPalm4MSA(varargin)
			import matfaust.ParamsFact
			MIN_NARGIN = 6;
			if(nargin < MIN_NARGIN)
				error(['matfaust.ParamsPalm4MSA() must receive at least ', int2str(MIN_NARGIN),' arguments.'])
			end
			num_facts = varargin{1};
			is_update_way_R2L = varargin{2};
			init_lambda = varargin{3};
			init_facts = varargin{4};
			constraints = varargin{5};
			stop_crit = varargin{6};
			% set default values
			step_size = ParamsFact.DEFAULT_STEP_SIZE;
			is_verbose = ParamsFact.DEFAULT_VERBOSITY;
			constant_step_size = ParamsFact.DEFAULT_CONSTANT_STEP_SIZE;
			if(nargin > MIN_NARGIN)
				step_size = varargin{MIN_NARGIN+1};
				if(nargin > MIN_NARGIN+1)
					constant_step_size = varargin{MIN_NARGIN+2};
					if(nargin > MIN_NARGIN+2)
						is_verbose = varargin{MIN_NARGIN+3};
					end
				end
			end
			% parent constructor handles verification for its own arguments
			p = p@matfaust.ParamsFact(num_facts, is_update_way_R2L, init_lambda, init_facts, ...
				constraints, step_size, constant_step_size, is_verbose);
			%p = p@matfaust.ParamsFact(varargin{:});
			if(~ isa(stop_crit, 'matfaust.StoppingCriterion'))
				error('matfaust.ParamsPalm4MSA argument (stop_crit) must be a matfaust.StoppingCriterion objects.')
			end
			p.stop_crit = stop_crit;
		end
	end
end
