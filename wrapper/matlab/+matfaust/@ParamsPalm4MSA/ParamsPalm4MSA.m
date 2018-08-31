classdef ParamsPalm4MSA < matfaust.ParamsFact
	properties (SetAccess = public)
		init_facts
		stop_crit
	end
	methods
		function p = ParamsPalm4MSA(varargin)
			import matfaust.ParamsFact
			MIN_NARGIN = 5;
			if(nargin < MIN_NARGIN)
				error(['matfaust.ParamsPalm4MSA() must receive at least ', int2str(MIN_NARGIN),' arguments.'])
			end
			num_facts = varargin{1};
			is_update_way_R2L = varargin{2};
			init_lambda = varargin{3};
			constraints = varargin{4};
			stop_crit = varargin{5};
			% set default values
			step_size = ParamsFact.DEFAULT_STEP_SIZE;
			is_verbose = ParamsFact.DEFAULT_VERBOSITY;
			constant_step_size = ParamsFact.DEFAULT_CONSTANT_STEP_SIZE;
			is_init_facts_to_default = nargin <= MIN_NARGIN;
			if(~ is_init_facts_to_default)
				init_facts = varargin{MIN_NARGIN+1};
				if(nargin > MIN_NARGIN+1)
					step_size = varargin{MIN_NARGIN+2};
					if(nargin > MIN_NARGIN+2)
						constant_step_size = varargin{MIN_NARGIN+3};
						if(nargin > MIN_NARGIN+3)
							is_verbose = varargin{MIN_NARGIN+4};
						end
					end
				end
			end

			% parent constructor handles verification for its own arguments
			p = p@matfaust.ParamsFact(num_facts, is_update_way_R2L, init_lambda, ...
				constraints, step_size, constant_step_size, is_verbose);
			if(is_init_facts_to_default || iscell(init_facts) && length(init_facts) == 0)
				init_facts = cell(num_facts, 1)
				if(is_update_way_R2L)
					zeros_id = num_facts
				else
					zeros_id = 1
				end
				for i=1:num_facts
					if(i ~= zeros_id)
						init_facts{i} = eye(constraints{i}.num_rows, constraints{i}.num_cols)
					end
				end
				init_facts{zeros_id} = ...
					zeros(constraints{zeros_id}.num_rows, constraints{zeros_id}.num_cols)
			elseif(~ iscell(init_facts)) % TODO: check init_facts length
				error('matfaust.ParamsFactPalm4MSA 4th argument (init_facts) must be a cell array.')
			else
				for i = 1:length(init_facts)
					if(~ ismatrix(init_facts{i}) || ~ isnumeric(init_facts{i}))
						error('matfaust.ParamsFactPalm4MSA 4th argument (init_facts) must contain matrices.')
						%TODO: check matrix dimensions
					end
				end
			end
			p.init_facts = init_facts
			if(~ isa(stop_crit, 'matfaust.StoppingCriterion'))
				error('matfaust.ParamsPalm4MSA argument (stop_crit) must be a matfaust.StoppingCriterion objects.')
			end
			p.stop_crit = stop_crit;
		end
	end
end
