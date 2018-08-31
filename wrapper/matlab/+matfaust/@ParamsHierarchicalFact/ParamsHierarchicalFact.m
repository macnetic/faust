classdef ParamsHierarchicalFact < matfaust.ParamsFact
	properties (SetAccess = public)
		stop_crits
		data_num_rows
		data_num_cols
		is_fact_side_left
	end
	properties(Constant)
		DEFAULT_IS_FACT_SIDE_LEFT = false
	end
	methods
		function p = ParamsHierarchicalFact(varargin)
			MIN_NARGIN = 7;
			if(nargin < MIN_NARGIN)
				error(['matfaust.ParamsHierarchicalFact() must receive at least',int2str(MIN_NARGIN),' arguments'])
			end
			num_facts = varargin{1};
			is_update_way_R2L = varargin{2};
			init_lambda = varargin{3};
			constraints = varargin{4};
			data_num_rows = floor(varargin{5});
			data_num_cols = floor(varargin{6});
			stop_crits = varargin{7};
			% set default values
			is_fact_side_left = matfaust.ParamsHierarchicalFact.DEFAULT_IS_FACT_SIDE_LEFT;
			step_size = matfaust.ParamsFact.DEFAULT_STEP_SIZE;
			is_verbose = matfaust.ParamsFact.DEFAULT_VERBOSITY;
			constant_step_size = matfaust.ParamsFact.DEFAULT_CONSTANT_STEP_SIZE;
			if(nargin > MIN_NARGIN+1)
				step_size = varargin{MIN_NARGIN+2};
				if(nargin > MIN_NARGIN+2)
					constant_step_size = varargin{MIN_NARGIN+3};
					if(nargin > MIN_NARGIN+3)
						is_verbose = varargin{MIN_NARGIN+4};
						if(nargin > MIN_NARGIN+4)
							is_fact_side_left = varargin{MIN_NARGIN+5};
						end
					end
				end
			end
			% parent constructor handles verification for its own arguments
			p = p@matfaust.ParamsFact(num_facts, is_update_way_R2L, init_lambda, ...
				constraints, step_size, constant_step_size, is_verbose);
			if(~ isscalar(data_num_rows) || ~ isreal(data_num_rows))
				error('matfaust.ParamsHierarchicalFact 5th argument (data_num_rows) must be an integer.')
			else
				data_num_rows = floor(data_num_rows);
			end
			if(~ isscalar(data_num_cols) || ~ isreal(data_num_cols))
				error('matfaust.ParamsHierarchicalFact 6th argument (data_num_cols) must be an integer.')
			else
				data_num_cols = floor(data_num_cols);
			end
			if(~ iscell(stop_crits))
				error('matfaust.ParamsHierarchicalFact 7th argument (stop_crits) must be a cell array.')
				if(length(stop_crits) ~= 2 )
					error('matfaust.ParamsHierarchicalFact 7th argument (stop_crits) must be a cell array of 2 elements.')
				end
				for i = 1:length(stop_crits)
					if(~ isa(stop_crits{i}, matfaust.StoppingCriterion))
						error('matfaust.ParamsHierarchicalFact 7th argument (stop_crits) must contain matfaust.StoppingCriterion objects.')
					end
				end
			end
			if(~ islogical(is_fact_side_left))
				error('matfaust.ParamsHierarchicalFact 8th argument (is_fact_side_left) must be logical.')
			end
			p.stop_crits = stop_crits;
			p.is_fact_side_left = is_fact_side_left;
			p.data_num_rows = data_num_rows;
			p.data_num_cols = data_num_cols;
			p.is_fact_side_left = is_fact_side_left;
		end
	end
end
