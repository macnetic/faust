%% class ParamsHierarchicalFact
%%
classdef ParamsHierarchicalFact < matfaust.factparams.ParamsFact
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
			MIN_NARGIN = 6;
			if(nargin < MIN_NARGIN)
				error(['matfaust.factparams.ParamsHierarchicalFact() must receive at least',int2str(MIN_NARGIN),' arguments'])
			end
			num_facts = varargin{1};
			is_update_way_R2L = varargin{2};
			init_lambda = varargin{3};
			fact_constraints = varargin{4};
			res_constraints = varargin{5};
			constraints = {fact_constraints{:}, res_constraints{:}};
			% data_num_rows/data_num_cols are set by FaustFactory.fact_hierarchical()
			stop_crits = varargin{6};
			% set default values
			is_fact_side_left = matfaust.factparams.ParamsHierarchicalFact.DEFAULT_IS_FACT_SIDE_LEFT;
			step_size = matfaust.factparams.ParamsFact.DEFAULT_STEP_SIZE;
			is_verbose = matfaust.factparams.ParamsFact.DEFAULT_VERBOSITY;
			constant_step_size = matfaust.factparams.ParamsFact.DEFAULT_CONSTANT_STEP_SIZE;
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
			p = p@matfaust.factparams.ParamsFact(num_facts, is_update_way_R2L, init_lambda, ...
				constraints, step_size, constant_step_size, is_verbose);
			if(~ iscell(stop_crits))
				error('matfaust.factparams.ParamsHierarchicalFact 6th argument (stop_crits) must be a cell array.')
				if(length(stop_crits) ~= 2 )
					error('matfaust.factparams.ParamsHierarchicalFact 6th argument (stop_crits) must be a cell array of 2 elements.')
				end
				for i = 1:length(stop_crits)
					if(~ isa(stop_crits{i}, matfaust.factparams.StoppingCriterion))
						error('matfaust.factparams.ParamsHierarchicalFact 6th argument (stop_crits) must contain matfaust.factparams.StoppingCriterion objects.')
					end
				end
			end
			if(~ islogical(is_fact_side_left))
				error('matfaust.factparams.ParamsHierarchicalFact 11th argument (is_fact_side_left) must be logical.')
			end
			p.stop_crits = stop_crits;
			p.is_fact_side_left = is_fact_side_left;
			% auto-deduced to-factorize-matrix dim. sizes
			p.data_num_rows = p.constraints{1}.num_rows;
			p.data_num_cols = p.constraints{end}.num_cols;
		end



	end
end
