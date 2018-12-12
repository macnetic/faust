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
		IDX_IS_FACT_SIDE_LEFT = 1
		OPT_ARG_NAMES2 = { 'is_fact_side_left' }
	end
	methods
		function p = ParamsHierarchicalFact(fact_constraints, res_constraints, stop_crit1, stop_crit2, varargin)
			import matfaust.factparams.*
			if(~ iscell(fact_constraints))
				error('fact_constraints (argument 1) must be a cell array')
			end
			if(~ iscell(res_constraints))
				error('res_constraints (argument 2) must be a cell array')
			end
			if(length(fact_constraints) ~= length(res_constraints))
				error('lengths of fact_constraints and res_constraints must be equal.')
			end
			if(~ isa(stop_crit1, 'StoppingCriterion'))
				error('stop_crit1 (argument 3) must a StoppingCriterion')
			end
			if(~ isa(stop_crit2, 'StoppingCriterion'))
				error('stop_crit2 (argument 4) must a StoppingCriterion')
			end
			constraints = {fact_constraints{:}, res_constraints{:}};
			stop_crits = {stop_crit1, stop_crit2};
			% infer number of factors from constraints
			num_facts = length(fact_constraints)+1;
			parent_args = {};
			opt_arg_map = containers.Map();
			if(length(varargin) > 0)
				% retrieve all optional argument key-value pairs
				opt_arg_names = {ParamsFact.OPT_ARG_NAMES, ParamsHierarchicalFact.OPT_ARG_NAMES2};
				opt_arg_names = {opt_arg_names{1}{:}, opt_arg_names{2}{:}};
				ParamsFact.parse_opt_args(varargin, opt_arg_names, opt_arg_map)
				% gather all parent argument key-value pairs

				for i=1:length(ParamsFact.OPT_ARG_NAMES)
					if(opt_arg_map.isKey(ParamsFact.OPT_ARG_NAMES{i}))
						parent_args = [ parent_args, {ParamsFact.OPT_ARG_NAMES{i}}, {opt_arg_map(ParamsFact.OPT_ARG_NAMES{i}) }];
					end
				end
				% parent constructor handles verification for its own arguments
			end
			p = p@matfaust.factparams.ParamsFact(num_facts, constraints, parent_args{:});
			% data_num_rows/data_num_cols are set by FaustFactory.fact_hierarchical()
			% set default values
			is_fact_side_left = ParamsHierarchicalFact.DEFAULT_IS_FACT_SIDE_LEFT;
			if(opt_arg_map.isKey(ParamsHierarchicalFact.OPT_ARG_NAMES2{p.IDX_IS_FACT_SIDE_LEFT}))
				is_fact_side_left = opt_arg_map(ParamsHierarchicalFact.OPT_ARG_NAMES2{p.IDX_IS_FACT_SIDE_LEFT});
			end
			if(~ islogical(is_fact_side_left))
				error('matfaust.factparams.ParamsHierarchicalFact: is_fact_side_left argument must be logical.')
			end
			p.stop_crits = stop_crits;
			p.is_fact_side_left = is_fact_side_left;
			% auto-deduced to-factorize-matrix dim. sizes
			p.data_num_rows = p.constraints{1}.num_rows;
			p.data_num_cols = p.constraints{end}.num_cols;
		end



	end
end
