
% =========================================================
%> The parent class to set input parameters for the hierarchical factorization algorithm.
% =========================================================
classdef ParamsHierarchical < matfaust.factparams.ParamsFact
	properties (SetAccess = public)
		stop_crits
		data_num_rows
		data_num_cols
		is_fact_side_left
	end
	properties(Constant, SetAccess = protected, Hidden)
		DEFAULT_IS_FACT_SIDE_LEFT = false
		IDX_IS_FACT_SIDE_LEFT = 1
		OPT_ARG_NAMES2 = { 'is_fact_side_left' }
	end
	methods
		function p = ParamsHierarchical(fact_constraints, res_constraints, stop_crit1, stop_crit2, varargin)
			import matfaust.factparams.*
			if(iscell(fact_constraints))
				for i=1:length(fact_constraints)
					if(isa(fact_constraints{i}, 'matfaust.proj.proj_gen'))
						fact_constraints{i} = fact_constraints{i}.constraint;
					end
				end
			end
			if(iscell(res_constraints))
				for i=1:length(res_constraints)
					if(isa(res_constraints{i}, 'matfaust.proj.proj_gen'))
						res_constraints{i} = res_constraints{i}.constraint;
					end
				end
			end
			if(isa(fact_constraints, 'ConstraintList'))
				fact_constraints = fact_constraints.clist;
			end
			if(isa(res_constraints, 'ConstraintList'))
				res_constraints = res_constraints.clist;
			end
			if(~ iscell(fact_constraints))
				error('fact_constraints (argument 1) must be a cell array of matfaust.factparams.ConstraintGeneric or matfaust.proj.proj_gen.')
			end
			if(~ iscell(res_constraints))
				error('res_constraints (argument 2) must be a cell array of matfaust.factparams.ConstraintGeneric or matfaust.proj.proj_gen.')
			end
			if(length(fact_constraints) ~= length(res_constraints))
				error('lengths of fact_constraints and res_constraints must be equal.')
			end
			if(~ isa(stop_crit1, 'StoppingCriterion'))
				error('stop_crit1 (argument 3) must be a StoppingCriterion')
			end
			if(~ isa(stop_crit2, 'StoppingCriterion'))
				error('stop_crit2 (argument 4) must a StoppingCriterion')
			end

			parent_args = {};
			opt_arg_map = containers.Map();
			if(length(varargin) > 0)
				% retrieve all optional argument key-value pairs
				opt_arg_names = {ParamsFact.OPT_ARG_NAMES, ParamsHierarchical.OPT_ARG_NAMES2};
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
			% set default values
			is_fact_side_left = ParamsHierarchical.DEFAULT_IS_FACT_SIDE_LEFT;
			if(opt_arg_map.isKey(ParamsHierarchical.OPT_ARG_NAMES2{ParamsHierarchical.IDX_IS_FACT_SIDE_LEFT}))
				is_fact_side_left = opt_arg_map(ParamsHierarchical.OPT_ARG_NAMES2{ParamsHierarchical.IDX_IS_FACT_SIDE_LEFT});
			end
			if(~ islogical(is_fact_side_left))
				error('matfaust.factparams.ParamsHierarchical: is_fact_side_left argument must be logical.')
			end
			if(is_fact_side_left)
				constraints = {res_constraints{:}, fact_constraints{:}};
			else
				constraints = {fact_constraints{:}, res_constraints{:}};
			end
			stop_crits = {stop_crit1, stop_crit2};
			% infer number of factors from constraints
			num_facts = length(fact_constraints)+1;
			p = p@matfaust.factparams.ParamsFact(num_facts, constraints, parent_args{:});
			p.stop_crits = stop_crits;
			p.is_fact_side_left = is_fact_side_left;
			% auto-deduced to-factorize-matrix dim. sizes
			if(p.is_fact_side_left)
				p.data_num_rows = res_constraints{end}.num_rows;
				p.data_num_cols = fact_constraints{1}.num_cols;
			else
				p.data_num_rows = p.constraints{1}.num_rows;
				p.data_num_cols = p.constraints{end}.num_cols;
			end
		end

		function bool = is_mat_consistent(this, M)
			if(~ ismatrix(M))
				error('M must be a matrix.')
			else
				% no need to take care of is_fact_side_left
				% because data_num_rows and data_num_cols have been infered according to the constraints and is_fact_side_left
				s = size(M);
				bool = all(s == [this.data_num_rows, this.data_num_cols]);
			end
		end

	end
end
