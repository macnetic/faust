
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
		% =========================================================
		%>	@brief Constructor.
		%>
		%>	@param fact_constraints a ConstraintList object or a list of matfaust.proj.proj_gen objects to define the constraints of the main factor at each level of the factorization hierarchy (the first one for the first factorization and so on).
		%>	@param res_constraints a ConstraintList object or a list of matfaust.proj.proj_gen objects to define the constraints to apply to the residual factor at each level of the factorization hierarchy (the first one for the first factorization and so on).
		%>	@param stop_crit1 	a matfaust.factparams.StoppingCriterion instance which defines the algorithm stopping criterion for the local optimization of the 2 terms of the last factorization (a main factor and a residual).
		%>	@param stop_crit2 a matfaust.factparams.StoppingCriterion instance which defines the algorithm stopping criterion for the global optimization.
		%>	@param 'is_update_way_R2L', bool 	if true matfaust.fact.palm4msa (called for each optimization stage) will update factors from the right to the left, otherwise it's done in reverse order.
		%>	@param 'init_lambda', real the scale scalar initial value for the global optimization (by default the value is one). It applies only to local optimization at each iteration (the global optimization lambda is updated consequently). 
		%>	@param 'step_size', real the initial step of the PALM descent for both local and global optimization stages. 
		%>	@param 'constant_step_size', bool if true the step_size keeps constant along the algorithm iterations otherwise it is updated before every factor update. 
		%>	@param 'is_fact_side_left', bool if true the leftmost factor is factorized, otherwise it's the rightmost. 
		%>	@param 'is_verbose', bool true to enable the verbose mode. 
		%>	@param 'factor_format', str (optional) 'dynamic' (by default), 'dense', or 'sparse'. If 'dense' or 'sparse' then all factors will be respectively full arrays or sparse matrices. If 'dynamic' is used then the algorithm determines the format of each factor automatically in order to decrease the memory footprint of the Faust. This option is available only on the 2020 backend matfaust.fact.palm4msa, matfaust.fact.hierarchical or matfaust.fact.palm4msa_mhtp, matfaust.fact.hierarchical_mhtp.
		%>	@param 'packing_RL', bool true (by default) to pre-compute R and L products (only available with 2020 backend of pyfaust.fact.hierarchical).
		%>	@param 'no_normalization', bool false (by default), if true it disables the normalization of prox output matrix in PALM4MSA algorithm. Note that this option is experimental (only available with 2020 backend of pyfaust.fact.palm4msa).
		%>  @param 'no_lambda', bool: false (by default), if true it disables the lambda scalar factor in the PALM4MSA algorithm which consists basically to set it always to one (it also lowers the algorithm cost).
		%>	@param 'norm2_max_iter', int maximum number of iterations of power iteration algorithm (default to 100). Used for computing 2-norm.
		%>	@param 'norm2_threshold', real power iteration algorithm threshold (default to 1e-6). Used for computing 2-norm.
		%>	@param 'grad_calc_opt_mode', int the mode used for computing the PALM gradient. It can be one value among matfaust.factparams.ParamsFact.EXTERNAL_OPT, matfaust.factparams.ParamsFact.INTERNAL_OPT or matfaust.factparams.ParamsFact.DISABLED_OPT. This parameter is experimental, its value shouln't be changed. 
		% =========================================================
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
				% because data_num_rows and data_num_cols have been inferred according to the constraints and is_fact_side_left
				s = size(M);
				bool = all(s == [this.data_num_rows, this.data_num_cols]);
			end
		end

		% =========================================================
		%> @brief This function handles the conversion of the structure instance to the mex format used by wrapper implementation.
		%>
		%> @note this function is not intended to be used by users.
		%>
		%> @retval The mex structure instance as expected by the mex wrapper used by matfaust.fact.hierarchical.
		% =========================================================
		function mex_params = to_mex_struct(this, M)
			mex_constraints = cell(2, this.num_facts-1);
			M_is_single = strcmp('single', class(M));
			is_single = @(mat) strcmp('single', class(mat));
			%mex_fact_constraints = cell(1, this.num_facts-1)
			for i=1:this.num_facts-1
				cur_cell = cell(1, 4);
				cur_cell{1} = this.constraints{i}.name.conv2str();
				cur_cell{2} = this.constraints{i}.param;
				cur_cell{3} = this.constraints{i}.num_rows;
				cur_cell{4} = this.constraints{i}.num_cols;
				%mex_fact_constraints{i} = cur_cell;
				mex_constraints{1,i} = cur_cell;
			end
			%mex_residuum_constraints = cell(1, this.num_facts-1)
			for i=1:this.num_facts-1
				cur_cell = cell(1, 4);
				cur_cell{1} = this.constraints{i+this.num_facts-1}.name.conv2str();
				cur_cell{2} = this.constraints{i+this.num_facts-1}.param;
				cur_cell{3} = this.constraints{i+this.num_facts-1}.num_rows;
				cur_cell{4} = this.constraints{i+this.num_facts-1}.num_cols;
				%mex_residuum_constraints{i} = cur_cell;
				mex_constraints{2,i} = cur_cell;
			end
			% verify ConstraintMat type consistency with the factorized matrix
			for i=0:1
				if(i == 0)
					cons_factor = 'main factor';
				else
					cons_factor = 'residual factor';
				end
				for j=1:this.num_facts-1
					ci = (this.num_facts-1)*i+j;
					if(isa(this.constraints{ci}, 'matfaust.factparams.ConstraintMat'))
						if(~ M_is_single == is_single(this.constraints{ci}.param))
							error(['ParamsHierarchical.constraints, ' cons_factor ' constraint #' int2str(j) ' must be in single precision (resp. double) if the matrix to factorize is in single precision (resp. double).'])
						end
						if(~ isreal(M) == isreal(this.constraints{ci}.param))
							error(['ParamsHierarchical.constraints, ' cons_factor ' constraint #' int2str(j) ' must be real (resp. complex) if the matrix to factorize is real (resp. complex).'])
						end
					end

				end
			end
			mex_params = struct('nfacts', this.num_facts, 'cons', {mex_constraints}, 'niter1', this.stop_crits{1}.num_its,'niter2', this.stop_crits{2}.num_its, 'sc_is_criterion_error', this.stop_crits{1}.is_criterion_error, 'sc_error_treshold', this.stop_crits{1}.tol, 'sc_max_num_its', this.stop_crits{1}.maxiter, 'sc_is_criterion_error2', this.stop_crits{2}.is_criterion_error, 'sc_error_treshold2', this.stop_crits{2}.tol, 'sc_max_num_its2', this.stop_crits{2}.maxiter, 'nrow', this.data_num_rows, 'ncol', this.data_num_cols, 'fact_side', this.is_fact_side_left, 'update_way', this.is_update_way_R2L, 'verbose', this.is_verbose, 'init_lambda', this.init_lambda, 'factor_format', matfaust.factparams.ParamsFact.factor_format_str2int(this.factor_format), 'packing_RL', this.packing_RL, 'no_normalization', this.no_normalization, 'no_lambda', this.no_lambda, 'norm2_threshold', this.norm2_threshold, 'norm2_max_iter', this.norm2_max_iter, 'step_size', this.step_size, 'constant_step_size', this.constant_step_size);
			if(~ (islogical(this.use_MHTP) &&  this.use_MHTP == false))
				% use_MHTP must be a MHTPParams if not false (cf. ParamsFact)
				if(~ isa(this.use_MHTP, 'matfaust.factparams.MHTPParams'))
					error('use_MHTP is not a MHTPParams')
				end
				mhtp_p = this.use_MHTP;
				mex_params.mhtp_num_its = mhtp_p.num_its;
				mex_params.mhtp_constant_step_size = mhtp_p.constant_step_size;
				mex_params.mhtp_step_size = mhtp_p.step_size;
				mex_params.mhtp_palm4msa_period = mhtp_p.palm4msa_period;
				mex_params.mhtp_updating_lambda = mhtp_p.updating_lambda;
			end
		end

	end
end
