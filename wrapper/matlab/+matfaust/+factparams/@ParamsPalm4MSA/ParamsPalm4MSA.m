% =========================================================
%> The class is to set input parameters for the Palm4MSA algorithm.
% =========================================================
classdef ParamsPalm4MSA < matfaust.factparams.ParamsFact
	properties (SetAccess = public)
		init_facts
		stop_crit
	end
	properties (Constant, SetAccess = protected, Hidden)
		IDX_INIT_FACTS = 1
		OPT_ARG_NAMES2 = { 'init_facts' }
	end
	methods
		% =========================================================
		%>	@brief Constructor.
		%>
		%>	@param constraints a ConstraintList object or a cell array of matfaust.proj.proj_gen objects to define the constraints of the factors.
		%>	@param stop_crit a matfaust.factparams.StoppingCriterion instance which defines the algorithm stopping criterion.
		%>	@param 'init_facts', cell array of matrices (optional) if defined, matfaust.fact.palm4msa will initialize the factors with init_facts (by default the first factor to be updated is initialized to zero and the others to identity. Note that the so called first factor can be the rightmost or the leftmost depending on the is_update_way_R2L argument).
		%>	@param 'is_update_way_R2L', bool (optional) if true matfaust.fact.palm4msa will update factors from the right to the left, otherwise it's done in reverse order.
		%>	@param 'init_lambda', real (optional) the scale scalar initial value (by default the value is one).
		%>	@param 'step_size', real (optional) the initial step of the PALM descent.
		%>	@param 'constant_step_size', real if true the step_size keeps constant along the algorithm iterations otherwise it is updated before every factor update.
		%>	@param 'is_verbose', bool (optional) True to enable the verbose mode.
		%>	parameter is experimental, its value shouldn't be changed.
		%>	@param 'norm2_max_iter', real (optional) maximum number of iterations of power iteration algorithm. Used for computing 2-norm.
		%>	@param 'norm2_threshold', real (optional) power iteration algorithm threshold (default to 1e-6). Used for computing 2-norm.
		%>	@param 'grad_calc_opt_mode', int (optional) the mode used for computing the PALM gradient. It can be one value among matfaust.factparams.ParamsFact.EXTERNAL_OPT, matfaust.factparams.ParamsFact.INTERNAL_OPT or matfaust.factparams.ParamsFact.DISABLED_OPT. This
		%>	@param 'factor_format', str (optional) 'dynamic' (by default), 'dense', or 'sparse'. If 'dense' or 'sparse' then all factors will be respectively full arrays or sparse matrices. If 'dynamic' is used then the algorithm determines the format of each factor automatically in order to decrease the memory footprint of the Faust. This option is available only on the 2020 backend matfaust.fact.palm4msa, matfaust.fact.hierarchical or matfaust.fact.palm4msa_mhtp, matfaust.fact.hierarchical_mhtp.
		%>	@param 'packing_RL', bool true (by default) to pre-compute R and L products (only available with 2020 backend of pyfaust.fact.palm4msa).
		%>	@param 'no_normalization', bool False (by default), if True it disables the normalization of prox output matrix in PALM4MSA algorithm. Note that this option is experimental (only available with 2020 backend of pyfaust.fact.palm4msa).
		%>
		%>
		% =========================================================
		function p = ParamsPalm4MSA(constraints, stop_crit, varargin)
			import matfaust.factparams.*
			if(isa(constraints, 'ConstraintList'))
				constraints = constraints.clist;
			end
			if(~ iscell(constraints))
				error('constraints (argument 1) must be a cell array')
			end
			num_facts = length(constraints);
			if(~ isa(stop_crit, 'StoppingCriterion'))
				error('stop_crit (argument 2) must a StoppingCriterion')
			end
			parent_args = {};
			opt_arg_map = containers.Map();
			if(length(varargin) > 0)
				% retrieve all optional argument key-value pairs
				opt_arg_names = {ParamsFact.OPT_ARG_NAMES, ParamsPalm4MSA.OPT_ARG_NAMES2};
				opt_arg_names = {opt_arg_names{1}{:}, opt_arg_names{2}{:}};
				ParamsFact.parse_opt_args(varargin, opt_arg_names, opt_arg_map)
				% gather all parent argument key-value pairs
				for i=1:length(ParamsFact.OPT_ARG_NAMES)
					if(opt_arg_map.isKey(ParamsFact.OPT_ARG_NAMES{i}))
						parent_args = [ parent_args, {ParamsFact.OPT_ARG_NAMES{i}}, {opt_arg_map(ParamsFact.OPT_ARG_NAMES{i}) }];
					end
				end%
				% parent constructor handles verification for its own arguments
			end
			p = p@matfaust.factparams.ParamsFact(num_facts, constraints, parent_args{:});
			init_facts_name = p.OPT_ARG_NAMES{p.IDX_INIT_FACTS};
			try
				init_facts = opt_arg_map(init_facts_name);
			catch
				% arg int_facts not passed
			end
			if(~ exist(init_facts_name) || iscell(init_facts) && length(init_facts) == 0)
				init_facts = p.get_default_init_facts(num_facts);
			elseif(~ iscell(init_facts)) % TODO: check init_facts length
				error(['matfaust.factparams.ParamsFactPalm4MSA argument ' init_facts_name ' must be a cell array.'])
			else
				% check init_facts
				% TODO: check the number of init_facts
				for i = 1:length(init_facts)
					if(~ ismatrix(init_facts{i}) || ~ isnumeric(init_facts{i}))
						error(['matfaust.factparams.ParamsFactPalm4MSA ' init_facts_name ' argument must contain matrices.'])
						% matrix dimensions are tested later by fact_palm4msa() with is_mat_consistent()
						% TODO: add to is_mat_consistent() the checking of init_facts
					end
				end
			end
			p.init_facts = init_facts;
			if(~ isa(stop_crit, 'matfaust.factparams.StoppingCriterion'))
				error('matfaust.factparams.ParamsPalm4MSA argument (stop_crit) must be a matfaust.factparams.StoppingCriterion objects.')
			end
			p.stop_crit = stop_crit;
		end

		function mex_params = to_mex_struct(this, M)
			mex_constraints = cell(1, length(this.constraints));
			for i=1:length(this.constraints)
				cur_cell = cell(1, 4);
				cur_cell{1} = this.constraints{i}.name.conv2str();
				cur_cell{2} = this.constraints{i}.param;
				cur_cell{3} = this.constraints{i}.num_rows;
				cur_cell{4} = this.constraints{i}.num_cols;
				mex_constraints{i} = cur_cell;
			end
			% put mex_constraints in a cell array again because mex eats one level of array
			mex_params = struct('data', M, 'nfacts', this.num_facts, 'cons', {mex_constraints}, 'init_facts', {this.init_facts}, 'niter', this.stop_crit.num_its, 'sc_is_criterion_error', this.stop_crit.is_criterion_error, 'sc_error_treshold', this.stop_crit.tol, 'sc_max_num_its', this.stop_crit.maxiter, 'update_way', this.is_update_way_R2L, 'grad_calc_opt_mode', this.grad_calc_opt_mode, 'constant_step_size', this.constant_step_size, 'step_size', this.step_size, 'verbose', this.is_verbose, 'norm2_max_iter', this.norm2_max_iter, 'norm2_threshold', this.norm2_threshold, 'init_lambda', this.init_lambda, 'factor_format', matfaust.factparams.ParamsFact.factor_format_str2int(this.factor_format), 'packing_RL', this.packing_RL, 'no_normalization', this.no_normalization);
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
	methods
		function init_facts = get_default_init_facts(p, num_facts)
			init_facts = cell(num_facts, 1);
			if(p.is_update_way_R2L)
				zeros_id = num_facts;
			else
				zeros_id = 1;
			end
			for i=1:num_facts
				if(i ~= zeros_id)
					init_facts{i} = eye(p.constraints{i}.num_rows, p.constraints{i}.num_cols);
				end
			end
			init_facts{zeros_id} = ...
				zeros(p.constraints{zeros_id}.num_rows, p.constraints{zeros_id}.num_cols);
		end
	end
end
