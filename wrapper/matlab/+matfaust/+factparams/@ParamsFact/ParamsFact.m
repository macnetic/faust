
% =========================================================
%> The parent abstract class to represent the general factorization parameters.
% =========================================================
classdef (Abstract) ParamsFact
	properties (SetAccess = public)
		num_facts
		is_update_way_R2L
		init_lambda
		step_size
		constant_step_size
		constraints
		is_verbose
		grad_calc_opt_mode
		factor_format
		packing_RL
		norm2_max_iter
		norm2_threshold
		use_MHTP
		no_normalization
	end
	properties (Constant, SetAccess = protected, Hidden)
		DEFAULT_STEP_SIZE = 10^-16
		DEFAULT_VERBOSITY = false
		DEFAULT_CONSTANT_STEP_SIZE = false
		DEFAULT_INIT_LAMBDA = 1.0
		DEFAULT_IS_UPDATE_WAY_R2L = false
		DEFAULT_PACKING_RL = true
		DEFAULT_FACTOR_FORMAT = 'dynamic'
		DEFAULT_NO_NORMALIZATION = false
		%> @note 0 means the default value defined in C++ core (100).
		DEFAULT_NORM2_MAX_ITER = 0 % norm2 parameters to 0 in order to use default parameters from C++ core
		%> @note 0 means the default value defined in C++ core (1e-6).
		DEFAULT_NORM2_THRESHOLD = 0
		% constructor opt. arguments indices
		IDX_IS_UPDATE_WAY_R2L = 1
		IDX_INIT_LAMBDA = 2
		IDX_STEP_SIZE = 3
		IDX_CONSTANT_STEP_SIZE = 4
		IDX_VERBOSITY = 5
		IDX_GRAD_CALC_OPT_MODE = 6
		IDX_NORM2_MAX_ITER = 7
		IDX_NORM2_THRESHOLD = 8
		IDX_FACTOR_FORMAT = 9
		IDX_PACKING_RL = 10
		IDX_NO_NORMALIZATION = 11
		% flags to control the optimization of the multiplication L'(LSR)R' in PALM4MSA
		%> gradient product optimization is disabled
		DISABLED_OPT = 0
		%> gradient product optimization is handled internally
		INTERNAL_OPT = 1
		%> gradient product optimization is handled externally
		EXTERNAL_OPT = 2
		%> default is EXTERNAL_OPT
		DEFAULT_OPT = 2
		% the order of names matters and must respect the indices above
		OPT_ARG_NAMES = {'is_update_way_R2L', 'init_lambda', 'step_size', 'constant_step_size', 'is_verbose', 'grad_calc_opt_mode', 'norm2_max_iter', 'norm2_threshold', 'factor_format', 'packing_RL', 'no_normalization'}
	end
	methods
		function p = ParamsFact(num_facts, constraints, varargin)
			import matfaust.factparams.ParamsFact
			if(iscell(constraints))
				for i=1:length(constraints)
					if(isa(constraints{i}, 'matfaust.proj.proj_gen'))
						constraints{i} = constraints{i}.constraint;
					end
				end
			end
			% set default values
			is_update_way_R2L = ParamsFact.DEFAULT_IS_UPDATE_WAY_R2L;
			init_lambda = ParamsFact.DEFAULT_INIT_LAMBDA;
			step_size = ParamsFact.DEFAULT_STEP_SIZE;
			is_verbose = ParamsFact.DEFAULT_VERBOSITY;
			constant_step_size = ParamsFact.DEFAULT_CONSTANT_STEP_SIZE;
			grad_calc_opt_mode = ParamsFact.DEFAULT_OPT;
			norm2_threshold = ParamsFact.DEFAULT_NORM2_THRESHOLD;
			norm2_max_iter = ParamsFact.DEFAULT_NORM2_MAX_ITER;
			factor_format = ParamsFact.DEFAULT_FACTOR_FORMAT;
			packing_RL = ParamsFact.DEFAULT_PACKING_RL;
			no_normalization = ParamsFact.DEFAULT_NO_NORMALIZATION;
			% check mandatory arguments
			if(~ isscalar(num_facts) || ~ isreal(num_facts))
				error('matfaust.factparams.ParamsFact num_facts argument must be an integer.')
			else
				num_facts = floor(num_facts);
			end
			if(~ iscell(constraints))
				error(['matfaust.factparams.ParamsFact constraints argument must be a cell array.'])
			end
			for i = 1:length(constraints)
				if(isa(constraints{i}, 'matfaust.proj.proj_gen'))
					constraints{i} = constraints{i}.constraint;
				end
			end
			for i = 1:length(constraints) %ParamsFact.TODO: check constraints length in sub-class
				if(~ isa(constraints{i}, 'matfaust.factparams.ConstraintGeneric'))
					error(['matfaust.factparams.ParamsFact constraints argument must contain matfaust.factparams.ConstraintGeneric objects.'])
				end
			end
			% try to get optional arguments
			% construct a map of optional arguments from varargin
			% NOTE: it avoids a if imbrication of depth the number length(varargin)
			% NOTE: we can change the order of arguments in properties (ParamsFact.IDX), the code will stay valid
			opt_arg_map = containers.Map();
			ParamsFact.parse_opt_args(varargin, ParamsFact.OPT_ARG_NAMES, opt_arg_map)
			% now get values of opt args really passed
			% TODO: a function to hide the complexity returning the value arg or the default value if not set in map
			if(opt_arg_map.isKey(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_IS_UPDATE_WAY_R2L}))
				is_update_way_R2L = opt_arg_map(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_IS_UPDATE_WAY_R2L});
			end
			if(opt_arg_map.isKey(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_INIT_LAMBDA}))
				init_lambda = opt_arg_map(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_INIT_LAMBDA});
			end
			if(opt_arg_map.isKey(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_STEP_SIZE}))
				step_size = opt_arg_map(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_STEP_SIZE});
			end
			if(opt_arg_map.isKey(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_CONSTANT_STEP_SIZE}))
				constant_step_size = opt_arg_map(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_CONSTANT_STEP_SIZE});
			end
			if(opt_arg_map.isKey(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_VERBOSITY}))
				is_verbose = opt_arg_map(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_VERBOSITY});
			end
			if(opt_arg_map.isKey(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_GRAD_CALC_OPT_MODE}))
				grad_calc_opt_mode = opt_arg_map(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_GRAD_CALC_OPT_MODE});
			end
			if(opt_arg_map.isKey(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_NORM2_MAX_ITER}))
				norm2_max_iter = opt_arg_map(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_NORM2_MAX_ITER});
			end
			if(opt_arg_map.isKey(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_NORM2_THRESHOLD}))
				norm2_threshold = opt_arg_map(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_NORM2_THRESHOLD});
			end
			if(opt_arg_map.isKey(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_PACKING_RL}))
				packing_RL = opt_arg_map(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_PACKING_RL});
			end
			if(opt_arg_map.isKey(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_FACTOR_FORMAT}))
				factor_format = opt_arg_map(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_FACTOR_FORMAT});
			end
			if(opt_arg_map.isKey(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_NO_NORMALIZATION}))
				no_normalization = opt_arg_map(ParamsFact.OPT_ARG_NAMES{ParamsFact.IDX_NO_NORMALIZATION});
			end
			% then check validity of opt args (it's useless for default values but it's not too costfull)
			% TODO: group argument verifs by type in loops
			if(~ islogical(is_update_way_R2L))
				error(['matfaust.factparams.ParamsFact ', p.OPT_ARG_NAMES{ParamsFact.IDX_IS_UPDATE_WAY_R2L} ,' argument (is_update_way_R2L) must be logical.'])
			end
			if(~ isscalar(init_lambda))
				error(['matfaust.factparams.ParamsFact ', p.OPT_ARG_NAMES{ParamsFact.IDX_INIT_LAMBDA},' argument (init_lambda) must be a scalar.'])
			end
			if(~ isscalar(step_size))
				step_size
				error(['matfaust.factparams.ParamsFact ', p.OPT_ARG_NAMES{ParamsFact.IDX_STEP_SIZE}, ' argument (step_size) must be a real.'])
			end
			if(~ islogical(constant_step_size))
				error(['matfaust.factparams.ParamsFact ', p.OPT_ARG_NAMES{ParamsFact.IDX_CONSTANT_STEP_SIZE}, ' argument (constant_step_size) must be logical.'])
			end
			if(~ islogical(is_verbose))
				error(['matfaust.factparams.ParamsFact ', p.OPT_ARG_NAMES{ParamsFact.IDX_VERBOSITY},' argument (is_verbose) must be logical.'])
			end
			if(~ isscalar(grad_calc_opt_mode) && (grad_calc_opt_mode == ParamsFact.INTERNAL_OPT || grad_calc_opt_mode == ParamsFact.EXTERNAL_OPT || grad_calc_opt_mode == ParamsFact.DISABLED_OPT))
				error(['matfaust.factparams.ParamsFact ', p.OPT_ARG_NAMES{ParamsFact.IDX_GRAD_CALC_OPT_MODE},' argument (grad_calc_opt_mode) must be an integer equal to ParamsFact.INTERNAL_OPT, ParamsFact.EXTERNAL_OPT, or ParamsFact.DISABLED_OPT.'])
			end
			if(~ ischar(factor_format) && ~ isstring(factor_format) || ~ any(strcmp(factor_format, {'dense', 'sparse', 'dynamic'})))
				error(['matfaust.factparams.ParamsFact ', p.OPT_ARG_NAMES{ParamsFact.IDX_FACTOR_FORMAT},' argument must be a str in ''dense'', ''sparse'' or ''dynamic''.'])
			end
			if(~ islogical(packing_RL))
				error(['matfaust.factparams.ParamsFact ', p.OPT_ARG_NAMES{ParamsFact.IDX_PACKING_RL},' argument must be logical.'])
			end
			if(~ islogical(no_normalization))
				error(['matfaust.factparams.ParamsFact ', p.OPT_ARG_NAMES{ParamsFact.IDX_NO_NORMALIZATION},' argument must be logical.'])
			end
			if(~ isscalar(norm2_max_iter) || floor(norm2_max_iter) < norm2_max_iter)
				norm2_max_iter
				error(['matfaust.factparams.ParamsFact ', p.OPT_ARG_NAMES{ParamsFact.IDX_NORM2_MAX_ITER}, ' argument must be a int.'])
			end
 			if(~ isscalar(norm2_threshold))
				norm2_threshold
				error(['matfaust.factparams.ParamsFact ', p.OPT_ARG_NAMES{ParamsFact.IDX_NORM2_THRESHOLD}, ' argument must be a real.'])
			end
			p.num_facts = num_facts;
			p.is_update_way_R2L = is_update_way_R2L;
			p.init_lambda = init_lambda;
			p.constraints = constraints;
			p.step_size = step_size;
			p.is_verbose = is_verbose;
			p.constant_step_size = constant_step_size;
			p.grad_calc_opt_mode = grad_calc_opt_mode;
			p.factor_format = ParamsFact.factor_format_str2int(factor_format);
			p.packing_RL = packing_RL;
			p.no_normalization = no_normalization;
			p.norm2_max_iter = norm2_max_iter;
			p.norm2_threshold = norm2_threshold;
			p.use_MHTP = false; % by default no MHTP in PALM4MSA, neither in hierarchical fact.
		end

		function bool = is_mat_consistent(this, M)
			if(~ ismatrix(M))
				error('M must be a matrix.')
			else
				s = size(M);
				bool = s(1) == this.constraints{1}.num_rows && s(2) == this.constraints{end}.num_cols;
			end
		end

	end
	methods(Static)
		function ff_int = factor_format_str2int(factor_format)
			if(ischar(factor_format) || isstring(factor_format))
				if(strcmp('dense', factor_format))
					ff_int = 0;
				elseif(strcmp('sparse', factor_format))
					ff_int = 1;
				elseif(strcmp('dynamic', factor_format))
					ff_int = 2;
				end
			elseif(isreal(factor_format) && factor_format >= 0 && factor_format <=2)
				% already a int or real to truncate
				ff_int = floor(factor_format);
			else
				error('factor_format should be a char array: ''dense'', ''sparse'' or ''dynamic''')
			end
		end
		function parse_opt_args(cell_args, opt_arg_names, opt_arg_map)
			if(~ iscell(cell_args))
				error('cell_args must be a cell array')
			end
			nopt_args = length(cell_args);
			if(nopt_args > 0)
				i=1;
				while(i<=nopt_args)
					i_arg_matched = false;
					if(~ isstr(cell_args{i}))
						unrecognized_argument = cell_args{i};
						error('above value is not a keyword argument. It must be a char array.')
					end
					for j=1:length(opt_arg_names)
						if(strcmp(cell_args{i},opt_arg_names{j}))
							i_arg_matched = true;
							if(length(cell_args) >= i)
								opt_arg_map(opt_arg_names{j}) = cell_args{i+1};
							else
								error([ 'the keyword argument ' opt_arg_names{j} ' must be followed by a value argument'])
							end
						end
					end
					if(~ i_arg_matched)
						unrecognized_argument = cell_args{i};
						error(['Value ' unrecognized_argument ' is not a recognized keyword argument.' ])
					end
					i = i + 2;
				end
			end

		end
	end
end
