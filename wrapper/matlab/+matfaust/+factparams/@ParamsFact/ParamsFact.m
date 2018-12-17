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
	properties (Constant, SetAccess = protected, Hidden)
		DEFAULT_STEP_SIZE = 10^-16
		DEFAULT_VERBOSITY = false
		DEFAULT_CONSTANT_STEP_SIZE = false
		DEFAULT_INIT_LAMBDA = 1.0
		DEFAULT_IS_UPDATE_WAY_R2L = false
		% constructor opt. arguments indices
		IDX_IS_UPDATE_WAY_R2L = 1
		IDX_INIT_LAMBDA = 2
		IDX_STEP_SIZE = 3
		IDX_CONSTANT_STEP_SIZE = 4
		IDX_VERBOSITY = 5
		% the order of names matters and must respect the indices above
		OPT_ARG_NAMES = {'is_update_way_R2L', 'init_lambda', 'step_size', 'constant_step_size', 'is_verbose'}
	end
	methods
		function p = ParamsFact(num_facts, constraints, varargin)
			import matfaust.factparams.ParamsFact
			% set default values
			is_update_way_R2L = ParamsFact.DEFAULT_IS_UPDATE_WAY_R2L;
			init_lambda = ParamsFact.DEFAULT_INIT_LAMBDA;
			step_size = ParamsFact.DEFAULT_STEP_SIZE;
			is_verbose = ParamsFact.DEFAULT_VERBOSITY;
			constant_step_size = ParamsFact.DEFAULT_CONSTANT_STEP_SIZE;
			% check mandatory arguments
			if(~ isscalar(num_facts) || ~ isreal(num_facts))
				error('matfaust.factparams.ParamsFact num_facts argument must be an integer.')
			else
				num_facts = floor(num_facts);
			end
			if(~ iscell(constraints))
				error(['matfaust.factparams.ParamsFact constraints argument must be a cell array.'])
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

			% then check validity of opt args (it's useless for default values but it's not too costfull)
			if(~ islogical(is_update_way_R2L))
				error('matfaust.factparams.ParamsFact ', p.OPT_ARG_NAMES{ParamsFact.IDX_IS_UPDATE_WAY_R2L} ,' argument (is_update_way_R2L) must be logical.')
			end
			if(~ isscalar(init_lambda))
				error('matfaust.factparams.ParamsFact ', p.OPT_ARG_NAMES{ParamsFact.IDX_INIT_LAMBDA},' argument (init_lambda) must be a scalar.')
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
			p.num_facts = num_facts;
			p.is_update_way_R2L = is_update_way_R2L;
			p.init_lambda = init_lambda;
			p.constraints = constraints;
			p.step_size = step_size;
			p.is_verbose = is_verbose;
			p.constant_step_size = constant_step_size;
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
						error(['Value above is an unrecognized keyword argument.'])
					end
					i = i + 2;
				end
			end

		end
	end
end
