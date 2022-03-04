% =========================================================
%> The simplified parameterization class for factorizing a rectangular matrix with the hierarchical factorization algorithm.
% =========================================================
classdef ParamsHierarchicalRectMat < matfaust.factparams.ParamsHierarchical
	properties(Constant, SetAccess = protected, Hidden)
		DEFAULT_RHO = 0.8
		DEFAULT_P_CONST_FACT = 1.4
		IDX_RHO = 1
		IDX_P = 2
		OPT_ARG_NAMES3 = { 'rho', 'P' }
	end
	methods
        %==================================
        %>  Constructor for the specialized parameterization used for example in the matfaust.demo.bsl (brain souce localization).
        %> For a better understanding you might refer to [1].
        %>
        %> [1] Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse
        %> approximations of matrices and applications", Journal of Selected
        %> Topics in Signal Processing, 2016. [https://hal.archives-ouvertes.fr/hal-01167948v1]
        %>
        %> @param m: the number of rows of the input matrix.
        %> @param n: the number of columns of the input matrix.
        %> @param j: the total number of factors.
        %> @param k: the integer sparsity per column (SPCOL, matfaust.proj.spcol) applied to the
        %> @param rightmost factor (index j-1) of shape (m, n).
        %> @param s: s*m is the integer sparsity targeted (SP, matfaust.proj.sp) for all the factors from the
        %> @param second (index 1) to index j-2. These factors are square of order n.
        %> @param rho: defines the integer sparsity (SP, matfaust.proj.sp) of the i-th residual (i=0:j-2): ceil(P*rho**i).
        %> @param P: (default value is ParamsHierarchicalRectMat.DEFAULT_P_CONST_FACT) defines the integer sparsity of the i-th residual (i=0:j-2): ceil(P*rho**i).
        %==================================
		function p = ParamsHierarchicalRectMat(m, n, j, k, s, varargin)
			import matfaust.factparams.*
			% verify arguments
			args = {m,n,j,k,s};
			anames = {'m','n','j','k','s'};
			for i=1:length(args)
				if(~ isscalar(args{i}) || ~ isreal(args{i}))
					error([ anames{i} ' must be a real scalar'])
				end
			end
			m = floor(m);
			n = floor(n);
			j = floor(j);
			k = floor(k);
			s = floor(s);
			opt_arg_map = containers.Map();
			opt_arg_names = ParamsHierarchicalRectMat.OPT_ARG_NAMES3;
			if(length(varargin) > 0)
				% retrieve all optional argument key-value pairs
				ParamsFact.parse_opt_args(varargin, opt_arg_names)
			end
			% set default values for optional args
			rho = ParamsHierarchicalRectMat.DEFAULT_RHO;
			P = ParamsHierarchicalRectMat.DEFAULT_P_CONST_FACT*m^2;
			idx_rho = ParamsHierarchicalRectMat.IDX_RHO;
			idx_P = ParamsHierarchicalRectMat.IDX_P;
			% get custom values if provided
			if(opt_arg_map.isKey(opt_arg_names{idx_rho}))
				rho = opt_arg_map(opt_arg_names{idx_rho});
			end
			if(opt_arg_map.isKey(opt_arg_names{idx_P}))
				P = opt_arg_map(opt_arg_names{idx_P});
			end
			S1_cons = ConstraintInt('spcol', m, n, k);
			S_cons = {S1_cons};
			for i=1:j-2
				S_cons = [ S_cons, {ConstraintInt('sp', m, m, s*m)} ];
			end
			R_cons = {};
			for i=1:j-1
				R_cons = [ R_cons, {ConstraintInt('sp', m, m, P*rho^(i-1))} ];
			end
			stop_crit = StoppingCriterion(30);
			p = p@matfaust.factparams.ParamsHierarchical(S_cons, R_cons, stop_crit,...
			stop_crit, 'is_update_way_R2L', true, 'is_fact_side_left', true);
		end
	end
	methods(Static)
		function sp = createParams(M, p)
			import matfaust.factparams.ParamsHierarchicalRectMat
			if(~ iscell(p))
				error('p must be a cell array')
			end
			if(~ ismatrix(M) || issparse(M))
				error('M must be a full matrix')
			end
			if(length(p) < 4)
				error('p must be of length 4')
			end
			sp = ParamsHierarchicalRectMat(size(M,1), size(M,2), p{2:end});
		end
	end
end
