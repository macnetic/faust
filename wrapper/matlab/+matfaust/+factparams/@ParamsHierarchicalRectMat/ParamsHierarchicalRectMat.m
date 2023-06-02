% =========================================================
%> The simplified parameterization class for factorizing a rectangular matrix with the hierarchical factorization algorithm.
% =========================================================
classdef ParamsHierarchicalRectMat < matfaust.factparams.ParamsHierarchical
	properties(Constant, SetAccess = protected, Hidden)
		DEFAULT_RHO = 0.8
		DEFAULT_P_CONST_FACT = 1.4
	end
	methods
        %==================================
        %>  Constructor for the specialized parameterization used for example in the matfaust.demo.bsl (brain souce localization).
        %> For a better understanding you might refer to [1].
	%>
        %> The figure below describes the sparsity of each factor of the Faust
        %> you'll obtain using pyfaust.fact.hierarchical with a
        %> ParamsHierarchicalRectMat instance.
	%>
	%> <img src="https://faust.inria.fr/files/2022/03/ParamsHierarchicalRectMat_nnz_figure.png" width="512" height="264" style="display:block;margin-left:auto;margin-right:auto"/>
	%>
	%> The resulting Faust.nnz_sum is: \f$ \lceil P m^2 \rho^{j-2} \rceil + (j-2) s m + k n \f$
        %>
        %> @param m: the number of rows of the input matrix.
        %> @param n: the number of columns of the input matrix.
        %> @param j: the total number of factors.
        %> @param k: the integer sparsity per column (SPCOL, matfaust.proj.spcol) applied to the
        %> rightmost factor (index j-1) of shape (m, n).
        %> @param s: s*m is the integer sparsity targeted (SP, matfaust.proj.sp) for all the factors from the
        %> second (index 1) to index j-2. These factors are square of order n.
        %> @param 'rho', real: defines the integer sparsity (SP, matfaust.proj.sp) of the i-th residual (i=0:j-2): ceil(P*m^2*rho^i).
        %> @param 'P', real: (default value is ParamsHierarchicalRectMat.DEFAULT_P_CONST_FACT) defines the integer sparsity of the i-th residual (i=0:j-2): ceil(P*m^2*rho^i).
	%>
	%> @b Example:
	%> @code
	%> >> import matfaust.factparams.*
	%> >> % set p1 with m, n, j, k parameters
	%> >> p1 = ParamsHierarchicalRectMat(32, 128, 8, 4, 2);
	%> >> % now with additional optional rho and P
	%> >> p2 =  ParamsHierarchicalRectMat(32, 128, 8, 4, 2, 'rho', .7, 'P', 1.5);
	%> @endcode
	%>
	%>
    %> [1] Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse
    %> approximations of matrices and applications", Journal of Selected
    %> Topics in Signal Processing, 2016. [https://hal.archives-ouvertes.fr/hal-01167948v1]
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

			% parse optional args rho and P
			parser = inputParser;
			valid_number = @(x) isnumeric(x) && isscalar(x) && isreal(x);
			addOptional(parser,'rho', ParamsHierarchicalRectMat.DEFAULT_RHO, valid_number);
			addOptional(parser,'P', ParamsHierarchicalRectMat.DEFAULT_P_CONST_FACT, valid_number);
			parse(parser, varargin{:})
			P = parser.Results.P;
			rho = parser.Results.rho;

			S1_cons = ConstraintInt('spcol', m, n, k);
			S_cons = {S1_cons};
			for i=1:j-2
				S_cons = [ S_cons, {ConstraintInt('sp', m, m, s*m)} ];
			end
			R_cons = {};
			for i=1:j-1
				R_cons = [ R_cons, {ConstraintInt('sp', m, m, ceil(P*m^2*rho^(i-1)))} ];
			end
			stop_crit = StoppingCriterion(30);
			p = p@matfaust.factparams.ParamsHierarchical(S_cons, R_cons, stop_crit,...
			stop_crit, 'is_update_way_R2L', true, 'is_fact_side_left', true);
		end
	end
	methods(Static)
		%========================
		%>
		%> Static member function to create a ParamsHierarchicalRectMat instance by a simplified parameterization expression.
		%>
		%> @param p: a list of the form {'rectmat', j, k, s} or {'rectmat', j, k, s, 'rho', rho, 'P', P} to create a parameter
		%> instance with the parameters j, k, s and optionally rho and P (see the class constructor
		%> ParamsHierarchicalRectMat.ParamsHierarchicalRectMat for their definitions).
		%>
		%> @b Example
		%> @code
		%> >> import matfaust.factparams.ParamsHierarchicalRectMat
		%> >> num_facts = 9;
		%> >> k = 10;
		%> >> s = 8;
		%> >> p = ParamsHierarchicalRectMat.createParams(rand(256, 1024), {'rectmat', num_facts, k, s});
		%> >> p2 = ParamsHierarchicalRectMat.createParams(rand(256, 1024), {'rectmat', num_facts, k, s, 'rho', 1.2, 'P', 2});
		%> @endcode
		%========================
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
