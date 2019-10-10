% =========================================================
%> @brief The simplified parameterization class for factorizing a square matrix (of order a power of two) with the hierarchical factorization algorithm.
%>
%> This type for parameters is typically used for Hadamard matrix factorization.
%>
%> <p>@b See @b also matfaust.demo.hadamard, matfaust.FaustFactory.fact_hierarchical</p>
% =========================================================
classdef ParamsHierarchicalSquareMat < matfaust.factparams.ParamsHierarchical
	methods
		function p = ParamsHierarchicalSquareMat(n)
			import matfaust.factparams.*
			d = 2^floor(n);
			fact_cons = {};
			res_cons = {};
			stop_crit = StoppingCriterion(30);
			% ENOTE: cell concatenation
			for i=0:(n-2)
				fact_cons = [ fact_cons, {ConstraintInt('splincol',d,d,2)} ];
				%fact_cons = [ fact_cons, {ConstraintInt('sp',d,d,2*d)} ];
			end
			for i=0:(n-2)
				res_cons = [ res_cons, {ConstraintInt('splincol',d,d,d/2^(i+1))} ];
				%res_cons = [ res_cons, {ConstraintInt('sp',d,d,d*d/2^(i+1))} ];
			end
			p = p@matfaust.factparams.ParamsHierarchical(fact_cons, res_cons, stop_crit,...
			stop_crit, 'is_update_way_R2L', true);
		end
	end
	methods(Static)
		function sp = createParams(M, p)
			import matfaust.factparams.ParamsHierarchicalSquareMat
			pot = log2(size(M,1));
			if(size(M,1) ~= size(M,2) || pot-floor(pot) > 0)
				error('M must be a square matrix of order a power of two.')
			end
			sp = ParamsHierarchicalSquareMat(pot);
		end
	end
end
