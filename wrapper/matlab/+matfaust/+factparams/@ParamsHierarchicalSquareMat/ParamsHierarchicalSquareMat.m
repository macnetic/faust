% =========================================================
%> @brief This class is deprecated, please use ParamsHierarchicalWHT instead.
%>
%> This class will be removed in a few minor versions of matfaust.
% =========================================================
classdef ParamsHierarchicalSquareMat < matfaust.factparams.ParamsHierarchicalWHT
	methods
		function p = ParamsHierarchicalSquareMat(n)
			import matfaust.factparams.*
%			d = 2^floor(n);
%			fact_cons = {};
%			res_cons = {};
%			stop_crit = StoppingCriterion(30);
%			% ENOTE: cell concatenation
%			for i=0:(n-2)
%				fact_cons = [ fact_cons, {ConstraintInt('splincol',d,d,2)} ];
%				%fact_cons = [ fact_cons, {ConstraintInt('sp',d,d,2*d)} ];
%			end
%			for i=0:(n-2)
%				res_cons = [ res_cons, {ConstraintInt('splincol',d,d,d/2^(i+1))} ];
%				%res_cons = [ res_cons, {ConstraintInt('sp',d,d,d*d/2^(i+1))} ];
%			end
%			p = p@matfaust.factparams.ParamsHierarchical(fact_cons, res_cons, stop_crit,...
%			stop_crit, 'is_update_way_R2L', true, 'packing_RL', false);
			p = p@matfaust.factparams.ParamsHierarchicalWHT(n)
			warning('ParamsHierarchicalSquareMat class is deprecated, please use ParamsHierarchicalWHT instead. ParamsHierarchicalSquareMat class will be removed in a few minor versions of matfaust.')
		end
	end
	methods(Static)
		function sp = createParams(M, p)
			import matfaust.factparams.ParamsHierarchicalWHT
			warning('ParamsHierarchicalSquareMat class is deprecated, please use ParamsHierarchicalWHT instead. ParamsHierarchicalSquareMat class will be removed in a few minor versions of matfaust.')
			sp = ParamsHierarchicalWHT.createParams(M, p);
		end
	end
end
