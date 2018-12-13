classdef ParamsHierarchicalFactSquareMat < matfaust.factparams.ParamsHierarchicalFact
	methods
		function p = ParamsHierarchicalFactSquareMat(n)
			import matfaust.factparams.*
			d = 2^floor(n);
			fact_cons = {};
			res_cons = {};
			stop_crit = StoppingCriterion(30);
			% ENOTE: cell concatenation
			for i=0:(n-2)
				fact_cons = [ fact_cons, {ConstraintInt('splincol',d,d,2)} ];
			end
			for i=0:(n-2)
				res_cons = [ res_cons, {ConstraintInt('splincol',d,d,d/2^(i+1))} ];
			end
			p = p@matfaust.factparams.ParamsHierarchicalFact(fact_cons, res_cons, stop_crit,...
			stop_crit, 'is_update_way_R2L', true);
		end
	end
end
