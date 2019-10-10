%==========================================================================================
%> @brief Approximates M by A S_1 ... S_n B using hierarchical.
%>
%> @Example
%> @code
%> import matfaust.fact.hierarchical
%> import matfaust.factparams.*
%>
%> p = ParamsHierarchical(…
%>     ConstraintList('spcol', 2, 10, 20, 'sp', 30, 10, 10), ConstraintList('sp', 4, 10, 20, 'splin', 5, 10, 10),…
%>     StoppingCriterion(50), StoppingCriterion(50),…
%>     'is_fact_side_left', true, 'is_verbose', false…
%>     );
%> M = rand(10,10);
%> A = matfaust.rand(10,10);
%> B = matfaust.rand(20, 10);
%> [F, lamdba, ~] = hierarchical_constends(M, p, A, B)
%>
%> assert(norm(A - factors(F,1))/norm(A) <= eps(double(1)))
%> assert(norm(B - factors(F,4))/norm(B) <= eps(double(1)))
%> @endcode
%>
%==========================================================================================
function varargout = hierarchical_constends(M, p, A, B)
	import matfaust.factparams.*
	import matfaust.fact.hierarchical
	if(~ ismatrix(A) || ~ ismatrix(B))
		error('A and B must be matrices.')
	end
	consA = ConstraintList('const', A, size(A, 1), size(A, 2));
	consB = ConstraintList('const', B, size(B, 1), size(B, 2));
	consts = p.constraints;
	nconsts = length(p.constraints);
	% consts: factor constraints + residuum constraints
	fac_cons = {};
	res_cons = {};
	for i=1:p.num_facts-1
		fac_cons = { fac_cons{:}, consts{i} };
	end
	for i=p.num_facts:length(consts)
		res_cons = { res_cons{:}, consts{i} };
	end
	assert(length(fac_cons) == length(res_cons))
	% add two constants factor constraints for A and B to the old constraints
	% According to the factorization direction, switch A and B positions
	if(p.is_fact_side_left)
		new_fact_cons = { consB.clist{:}, fac_cons{:} };
		new_res_cons = { res_cons{:}, consA.clist{:} };
	else
		new_fact_cons = { consA.clist{:}, fac_cons{:} };
		new_res_cons = { res_cons{:}, consB.clist{:} };
	end

	p = ParamsHierarchical(new_fact_cons, new_res_cons,...
		p.stop_crits{1}, p.stop_crits{2}, 'is_update_way_R2L', p.is_update_way_R2L, ...
		'init_lambda', p.init_lambda, 'step_size', p.step_size, 'constant_step_size', ...
		p.constant_step_size, 'is_verbose', p.is_verbose, 'is_fact_side_left', p.is_fact_side_left);
	[F, lambda, p] = hierarchical(M, p);
	f1 = factors(F, 1);
	f1 = f1 / lambda;
	nF = cell(1, numfactors(F));
	nF{1} = f1;
	for i=2:numfactors(F)
		nF{i} = factors(F, i);
	end
	nF{2} = nF{2}*lambda;
	F = matfaust.Faust(nF);
	varargout = {F, lambda, p};
end
