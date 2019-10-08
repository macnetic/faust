%==========================================================================================
%> @brief Approximates M by A S_1 … S_n B using fact_palm4msa.
%>
%>
%> @Example
%> @code
%> import matfaust.fact.fact_palm4msa
%> import matfaust.factparams.*
%>
%> p = ParamsPalm4MSA(…
%>     ConstraintList('spcol', 2, 10, 20, 'sp', 30, 20, 20),…
%>     StoppingCriterion(50), 'is_verbose', false);
%> M = rand(10,10);
%> A = matfaust.rand(10,10);
%> B = matfaust.rand(20, 10);
%> [F, lamdba] = fact_palm4msa_constends(M, p, A, B)
%>
%> assert(norm(A - factors(F,1))/norm(A) <= eps(double(1)))
%> assert(norm(B - factors(F,4))/norm(B) <= eps(double(1)))
%>
%> @endcode
%==========================================================================================
function [F, lambda] = fact_palm4msa_constends(M, p, A, varargin)
	import matfaust.factparams.*
	import matfaust.fact.fact_palm4msa
	if(~ ismatrix(A))
		error('A must be a matrix.')
	end
	consA = ConstraintList('const', A, size(A,1), size(A,2));
	new_consts = {};
	new_consts = [ {consA.clist{:}}, {p.constraints{:}} ];
	if(length(varargin) > 0)
		B = varargin{1};
		if(~ ismatrix(B))
			error('B must be a matrix.')
		end
		consB = ConstraintList('const', B, size(B,1), size(B,2));
		new_consts = [ new_consts, {consB.clist{:}} ];
	end
	new_consts = ConstraintList(new_consts{:});
	p = ParamsPalm4MSA(new_consts, p.stop_crit, 'is_update_way_R2L', p.is_update_way_R2L, ...
		'init_lambda', p.init_lambda, 'step_size', p.step_size, 'constant_step_size', ...
		p.constant_step_size, 'is_verbose', p.is_verbose);
	[F, lambda ] = fact_palm4msa(M, p);
	f1 = factors(F, 1);
	f1 = f1 / lambda;
	nF = cell(1, numfactors(F));
	nF{1} = f1;
	for i=2:numfactors(F)
		nF{i} = factors(F, i);
	end
	nF{2} = nF{2}*lambda;
	F = matfaust.Faust(nF);
end
