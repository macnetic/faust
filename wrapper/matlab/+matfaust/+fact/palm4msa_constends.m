% experimental block start
%==========================================================================================
%> @brief Approximates M by A S_1 … S_n B using palm4msa.
%>
%> @note Notice that A (resp. B) might be multiplied lambda (scale constant of PALM4MSA).
%> It happens if it is the smallest factor in memory of the output Faust
%> (so the first (resp. last) factor might differ to A (resp. B) by this multiplicative constant).
%>
%> @Example
%> @code
%>  >> import matfaust.fact.palm4msa_constends
%>  >> import matfaust.factparams.*
%>  >> rng(42)
%>  >> p = ParamsPalm4MSA(...
%>   .. ConstraintList('spcol', 2, 10, 20, 'sp', 30, 20, 20),...
%>   .. StoppingCriterion(50), 'is_verbose', false);
%>  >> M = rand(10,10);
%>  >> A = rand(10,10);
%>  >> B = rand(20, 10);
%>  >> [F, lamdba] = palm4msa_constends(M, p, A, B);
%>  >> assert(norm(A - factors(F,1))/norm(A) <= eps(double(1)))
%>  >> assert(norm(B - factors(F,4))/norm(B) <= eps(double(1)))
%>
%> @endcode
%==========================================================================================
function [F, lambda] = palm4msa_constends(M, p, A, varargin)
	import matfaust.factparams.*
	import matfaust.fact.palm4msa
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
	[F, lambda ] = palm4msa(M, p);
end
% experimental block end
