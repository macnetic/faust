% ==========================================================
%> @brief This class is a simple parameterization of PALM4MSA to factorize a Hadamard
%> matrix using the matfaust.proj.skperm proximity operator.
%>
%> @b Examples
%>
%> @code
%> import matfaust.factparams.ParamsPalm4msaWHT
%> import matfaust.fact.palm4msa
%> import matfaust.wht
%> d = 128;
%> H = full(wht(d));
%> p = ParamsPalm4msaWHT(d);
%> F = palm4msa(H, p);
%> err = norm(full(F)-H)/norm(H) % should be about 1e-16
%> @endcode
%>
%> Reference:
%>    [1] Quoc-Tung Le, Rémi Gribonval. Structured Support Exploration For
%>    Multilayer Sparse Matrix Fac- torization. ICASSP 2021 - IEEE International
%>    Conference on Acoustics, Speech and Signal Processing, Jun 2021, Toronto,
%>    Ontario, Canada. pp.1-5 <a href="https://hal.inria.fr/hal-03132013/document">hal-03132013</a>.
%>
%> See also matfaust.fact.palm4msa
%>
%>
% ==========================================================
classdef ParamsPalm4msaWHT < matfaust.factparams.ParamsPalm4MSA
	properties
	end
	methods
		% =========================================================
		%>
		%> @brief Constructor.
		%>
		% ==========================================================
		function p = ParamsPalm4msaWHT(matrix_size, varargin)
			import matfaust.factparams.*
			import matfaust.proj.skperm
			fac_projs = {}
			n = log2(matrix_size)
			for i=1:n
				fac_projs = {fac_projs{:}, skperm([matrix_size, matrix_size], 2, 'normalized', false, 'pos', false)}
			end

			stop_crit = StoppingCriterion(30)
			p = p@matfaust.factparams.ParamsPalm4MSA(fac_projs, stop_crit, 'is_update_way_R2L', false, 'packing_RL', false, varargin{:});
		end
	end
end
