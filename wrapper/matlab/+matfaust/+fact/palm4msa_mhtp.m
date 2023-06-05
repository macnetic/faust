%==========================================================================
%> @brief Runs the MHTP-PALM4MSA algorithm to factorize the matrix M.
%>
%> MHTP stands for Multilinear Hard Tresholding Pursuit. This is a generalization of the Bilinear HTP algorithm describe in [1].
%>
%> [1] Quoc-Tung Le, Rémi Gribonval. Structured Support Exploration For Multilayer Sparse Matrix Fac- torization. ICASSP 2021 - IEEE International Conference on Acoustics, Speech and Signal Processing, Jun 2021, Toronto, Ontario, Canada. pp.1-5. <a href="https://hal.inria.fr/hal-03132013/document">hal-03132013</a>
%>
%> @param M the dense matrix to factorize.
%> @param palm4msa_p the matfaust.factparams.ParamsPalm4MSA instance to define the PALM4MSA algorithm parameters.
%> @param mthp_p the matfaust.factparams.MHTPParams instance to define the MHTP algorithm parameters.
%> @param varargin: see matfaust.fact.palm4msa for the other parameters.
%>
%> @retval F the Faust object result of the factorization.
%> @retval [F, lambda] = palm4msa(M, p) when optionally getting lambda (scale).
%>
%> @b Example
%> @code
%> >> % in a matlab terminal
%> >> import matfaust.fact.palm4msa_mhtp
%> >> import matfaust.factparams.*
%> >> import matfaust.proj.*
%> >> M = rand(500,32);
%> >> projs = { splin([500,32], 5), normcol([32,32], 1.0)};
%> >> stop_crit = StoppingCriterion(200);
%> >> param = ParamsPalm4MSA(projs, stop_crit);
%> >> mhtp_param = MHTPParams('num_its', 60, 'palm4msa_period', 10);
%> >> G = palm4msa_mhtp(M, param, mhtp_param)
%>
%> G =
%>
%> Faust size 500x32, density 0.18225, nnz_sum 2916, 2 factor(s):
%> - FACTOR 0 (double) SPARSE, size 500x32, density 0.15625, nnz 2500
%> - FACTOR 1 (double) SPARSE, size 32x32, density 0.40625, nnz 416
%>
%> >>
%> @endcode
%>
%==========================================================================
function  [F,lambda] = palm4msa_mhtp(M, palm4msa_p, mhtp_p, varargin)
	palm4msa_p.use_MHTP = mhtp_p;
	[F, lambda] = matfaust.fact.palm4msa(M, palm4msa_p, 'backend', 2020, varargin{:});
end
