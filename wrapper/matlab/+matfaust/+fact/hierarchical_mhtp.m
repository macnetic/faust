%==========================================================================
%> @brief Runs the MHTP-PALM4MSA hierarchical factorization algorithm on the matrix M.
%>
%> This algorithm uses the MHTP-PALM4MSA (matfaust.fact.palm4msa_mhtp) instead of only PALM4MSA as matfaust.fact.hierarchical.
%>
%> @param M the dense matrix to factorize.
%> @param hierarchical_p is a set of factorization parameters. See matfaust.fact.hierarchical.
%> @param mhtp_p the matfaust.factparams.MHTPParams instance to define the MHTP algorithm parameters.
%> @param varargin: see matfaust.fact.hierarchical for the other parameters.
%>
%> @retval F The Faust object result of the factorization.
%> @retval [F, lambda] = palm4msa(M, p) when optionally getting lambda (scale).
%>
%>@b Example
%>@code
%> >> import matfaust.fact.hierarchical_mhtp
%> >> import matfaust.factparams.ParamsHierarchical
%> >> import matfaust.factparams.StoppingCriterion
%> >> import matfaust.factparams.MHTPParams
%> >> import matfaust.proj.*
%> >> M = rand(500,32);
%> >> fact_projs = { splin([500,32], 5), sp([32,32], 96), sp([32, 32], 96)};
%> >> res_projs = { normcol([32,32], 1), sp([32,32], 666), sp([32, 32], 333)};
%> >> stop_crit1 = StoppingCriterion(200);
%> >> stop_crit2 = StoppingCriterion(200);
%> >> % 50 iterations of MHTP will run every 100 iterations of PALM4MSA (each time PALM4MSA is called by the hierarchical algorithm)
%> >> mhtp_param = MHTPParams('num_its', 150, 'palm4msa_period', 100);
%> >> param = ParamsHierarchical(fact_projs, res_projs, stop_crit1, stop_crit2);
%> >> F = hierarchical_mhtp(M, param, mhtp_param)
%> Faust::hierarchical: 1/3
%> Faust::hierarchical: 2/3
%> Faust::hierarchical: 3/3
%>
%> F =
%>
%> Faust size 500x32, density 0.189063, nnz_sum 3025, 4 factor(s):
%> - FACTOR 0 (double) SPARSE, size 500x32, density 0.15625, nnz 2500
%> - FACTOR 1 (double) SPARSE, size 32x32, density 0.09375, nnz 96
%> - FACTOR 2 (double) SPARSE, size 32x32, density 0.09375, nnz 96
%> - FACTOR 3 (double) SPARSE, size 32x32, density 0.325195, nnz 333
%>
%> >>
%>
%>@endcode
%>
%==========================================================================
function  [F,lambda] = hierarchical_mhtp(M, hierarchical_p, mhtp_p, varargin)
	hierarchical_p.use_MHTP = mhtp_p;
	[F, lambda] = matfaust.fact.hierarchical(M, hierarchical_p, varargin{:}, 'backend', 2020);
end
