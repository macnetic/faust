%==========================================================================================
%> @brief Factorizes the matrix M with Palm4MSA algorithm using the parameters set in p.
%>
%>
%> @param M the dense matrix to factorize.
%> @param p the ParamsPalm4MSA instance to define the algorithm parameters.
%>
%> @retval F the Faust object result of the factorization.
%> @retval [F, lambda] = palm4msa(M, p) to optionally get lambda (scale).
%>
%> @b Example
%>
%> @code
%>  import matfaust.factparams.*
%>  import matfaust.fact.palm4msa
%>  M = rand(500, 32);
%>  cons = ConstraintList('splin', 5, 500, 32, 'normcol', 1, 32, 32);
%>  % or alternatively, using the projectors
%>  % import matfaust.proj.*
%>  % cons = {splin([500,32], 5), normcol([32,32], 1)};
%>  stop_crit = StoppingCriterion(200);
%>  params = ParamsPalm4MSA(cons, stop_crit, 'is_update_way_R2L', false, 'init_lambda', 1.0);
%>  F = palm4msa(M, params)
%> @endcode
%>
%> F =
%>
%> Faust size 500x32, density 0.22025, nnz_sum 3524, 2 factor(s):
%> - FACTOR 0 (real) SPARSE, size 500x32, density 0.15625, nnz 2500
%> - FACTOR 1 (real) SPARSE, size 32x32, density 1, nnz 1024
%>
%>
%==========================================================================================
function  [F,lambda] = palm4msa(M, p)
	import matfaust.Faust
	import matfaust.fact.check_fact_mat
	mex_constraints = cell(1, length(p.constraints));
	check_fact_mat('matfaust.fact.palm4msa', M)
	if(~ isa(p ,'matfaust.factparams.ParamsPalm4MSA'))
		error('p must be a ParamsPalm4MSA object.')
	end
	for i=1:length(p.constraints)
		cur_cell = cell(1, 4);
		cur_cell{1} = p.constraints{i}.name.conv2str();
		cur_cell{2} = p.constraints{i}.param;
		cur_cell{3} = p.constraints{i}.num_rows;
		cur_cell{4} = p.constraints{i}.num_cols;
		mex_constraints{i} = cur_cell;
	end
	if(~ p.is_mat_consistent(M))
		error('M''s number of columns must be consistent with the last residuum constraint defined in p. Likewise its number of rows must be consistent with the first factor constraint defined in p.')
	end
	% put mex_constraints in a cell array again because mex eats one level of array
	mex_params = struct('data', M, 'nfacts', p.num_facts, 'cons', {mex_constraints}, 'init_facts', {p.init_facts}, 'niter', p.stop_crit.num_its, 'sc_is_criterion_error', p.stop_crit.is_criterion_error, 'sc_error_treshold', p.stop_crit.error_treshold, 'sc_max_num_its', p.stop_crit.max_num_its, 'update_way', p.is_update_way_R2L, 'grad_calc_opt_mode', p.grad_calc_opt_mode, 'constant_step_size', p.constant_step_size, 'step_size', p.step_size, 'verbose', p.is_verbose);
	if(isreal(M))
		[lambda, core_obj] = mexPalm4MSAReal(mex_params);
	else
		[lambda, core_obj] = mexPalm4MSACplx(mex_params);
	end
	F = Faust(core_obj, isreal(M));
end
