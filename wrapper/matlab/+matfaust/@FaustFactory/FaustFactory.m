classdef FaustFactory
	properties (SetAccess = public)

	end
	methods(Static)

		%==========================================================================================
		%> @brief Factorizes the matrix M using the parameters set in p.
		%>
		%>
		%> @param m the dense matrix to factorize.
		%> @param p the ParamsPalm4MSA instance to define the algorithm parameters.
		%>
		%> @retval F the Faust object result of the factorization.
		%>
		%> @b Example
		%>
		%> @code
		%>  import matfaust.*
		%>  num_facts = 2
		%>  is_update_way_R2L = false
		%>  init_lambda = 1.0
		%>  init_facts = cell(2,1)
		%>  init_facts{1} = zeros(500,32)
		%>  init_facts{2} = eye(32)
		%>  M = rand(500, 32)
		%>  M = data
		%>  cons = cell(2,1)
		%>  cons{1} = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5);
		%>  cons{2} = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32, 32, 1.0);
		%>  stop_crit = StoppingCriterion(200);
		%>  params = ParamsPalm4MSA(num_facts, is_update_way_R2L, init_lambda, init_facts, cons, stop_crit);
		%>  F = FaustFactory.fact_palm4msa(M, params)
		%> @endcode
		%>
		%==========================================================================================
		function  F = fact_palm4msa(m, p)
			import matfaust.Faust
			mex_constraints = cell(1, length(p.constraints));
			for i=1:length(p.constraints)
				cur_cell = cell(1, 4);
				cur_cell{1} = p.constraints{i}.name.conv2str();
				cur_cell{2} = p.constraints{i}.param;
				cur_cell{3} = p.constraints{i}.num_rows;
				cur_cell{4} = p.constraints{i}.num_cols;
				mex_constraints{i} = cur_cell;
			end
			% put mex_constraints in a cell array again because mex eats one level of array
			mex_params = struct('data', m, 'nfacts', p.num_facts, 'cons', {mex_constraints}, 'init_facts', {p.init_facts}, 'niter', p.stop_crit.num_its, 'sc_is_criterion_error', p.stop_crit.is_criterion_error, 'sc_error_treshold', p.stop_crit.error_treshold, 'sc_max_num_its', p.stop_crit.max_num_its);
			[lambda, cell_facts] = mexPalm4MSA(mex_params);
			cell_facts{1} = cell_facts{1}*lambda;
			F = Faust(cell_facts);
		end

		%==========================================================================================
		%> @brief Factorizes the matrix M using the parameters set in p.
		%>
		%>
		%> @param m the dense matrix to factorize.
		%> @param p the ParamsHierarchicalFact instance to define the algorithm parameters.
		%>
		%> @retval F The Faust object result of the factorization.
		%>
		%> @Example
		%> @code
		%>  import matfaust.*;
		%>  num_facts = 4;
		%>  is_update_way_R2L = false;
		%>  init_lambda = 1.0;
		%>  init_facts = cell(num_facts,1);
		%>  init_facts{1} = zeros(500,32);
		%>  for i=2:num_facts
		%>  	init_facts{i} = zeros(32);
		%>  end
		%>  M = rand(500, 32);
		%>  cons = cell(4,1);
		%>  cons{1} = ConstraintInt(ConstraintName(ConstraintName.SPLIN), 500, 32, 5);
		%>  cons{2} = ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 96);
		%>  cons{3} = ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 96);
		%>  cons{4} = ConstraintReal(ConstraintName(ConstraintName.NORMCOL), 32, 32, 1);
		%>  cons{5} =  ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 666);
		%>  cons{6} =  ConstraintInt(ConstraintName(ConstraintName.SP), 32, 32, 333);
		%>  stop_crit = StoppingCriterion(200);
		%>  stop_crit2 = StoppingCriterion(200);
		%>  params = ParamsHierarchicalFact(num_facts, is_update_way_R2L, init_lambda, init_facts, cons, size(M,1), size(M,2), {stop_crit, stop_crit2});
		%>  F = FaustFactory.fact_hierarchical(M, params)
		%>  @endcode
		%>  Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorisation 1/3<br/>
		%>  Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorisation 2/3<br/>
		%>  Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorisation 3/3<br/>
		%>
		%>  F = 
		%>
		%>  Faust size 500x32, density 0.189063, nnz_sum 3025, 4 factor(s): 
		%>  - FACTOR 0 (real) DENSE, size 500x32, density 0.15625, nnz 2500
		%>  - FACTOR 1 (real) DENSE, size 32x32, density 0.09375, nnz 96
		%>  - FACTOR 2 (real) DENSE, size 32x32, density 0.09375, nnz 96
		%>  - FACTOR 3 (real) DENSE, size 32x32, density 0.325195, nnz 333
		%>
		%==========================================================================================
		function F = fact_hierarchical(m, p)
			import matfaust.Faust
			mex_constraints = cell(2, p.num_facts-1);
			%mex_fact_constraints = cell(1, p.num_facts-1)
			for i=1:p.num_facts-1
				cur_cell = cell(1, 4);
				cur_cell{1} = p.constraints{i}.name.conv2str();
				cur_cell{2} = p.constraints{i}.param;
				cur_cell{3} = p.constraints{i}.num_rows;
				cur_cell{4} = p.constraints{i}.num_cols;
				%mex_fact_constraints{i} = cur_cell;
				mex_constraints{1,i} = cur_cell;
			end
			%mex_residuum_constraints = cell(1, p.num_facts-1)
			for i=1:p.num_facts-1
				cur_cell = cell(1, 4);
				cur_cell{1} = p.constraints{i+p.num_facts-1}.name.conv2str();
				cur_cell{2} = p.constraints{i+p.num_facts-1}.param;
				cur_cell{3} = p.constraints{i+p.num_facts-1}.num_rows;
				cur_cell{4} = p.constraints{i+p.num_facts-1}.num_cols;
				%mex_residuum_constraints{i} = cur_cell;
				mex_constraints{2,i} = cur_cell;
			end
			mex_params = struct('data', m, 'nfacts', p.num_facts, 'cons', {mex_constraints}, 'init_facts', {p.init_facts}, 'niter1', p.stop_crits{1}.num_its,'niter2', p.stop_crits{2}.num_its, 'sc_is_criterion_error', p.stop_crits{1}.is_criterion_error, 'sc_error_treshold', p.stop_crits{1}.error_treshold, 'sc_max_num_its', p.stop_crits{1}.max_num_its, 'sc_is_criterion_error2', p.stop_crits{2}.is_criterion_error, 'sc_error_treshold2', p.stop_crits{2}.error_treshold, 'sc_max_num_its2', p.stop_crits{2}.max_num_its, 'nrow', p.data_num_rows, 'ncol', p.data_num_cols, 'fact_side', p.is_fact_side_left);
			[lambda, cell_facts] = mexHierarchical_fact(m, mex_params);
			cell_facts{1} = cell_facts{1}*lambda;
			F = Faust(cell_facts);
		end
	end
end
