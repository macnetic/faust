% =========================================================
%> A ParamsHierarchical specialization for which there is no residual constraints.
% =========================================================
classdef ParamsHierarchicalNoResCons < matfaust.factparams.ParamsHierarchical
	methods
		% =========================================================
		%>	@brief Constructor.
		%>
		%>	@param fact_constraints a ConstrainstList or a cell array of matfaust.proj.proj_gen or matfaust.factparams.ConstraintGeneric that define the structure of the pyfaust.fact.hierarchical resulting Faust factors in the same order if is_fact_side_left==true, in the reverse order otherwise.
		%>	@param stop_crit1  cf. ParamsHierarchical.ParamsHierarchical
		%>	@param stop_crit2 cf. ParamsHierarchical.ParamsHierarchical
		%>	@param 'is_update_way_R2L', bool cf. ParamsHierarchical.ParamsHierarchical
		%>	@param 'init_lambda', real cf. ParamsHierarchical.ParamsHierarchical
		%>	@param 'step_size', real cf. ParamsHierarchical.ParamsHierarchical
		%>	@param 'constant_step_size', bool cf. ParamsHierarchical.ParamsHierarchical
		%>	@param 'is_fact_side_left', bool cf. ParamsHierarchical.ParamsHierarchical
		%>	@param 'is_verbose', bool cf. ParamsHierarchical.ParamsHierarchical
		%>	@param 'factor_format', str cf. ParamsHierarchical.ParamsHierarchical
		%>	@param 'packing_RL', bool cf. ParamsHierarchical.ParamsHierarchical
		%>	@param 'no_normalization', bool cf. ParamsHierarchical.ParamsHierarchical
		%>  @param 'no_lambda', bool cf. ParamsHierarchical.ParamsHierarchical
		%>	@param 'norm2_max_iter', int cf. ParamsHierarchical.ParamsHierarchical
		%>	@param 'norm2_threshold', real cf. ParamsHierarchical.ParamsHierarchical
		%>	@param 'grad_calc_opt_mode', int cf. ParamsHierarchical.ParamsHierarchical
		%>
		%>
		%> @b Example
		%> This example shows two parametrizations that are equivalent. The first one, p1, is defined trough a ParamsHierarchical instance while the second one, p2, is defined using a ParamsHierarchicalNoResCons instance.
		%> @code
		%> import matfaust.proj.*
		%> import matfaust.factparams.*
		%> import matfaust.fact.hierarchical
		%> import matfaust.wht
		%> H = full(wht(32));
		%> d = size(H, 1);
		%> n = ceil(log2(d));
		%> res_projs = {};
		%> fac_projs = {};
		%> for i=1:n
		%> 	if i == n
		%> 		res_projs = { res_projs{:}, skperm([d,d], ceil(d/2^i), 'normalized', true)};
		%> 	else
		%> 		res_projs = { res_projs{:}, proj_id([d,d])};
		%> 	end
		%> 	fac_projs = {fac_projs{:}, skperm([d, d], 2, 'normalized', true)};
		%> end
		%> stop_crit = StoppingCriterion(30);
		%> p1 = ParamsHierarchical(fac_projs, res_projs, stop_crit, stop_crit, 'is_update_way_R2L', true, 'packing_RL', false);
		%> disp("factorizing with p1 (ParamsHierarchical) into Faust F1")
		%> F1 = hierarchical(H, p1, 'backend', 2020)
		%> F1_error = norm(full(F1)-H)/norm(H)
		%> simple_projs = {fac_projs{:}, res_projs{end}};
		%> p2 = ParamsHierarchicalNoResCons(simple_projs, stop_crit, stop_crit, 'is_update_way_R2L', true, 'packing_RL', false);
		%> disp("factorizing with p2 (ParamsHierarchical) into Faust F2")
		%> F2 = hierarchical(H, p2, 'backend', 2020)
		%> F2_error = norm(full(F2)-H)/norm(H)
		%> @endcode
		%>
		%> Output:
		%> @code
		%> factorizing with p1 (ParamsHierarchical) into Faust F1
		%> Faust::hierarchical: 1/5
		%> Faust::hierarchical: 2/5
		%> Faust::hierarchical: 3/5
		%> Faust::hierarchical: 4/5
		%> Faust::hierarchical: 5/5
		%>
		%> F1 = 
		%>
		%> Faust size 32x32, density 0.34375, nnz_sum 352, 6 factor(s): 
		%> - FACTOR 0 (double) SPARSE, size 32x32, density 0.0625, nnz 64
		%> - FACTOR 1 (double) SPARSE, size 32x32, density 0.0625, nnz 64
		%> - FACTOR 2 (double) SPARSE, size 32x32, density 0.0625, nnz 64
		%> - FACTOR 3 (double) SPARSE, size 32x32, density 0.0625, nnz 64
		%> - FACTOR 4 (double) SPARSE, size 32x32, density 0.0625, nnz 64
		%> - FACTOR 5 (double) SPARSE, size 32x32, density 0.03125, nnz 32
		%>
		%> F1_error =
		%>
		%>      0
		%>
		%> factorizing with p2 (ParamsHierarchical) into Faust F2
		%> Faust::hierarchical: 1/5
		%> Faust::hierarchical: 2/5
		%> Faust::hierarchical: 3/5
		%> Faust::hierarchical: 4/5
		%> Faust::hierarchical: 5/5
		%>
		%> F2 = 
		%>
		%> Faust size 32x32, density 0.34375, nnz_sum 352, 6 factor(s): 
		%> - FACTOR 0 (double) SPARSE, size 32x32, density 0.0625, nnz 64
		%> - FACTOR 1 (double) SPARSE, size 32x32, density 0.0625, nnz 64
		%> - FACTOR 2 (double) SPARSE, size 32x32, density 0.0625, nnz 64
		%> - FACTOR 3 (double) SPARSE, size 32x32, density 0.0625, nnz 64
		%> - FACTOR 4 (double) SPARSE, size 32x32, density 0.0625, nnz 64
		%> - FACTOR 5 (double) SPARSE, size 32x32, density 0.03125, nnz 32
		%>
		%> F2_error =
		%>
		%>           0
		%>
		%> @endcode
		%>
		% =========================================================
		function p = ParamsHierarchicalNoResCons(fact_constraints, stop_crit1, stop_crit2, varargin)
			import matfaust.factparams.*
			import matfaust.proj.proj_id
			% copy fact_constraints because it'll be modified
			% TODO: ParamsFact.get_constraints
			clist = ParamsFact.get_constraints(fact_constraints);
			fact_constraints = cell(1, length(clist));
			for i=1:length(clist)
				fact_constraints{i} = clist{i};
			end
			is_fact_side_left = false;
			%TODO: is_fact_side_left should be parsed by parent class ctor (refactor parsing in a parent clas function)
			for i=1:length(varargin)
				if(strcmp(varargin{i}, 'is_fact_side_left'))
					if(length(varargin) == i+1)
						error('A value argument is missing after is_fact_side_left')
					end
					if(~ islogical(varargin{i+1}))
						error('The value of is_fact_side_left must be a logical')
					end
					is_fact_side_left = varargin{i+1};
					break
				end
			end
			% there is as many residuals as factors
			res_constraints = {};
			for i=1:length(fact_constraints)-2
				if is_fact_side_left
					res_constraints = {res_constraints{:}, proj_id([fact_constraints{end}.num_rows, fact_constraints{i}]).constraint};
				else
					res_constraints = {res_constraints{:}, proj_id([fact_constraints{i}.num_cols, fact_constraints{end}.num_cols]).constraint};
				end
			end
			res_constraints = {res_constraints{:}, fact_constraints{end}};
			fact_constraints = {fact_constraints{1, 1:end-1}};
			p = p@matfaust.factparams.ParamsHierarchical(fact_constraints, res_constraints, stop_crit1, stop_crit2, varargin{:});
		end

	end
end
