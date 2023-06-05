% experimental block start
%===================================================================================
%> @brief Computes the FGFT for the Fourier matrix U which should be the eigenvectors of the Laplacian Lap.
%>
%> @note this algorithm is a variant of fact.hierarchical.
%>
%> @param Lap The Laplacian matrix.
%> @param U The Fourier matrix.
%> @param p The PALM hierarchical algorithm parameters.
%> @param init_D The initial diagonal vector. If none it will be the ones() vector by default.
%>
%> @retval [Uhat, Dhat, lambda, p]
%> - Uhat: the Faust factorization of U.
%> - Dhat: the diagonal matrix approximation of eigenvalues.
%> - lambda: see fact.hierarchical
%> - p: see fact.hierarchical
%>
%>
%> @b Example
%> @code
%>  >> import matfaust.*
%>  >> import matfaust.factparams.*
%>  >> import matfaust.fact.fgft_palm
%>  >> % get the Laplacian
%>  >> load('Laplacian_128_ring.mat');
%>  >> [U, D] = eig(Lap);
%>  >> [D, I] = sort(diag(D));
%>  >> D = diag(D);
%>  >> U = U(:,I);
%>  >> dim = size(Lap, 1);
%>  >> nfacts = round(log2(dim)-3);
%>  >> over_sp = 1.5; % sparsity overhead
%>  >> dec_fact = .5; % decrease of the residuum sparsity
%>  >> % define the sparsity constraints for the factors
%> >> fact_cons = {};
%> >> res_cons = {};
%> >> for i=1:nfacts
%> ..    fact_cons = [ fact_cons {ConstraintInt('sp', dim, dim, min(round(dec_fact^i*dim^2*over_sp), size(Lap,1)))} ];
%> ..    res_cons = [ res_cons {ConstraintInt('sp', dim, dim, min(round(2*dim*over_sp), size(Lap, 1)))} ];
%> .. end
%> >> % set the parameters for the PALM hierarchical algo.
%> >> params = ParamsHierarchical(fact_cons, res_cons, StoppingCriterion(50), StoppingCriterion(100), 'step_size', 1e-6, 'constant_step_size', true, 'init_lambda', 1.0, 'is_fact_side_left', false);
%> >> %% compute FGFT for Lap, U, D
%> >> init_D_diag = diag(D);
%> >> [Uhat, Dhat, lambda, ~ ] = fgft_palm(U, Lap, params, init_D_diag);
%> Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 1/4
%> Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 2/4
%> Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 3/4
%> Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 4/4
%>
%> >> err_Lap = norm(Uhat*full(Dhat)*Uhat'-Lap, 'fro') / norm(Lap, 'fro') % doctest: +ELLIPSIS
%>
%>     err_Lap =
%>
%>      	0.9...
%>
%> >>
%> @endcode
%>
%> <p> @b See @b also fact.fact_hierarchical, fact.eigtj, fact.fgft_givens
%>
%> @b References:
%> - [1]   Le Magoarou L., Gribonval R. and Tremblay N., "Approximate fast
%> graph Fourier transforms via multi-layer sparse approximations",
%> IEEE Transactions on Signal and Information Processing
%> over Networks 2018, 4(2), pp 407-420
%> <https://hal.inria.fr/hal-01416110>
%> - [2] Le Magoarou L. and Gribonval R., "Are there approximate Fast
%> Fourier Transforms on graphs ?", ICASSP, 2016.  <https://hal.inria.fr/hal-01254108>
%>
%===================================================================================
function varargout = fgft_palm(U, Lap, p, varargin)
	import matfaust.Faust
	import matfaust.factparams.*
	import matfaust.fact.check_fact_mat
	if(~ ismatrix(U) || ~ isnumeric(U) || ~ ismatrix(Lap) || ~ isnumeric(Lap))
		error('U and Lap must be real or complex matrices.')
	elseif(any(size(U) ~= size(Lap)) || any(size(Lap,1) ~= size(Lap,2)))
		error('U and Lap must be square matrices of same size.')
	end
	% TODO: refactor with fact_hierarchical
	if(length(varargin) == 1)
		init_D = varargin{1};
		if(~ ismatrix(init_D) || ~ isnumeric(init_D))
			error('fgft_palm arg. 4 (init_D) must be a matrix')
		end
		if(size(init_D,1) ~= size(U,1))
			error('fgft_palm arg. 4 (init_D) must be a diagonal vector of size == size(U,1).')
		end
	elseif(length(varargin) > 1)
		error('fgft_palm, too many arguments.')
	else % nargin == 0
		init_D = ones(size(U,1),1);
		if(~ isreal(U))
			init_D = complex(init_D);
		end
	end
	check_fact_mat('fgft_palm', U)
	if(~ isa(p, 'ParamsHierarchical') && ParamsFactFactory.is_a_valid_simplification(p))
		p = ParamsFactFactory.createParams(U, p);
	end
	mex_constraints = cell(2, p.num_facts-1);
	if(~ isa(p ,'ParamsHierarchical'))
		error('p must be a ParamsHierarchical object.')
	end
	% main factor constraints
	for i=1:p.num_facts-1
		cur_cell = cell(1, 6);
		cur_cell{1} = p.constraints{i}.name.conv2str();
		cur_cell{2} = p.constraints{i}.param;
		cur_cell{3} = p.constraints{i}.num_rows;
		cur_cell{4} = p.constraints{i}.num_cols;
		cur_cell{5} = p.constraints{i}.normalized;
		cur_cell{6} = p.constraints{i}.pos;
		mex_constraints{1,i} = cur_cell;
	end
	% residual factors constraints
	for i=1:p.num_facts-1
		cur_cell = cell(1, 6);
		cur_cell{1} = p.constraints{i+p.num_facts-1}.name.conv2str();
		cur_cell{2} = p.constraints{i+p.num_facts-1}.param;
		cur_cell{3} = p.constraints{i+p.num_facts-1}.num_rows;
		cur_cell{4} = p.constraints{i+p.num_facts-1}.num_cols;
		cur_cell{5} = p.constraints{i}.normalized;
		cur_cell{6} = p.constraints{i}.pos;
		mex_constraints{2,i} = cur_cell;
	end
	if(~ p.is_mat_consistent(U))
		error('U''s number of columns must be consistent with the last residuum constraint defined in p. Likewise its number of rows must be consistent with the first factor constraint defined in p.')
	end
	% the setters for num_rows/cols verifies consistency with constraints
	% TODO: should use ParamsHierarchical to_mex_struct
	mex_params = struct('nfacts', p.num_facts, 'cons', {mex_constraints}, 'niter1', p.stop_crits{1}.num_its, 'niter2', p.stop_crits{2}.num_its, 'sc_is_criterion_error', p.stop_crits{1}.is_criterion_error, 'sc_error_treshold', p.stop_crits{1}.tol, 'sc_max_num_its', p.stop_crits{1}.maxiter, 'sc_is_criterion_error2', p.stop_crits{2}.is_criterion_error, 'sc_error_treshold2', p.stop_crits{2}.tol, 'sc_max_num_its2', p.stop_crits{2}.maxiter, 'nrow', p.data_num_rows, 'ncol', p.data_num_cols, 'fact_side', p.is_fact_side_left, 'update_way', p.is_update_way_R2L, 'init_D', init_D, 'verbose', p.is_verbose, 'init_lambda', p.init_lambda, 'grad_calc_opt_mode', p.grad_calc_opt_mode, 'step_size', p.step_size);
	if(isreal(U))
		[lambda, core_obj, Ddiag] = mexHierarchical_factReal(U, mex_params, Lap);
	else
		[lambda, core_obj, Ddiag] = mexHierarchical_factCplx(U, mex_params, Lap);
	end
	D = sparse(diag(Ddiag));
	F = Faust(core_obj, isreal(U), 'cpu', 'double', true); % 4th arg to copy factors
	varargout = {F, D, lambda, p};
end
% experimental block end
