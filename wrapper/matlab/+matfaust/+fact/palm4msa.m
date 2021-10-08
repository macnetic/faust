%==========================================================================================
%> @brief Factorizes the matrix M with Palm4MSA algorithm using the parameters set in p.
%>
%>
%> @param M the dense matrix to factorize. The class(M) value can be double or single (only if isreal(M) is true). This class has a major impact on performance.
%> @param p the matfaust.factparams.ParamsPalm4MSA instance to define the algorithm parameters.
%> @param 'backend',int (optional) the backend to use (the C++ implementation). Must be 2016 (the default) or 2020 (which should be faster for most of the factorizations).
%> @param 'gpu', bool (optional) set to true to execute the algorithm using the GPU implementation. This option is only available when backend==2020.
%>
%> @retval F the Faust object result of the factorization.
%> @retval [F, lambda] = palm4msa(M, p) when optionally getting lambda (scale).
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
%>  F = palm4msa(M, params, 'backend', 2016)
%> @endcode
%>
%> F =
%>
%> Faust size 500x32, density 0.22025, nnz_sum 3524, 2 factor(s):
%> - FACTOR 0 (real) SPARSE, size 500x32, density 0.15625, nnz 2500
%> - FACTOR 1 (real) DENSE, size 32x32, density 1, nnz 1024
%>
%>
%==========================================================================================
function  [F,lambda] = palm4msa(M, p, varargin)
	import matfaust.Faust
	import matfaust.fact.check_fact_mat
	check_fact_mat('matfaust.fact.palm4msa', M)
	if(~ isa(p ,'matfaust.factparams.ParamsPalm4MSA'))
		error('p must be a ParamsPalm4MSA object.')
	end
	if(~ p.is_mat_consistent(M))
		error('M''s number of columns must be consistent with the last residuum constraint defined in p. Likewise its number of rows must be consistent with the first factor constraint defined in p.')
	end
	mex_params = p.to_mex_struct(M);
	backend = 2016;
	nargin = length(varargin);
	gpu = false;
	is_float = strcmp(class(M), 'single');
	if(is_float)
		dtype = 'float';
	else
		dtype = 'double'; % also for complex double
	end
	if(nargin > 0)
		for i=1:nargin
			switch(varargin{i})
				case 'backend'
					if(nargin < i+1)
						error('keyword argument ''backend'' must be followed by 2016 or 2020')
					else
						backend = varargin{i+1};
					end
				case 'gpu'
					if(nargin == i || ~ islogical(varargin{i+1}))
						error('gpu keyword argument is not followed by a logical')
					else
						gpu = varargin{i+1};
					end
			end
		end
		if(~ (isscalar(backend) && floor(backend) == backend) || backend ~= 2016 && backend ~= 2020)
			backend
			error('backend must be a int equal to 2016 or 2020')
		end
		if(backend ~= 2020 && gpu == true)
			error('GPU implementation is only available for 2020 backend.')
		end
	end
	if(backend == 2016)
		if(isreal(M))
			if(is_float)
				[lambda, core_obj] = mexPalm4MSARealFloat(mex_params);
			else
				[lambda, core_obj] = mexPalm4MSAReal(mex_params);
			end
		else
			[lambda, core_obj] = mexPalm4MSACplx(mex_params);
		end
	elseif(backend == 2020)
		% no need to keep the ParamsPalm4MSA extracted/generated cell for init_facts
		% mex_params = rmfield(mex_params, 'init_facts')
		if(isreal(M))
			init_faust = matfaust.Faust(p.init_facts, 'dtype', dtype);
			if(gpu)
				if(is_float)
					[lambda, core_obj] = mexPALM4MSA2020_gpu2RealFloat(mex_params, get_handle(init_faust));
				else
					[lambda, core_obj] = mexPALM4MSA2020_gpu2Real(mex_params, get_handle(init_faust));
				end
			else
				if(is_float)
					[lambda, core_obj] = mexPALM4MSA2020RealFloat(mex_params, get_handle(init_faust));
				else
					[lambda, core_obj] = mexPALM4MSA2020Real(mex_params, get_handle(init_faust));
				end
			end
		else
			init_faust = complex(matfaust.Faust(p.init_facts));
			if(gpu)
				[lambda, core_obj] = mexPALM4MSA2020_gpu2Cplx(mex_params, get_handle(init_faust));
			else
				[lambda, core_obj] = mexPALM4MSA2020Cplx(mex_params, get_handle(complex(init_faust)));
			end
		end
	end
	F = Faust(core_obj, isreal(M), 'cpu', dtype);
end
