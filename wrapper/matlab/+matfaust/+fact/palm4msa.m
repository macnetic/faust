%==========================================================================================
%> @brief Factorizes the matrix M with Palm4MSA algorithm using the parameters set in p.
%>
%>
%> @param M the dense matrix to factorize.
%> @param p the ParamsPalm4MSA instance to define the algorithm parameters.
%> @param 'backend',int (optional) the backend (the C++ implementation) chosen. Must be 2016 (the default) or 2020 (which should be quicker for certain configurations - e.g. factorizing a Hadamard matrix).
%> @param 'gpu', bool (optional) set to true to execute the algorithm using the GPU implementation. This options is only available when backend==2020.
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
			[lambda, core_obj] = mexPalm4MSAReal(mex_params);
		else
			[lambda, core_obj] = mexPalm4MSACplx(mex_params);
		end
	elseif(backend == 2020)
		if(isreal(M))
			if(gpu)
				[lambda, core_obj] = mexPALM4MSA2020_gpu2Real(mex_params);
			else
				[lambda, core_obj] = mexPALM4MSA2020Real(mex_params);
			end
		else
			error('backend 2020 doesn''t handle yet the complex matrices')
		end
	end
	F = Faust(core_obj, isreal(M));
end
