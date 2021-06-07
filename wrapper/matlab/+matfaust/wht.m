%==========================================================================================
%> @brief Constructs a Faust implementing the Walsh-Hadamard Transform of order n.
%>
%> The resulting Faust has log2(n) sparse factors of order n, each one having 2 non-zero elements
%> per row and per column.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b H = wht(n) <br/>
%> &nbsp;&nbsp;&nbsp; @b H = wht(n, 'normed', bool)
%> &nbsp;&nbsp;&nbsp; @b H = wht(n, 'normed', bool, 'dev', str) str might be 'cpu' or 'gpu'.
%>
%> @param n the power of two exponent for a Hadamard matrix of order n and a factorization into log2(n) factors.
%> @param 'normed',bool: (optional) true (by default) to normalize the returned Faust as if Faust.normalize() was called, false otherwise.
%> @param 'dev',str: (optional) 'cpu' to create a CPU Faust (default choice) and 'gpu' for a GPU Faust.
%>
%> @retval H the Faust implementing the Hadamard transform of dimension n.
%>
%> @b Example
%> @code
%> % in a matlab terminal
%> >> import matfaust.*
%> >> H = wht(1024) % is equal to
%> >> H = normalize(wht(1024, 'normed', false, 'dev', 'cpu'))
%> @endcode
%>
%>
%>H =
%>
%>Faust size 1024x1024, density 0.0195312, nnz_sum 20480, 10 factor(s):
%>- FACTOR 0 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>- FACTOR 1 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>- FACTOR 2 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>- FACTOR 3 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>- FACTOR 4 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>- FACTOR 5 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>- FACTOR 6 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>- FACTOR 7 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>- FACTOR 8 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>- FACTOR 9 (real) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%>
%==========================================================================================
function H = wht(n, varargin)
	% check n (must be integer > 0)
	if(~ isreal(n) || n < 0 || abs(n-floor(n)) > 0)
		error('n must be an integer greater than zero')
	end
	log2n = floor(log2(n));
	if(2^log2n < n)
		error('n must be a power of 2')
	end
	if(log2n > 31)
		error('Can''t handle a Hadamard Faust of order larger than 2^31')
	end
	normed = true; % normalization by default
	dev = 'cpu';
	argc = length(varargin);
	if(argc > 0)
		for i=1:2:argc
			if(argc > i)
				% next arg (value corresponding to the key varargin{i})
				tmparg = varargin{i+1};
			end
			switch(varargin{i})
				case 'normed'
					if(argc == i || ~ islogical(tmparg))
						error('normed keyword argument is not followed by a logical')
					else
						normed = tmparg;
					end
				case 'dev'
					if(argc == i || ~ strcmp(tmparg, 'cpu') && ~ startsWith(tmparg, 'gpu'))
						error('dev keyword argument is not followed by a valid value: cpu, gpu*.')
					else
						dev = tmparg;
					end
				otherwise
					if((isstr(varargin{i}) || ischar(varargin{i}))  && ~ strcmp(tmparg, 'cpu') && ~ startsWith(tmparg, 'gpu'))
						error([ tmparg ' unrecognized argument'])
					end
			end
		end
	end
	if(strcmp(dev, 'cpu'))
		core_obj = mexFaustReal('hadamard', log2n, normed);
	else
		core_obj = mexFaustGPUReal('hadamard', log2n, normed);
	end
	is_real = true;
	e = MException('FAUST:OOM', 'Out of Memory');
	if(core_obj == 0)
		throw(e)
	end
	H = matfaust.Faust(core_obj, is_real, dev);
end
