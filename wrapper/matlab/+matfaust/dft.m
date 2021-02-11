%==========================================================================================
%> @brief Constructs a Faust whose the full matrix is the Discrete Fourier Transform square matrix of order n.
%>
%> The factorization algorithm used is Cooley-Tukey (FFT).
%>
%> The resulting Faust is complex and has log2(n)+1 sparse factors whose the log2(n) first
%> have 2 nonzero elements per row and per column. The last factor is a permutation matrix.
%>
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b F = dft(n) <br/>
%> &nbsp;&nbsp;&nbsp; @b F = dft(n, normed)
%>
%> @param n: the power of two for a FFT of order n and a factorization in log2(n)+1 factors.
%> @param normed: (optional) true (by default) to normalize the returned Faust as if Faust.normalize() was called, false otherwise.
%>
%>
%> @retval F the Faust implementing the FFT transform of dimension n.
%>
%> @b Example
%> @code
%> % in a matlab terminal
%> >> import matfaust.*
%> >> F = dft(1024) % is equal to
%> >> F = normalize(dft(1024))
%> @endcode
%>
%>
%> F =
%>
%> Faust size 1024x1024, density 0.0205078, nnz_sum 21504, 11 factor(s):
%> - FACTOR 0 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%> - FACTOR 1 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%> - FACTOR 2 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%> - FACTOR 3 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%> - FACTOR 4 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%> - FACTOR 5 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%> - FACTOR 6 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%> - FACTOR 7 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%> - FACTOR 8 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%> - FACTOR 9 (complex) SPARSE, size 1024x1024, density 0.00195312, nnz 2048
%> - FACTOR 10 (complex) SPARSE, size 1024x1024, density 0.000976562, nnz 1024
%==========================================================================================
function F = dft(n, varargin)
	% check n (must be integer > 0)
	if(~ isreal(n) || n < 0 || abs(n-floor(n)) > 0)
		error('n must be an integer greater than zero')
	end
	log2n = floor(log2(n));
	if(2^log2n < n)
		error('n must be a power of 2')
	end
	if(log2n>31)
		error('Can''t handle a FFT Faust of order larger than 2^31')
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
						per_row = tmparg;
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
		core_obj = mexFaustCplx('fourier', log2n, normed);
	else
		core_obj = mexFaustGPUCplx('fourier', log2n, normed);
	end
	is_real = false;
	e = MException('FAUST:OOM', 'Out of Memory');
	if(core_obj == 0)
		throw(e)
	end
	F = matfaust.Faust(core_obj, is_real, dev);
end
