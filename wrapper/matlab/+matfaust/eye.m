%==========================================================================================
%> @brief Identity Faust.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b eye(m,n) or eye([m,n]) forms a M-by-N Faust F = Faust(speye(M,N)).<br/>
%> &nbsp;&nbsp;&nbsp; @b eye(m) is a short for eye(m,m).<br/>
%> &nbsp;&nbsp;&nbsp; @b eye(S, 'complex') with S the size, does the same as above but returns a complex Faust.<br/>
%> &nbsp;&nbsp;&nbsp; @b eye(S, 'complex', 'dev', 'gpu') or eye(S, 'dev', 'gpu') same as above but creates the Faust on GPU.<br/>
%>
%> @param 'dev',str 'gpu or 'cpu' to create the Faust on CPU or GPU (by default on CPU).
%> @param 'dtype',str 'double' (by default) or 'float' to select the scalar type used for the Faust generated.
%>
%> @b Example
%>
%> @code
%> % in a matlab terminal
%>>> matfaust.eye(4)
%>
%>ans =
%>
%>Faust size 4x4, density 0.25, nnz_sum 4, 1 factor(s):
%>- FACTOR 0 (real) SPARSE, size 4x4, density 0.25, nnz 4
%>identity matrix flag
%>
%>>> full(matfaust.eye(4))
%>
%>ans =
%>
%>     1     0     0     0
%>     0     1     0     0
%>     0     0     1     0
%>     0     0     0     1
%>
%>>> full(matfaust.eye(4,5))
%>
%>ans =
%>
%>     1     0     0     0     0
%>     0     1     0     0     0
%>     0     0     1     0     0
%>     0     0     0     1     0
%>
%>>> full(matfaust.eye([5,4]))
%>
%>ans =
%>
%>     1     0     0     0
%>     0     1     0     0
%>     0     0     1     0
%>     0     0     0     1
%>     0     0     0     0
%>
%>>> full(matfaust.eye([5,4],'complex'))
%>
%>ans =
%>
%>   1.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
%>   0.0000 + 0.0000i   1.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
%>   0.0000 + 0.0000i   0.0000 + 0.0000i   1.0000 + 0.0000i   0.0000 + 0.0000i
%>   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   1.0000 + 0.0000i
%>   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
%>
%>>> full(matfaust.eye([4],'complex'))
%>
%>ans =
%>
%>   1.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
%>   0.0000 + 0.0000i   1.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
%>   0.0000 + 0.0000i   0.0000 + 0.0000i   1.0000 + 0.0000i   0.0000 + 0.0000i
%>   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   1.0000 + 0.0000i
%> @endcode
%>
%==========================================================================================
function F = eye(varargin)
	if(nargin < 1)
		error('First argument is mandatory')
	end
	import matfaust.Faust
	if(ismatrix(varargin{1}))
		shape = varargin{1};
		ndim = size(shape,2);
		nrows = size(shape,1);
		if(ndim > 2)
			error('N-dimensional arrays are not supported.')
		elseif(nrows > 1)
			error('Size vector should be a row vector with real elements.')
		elseif(ndim == 2)
			m = shape(1);
			n = shape(2);
		elseif(ndim == 1)
			m = varargin{1};
			if(nargin > 1 && isnumeric(varargin{2}))
				n = varargin{2};
			else
				n = m;
			end
		else
			error('Size vector should be a row vector with real elements.')
		end
	else
		error('Size inputs must be numeric.')
	end
	dev = 'cpu';
	argc = length(varargin);
	is_real = true;
	dtype = 'double';
	if(argc > 0)
		i = 2;
		while(i <= argc)
			if(argc > i)
				% next arg (value corresponding to the key varargin{i})
				tmparg = varargin{i+1};
			end
			switch(varargin{i})
				case 'complex'
					is_real = false;
				case 'real'
					is_real = true;
				case 'dev'
					if(argc == i || ~ strcmp(tmparg, 'cpu') && ~ startsWith(tmparg, 'gpu'))
						error('dev keyword argument is not followed by a valid value: cpu, gpu*.')
					else
						dev = tmparg;
					end
					i = i + 1;
				case 'dtype'
					if(argc == i || ~ strcmp(tmparg, 'float') && ~ strcmp(tmparg, 'double'))
						error('dtype keyword argument is not followed by a valid value: float, double.')
					else
						dtype = tmparg;
					end
					i = i + 1;
				otherwise
					if((isstr(varargin{i}) || ischar(varargin{i}))  && ~ strcmp(tmparg, 'cpu') && ~ startsWith(tmparg, 'gpu') && ~ strcmp(tmparg, 'complex') && ~ strcmp(tmparg, 'float') && ~ strcmp(tmparg, 'double'))
						error([ tmparg ' unrecognized argument'])
					end
			end
			i = i + 1;
		end
	end
	if(strcmp(dev, 'cpu'))
		if(is_real)
			if(strcmp(dtype, 'float'))
				core_obj = mexFaustRealFloat('eye', m, n);
			else
				core_obj = mexFaustReal('eye', m, n);
			end
		else
			core_obj = mexFaustCplx('eye', m, n);
		end
	else
		if(is_real)
			if(strcmp(dtype, 'float'))
				core_obj = mexFaustGPURealFloat('eye', m, n);
			else
				core_obj = mexFaustGPUReal('eye', m, n);
			end
		else
			core_obj = mexFaustGPUCplx('eye', m, n);
		end
	end
	e = MException('FAUST:OOM', 'Out of Memory');
	if(core_obj == 0)
		throw(e)
	end
	F = matfaust.Faust(core_obj, is_real, dev, dtype);
end
