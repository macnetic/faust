%==========================================================================================
%> @brief Faust identity.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b eye(m,n) or eye([m,n]) forms a M-by-N Faust F = Faust(speye(M,N)).<br/>
%> &nbsp;&nbsp;&nbsp; @b eye(m) is a short for eye(m,n).<br/>
%> &nbsp;&nbsp;&nbsp; @b eye(S, 'complex') or eye(S, 'complex') or eye(S, 'complex') with S the size, does the same as above but returns a complex Faust.</br>
%>
%> @b Example
%> @code
%> % in a matlab terminal
%>>> matfaust.eye(4)
%>
%>ans =
%>
%>Faust size 4x4, density 0.25, nnz_sum 4, 1 factor(s):
%>- FACTOR 0 (real) SPARSE, size 4x4, density 0.25, nnz 4
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
	la = varargin{nargin};
	if(nargin ~= 1 && ~ isnumeric(la) && (ischar(la) || ischar(cell2mat(la))))
		% F = Faust(sparse(1:m, 1:n, 1+eps(1)*j)); % hack to avoid passing through a full matrix
		if(strcmp(la,'complex'))
			F = Faust(eye(m,n,'like', sparse(1,1,1+i)));
		elseif(strcmp(la, 'real'))
			F = Faust(speye(m,n));
		else
			if(iscell(la))
				la = cell2mat(la)
			end
			error(['Unknown option: ' la])
		end
	else
		F = Faust(speye(m,n));
	end
end
