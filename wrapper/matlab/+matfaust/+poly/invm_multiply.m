% experimental block start

%==============================================
%> @brief Computes an approximate of the action of the matrix inverse of A on B using Chebyshev polynomials.
%>
%> @param A the operator whose inverse is of interest (must be a symmetric positive definite sparse matrix).
%> @param B the dense matrix or vector to be multiplied by the matrix inverse of A.
%> @param 'rel_err', rel_err (optional): the targeted relative error between the approximate of the action and the action itself (if you were to compute it with inv(A)*x).
%> @param 'max_K', int (optional): the maximum degree of Chebyshev polynomial to use (useful to limit memory consumption).
%> @param 'dev', str (optional) 'cpu' or 'gpu', selects the device to use (currently only 'cpu' is supported).
%>
%>
%> @b Example
%> @code
%> >> import matfaust.poly.*
%> >> A = sprand(64, 64, .1);
%> >> A = A*A';
%> >> B = rand(size(A, 2), 2);
%> >> AinvB = invm_multiply(A, B, 'rel_err', 1e-6);
%> >> AinvB_ref = inv(full(A))*B;
%> >> err = norm(AinvB-AinvB_ref)/norm(AinvB)
%> @endcode
%>
%>err =
%>
%>   8.5162e-11
%>
%==============================================
function AinvB = invm_multiply(A, B, varargin)
	dev = 'cpu';
	rel_err = 1e-6;
	max_K = inf;
	poly_meth = 1;
	argc = length(varargin);
	if(argc > 0)
		for i=1:2:argc
			if(argc > i)
				% next arg (value corresponding to the key varargin{i})
				tmparg = varargin{i+1};
			end
			switch(varargin{i})
				case 'rel_err'
					if(argc == i || ~ isscalar(tmparg) || ~ isreal(tmparg) || tmparg < 0)
						error('rel_err argument must be followed by a positive real value')
					else
						rel_err = tmparg;
					end
				case 'dev' % not used yet
					if(argc == i || ~ strcmp(tmparg, 'cpu') && ~ startsWith(tmparg, 'gpu'))
						error('dev keyword argument is not followed by a valid value: cpu, gpu*.')
					else
						dev = tmparg;
					end
				case 'max_K'
					if(argc == i || ~ isscalar(tmparg) || ~isreal(tmparg) || tmparg < 0 || tmparg-floor(tmparg) > 0)
						error('max_K argument must be followed by an integer')
					else
						max_K = tmparg;
					end
				case 'poly_meth'
					if(argc == i || ~ isscalar(tmparg) || ~isreal(tmparg) || tmparg < 0 || tmparg-floor(tmparg) > 0)
						error('poly_meth argument must be followed by an integer')
					else
						poly_meth = tmparg;
					end
				otherwise
					if((isstr(varargin{i}) || ischar(varargin{i}))  && ~ strcmp(tmparg, 'cpu') && ~ startsWith(tmparg, 'gpu'))
						error([ tmparg ' unrecognized argument'])
					end
			end
		end
	end
	if (~ issparse(A))
		error('A must be a sparse matrix.')
	end
	if (size(A,1) ~= size(A,2))
		error('A must be symmetric positive definite')
	end
	Id = speye(size(A));
	b = eigs(A, 1);
	B_ = b*Id-A;
	b_ = eigs(B_, 1);
	a = b-b_;
	if(a <= 0)
		error('A is a singular matrix or its spectrum contains negative values')
	end
	m = (a + b) / 2;
	c = (b - a) / (b + a);
	g = abs(1 / c + sqrt(1/c^2 - 1));
	K = min(max_K, floor(((log(1/rel_err) + log(2/(m*sqrt(1-c^2))) - log(g-1)) / log(g))));
	Abar = 2*A/(b-a) - (b+a)*Id/(b-a);
	T = matfaust.poly.basis(Abar, K, 'chebyshev');
	coeffs = zeros(K+1, 1);
	for k=0:K
		coeffs(k+1) = 2 / (m*sqrt(1-c^2)) * (-1)^k * g^(-k);
	end
	coeffs(1) = coeffs(1)*.5;
	if(poly_meth == 2)
		TB = T*B
		AinvB = matfaust.poly.poly(coeffs, TB);
	elseif(poly_meth == 1)
		AinvB = matfaust.poly.poly(coeffs, T, 'X', B);
	else
		error('poly_meth must be 1 or 2')
	end
end
% experimental block end
