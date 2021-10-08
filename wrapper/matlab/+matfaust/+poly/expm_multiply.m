
%==============================================
%> @brief Computes an approximate of the action of the matrix exponential of A on B using series of Chebyshev polynomials.
%>
%> @param A the operator whose exponential is of interest (must be a symmetric positive definite sparse matrix).
%> @param B the matrix  or vector to be multiplied by the matrix exponential of A.
%> @param t (real array) time points.
%> @param 'K', integer (default value is 10) the greatest polynomial degree of the Chebyshev polynomial basis. The greater it is, the better is the approximate accuracy but note that a larger K increases the computational cost.
%> @param 'tradeoff', str (optional): 'memory' or 'time' to specify what matters the most: a small memory footprint or a small time of execution. It changes the implementation of pyfaust.poly.poly used behind. It can help when the memory size is limited relatively to the value of rel_err or the size of A and B.
%> @param 'dev', str (optional): the computing device ('cpu' or 'gpu').
%> @param 'dtype', str (optional): to decide in which data type the resulting array C will be encoded ('float' or 'double' by default).
%>
%>
%> @retval C the approximate of e^{t_k A} B. C is a tridimensional array of size (sizef(A,1), size(B,2), size(t, 1)), each slice C(:,:,i) is the action of the matrix exponentatial of A on B according to the time point t(i).
%>
%> @b Example
%> @code
%> % in a matlab terminal
%> >> import matfaust.poly.expm_multiply
%> >> L = sprand(5, 5, .2);
%> >> L = L*L';
%> >> x = rand(size(L,2), 1);
%> >> t = linspace(-.5, -.1, 3);
%> >> C = expm_multiply(L, x, t);
%> @endcode
%>
%> C =
%>
%>     0.1401    0.4218    0.9157    0.3332    0.2658
%>     0.1408    0.4218    0.9157    0.4620    0.4532
%>     0.1415    0.4218    0.9157    0.6580    0.7508
%==============================================
function C = expm_multiply(A, B, t, varargin)
	dev = 'cpu';
	rel_err = 1e-6; %TODO: to use later when capable to compute K depending on rel_err
	K = 10;
	argc = length(varargin);
	dev = 'cpu';
	tradeoff = 'time';
	dtype = 'double';
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
				case 'K'
					if(argc == i || ~ isscalar(tmparg) || ~isreal(tmparg) || tmparg < 0 || tmparg-floor(tmparg) > 0)
						error('K argument must be followed by an integer')
					else
						K = tmparg;
					end
				case 'tradeoff'
					if(argc == i || ~ isstr(tmparg))
						error('tradeoff argument must be followed by ''memory'' or ''time''')
					else
						tradeoff = tmparg;
					end
				case 'group_coeffs'
					if(argc == i || ~ islogical(tmparg))
						error('group_coeffs argument must be followed by a logical (true or false)')
					else
						group_coeffs = tmparg;
					end
				case 'dtype'
					if(argc == i || ~ strcmp(tmparg, 'float') && ~ startsWith(tmparg, 'double'))
						error('dtype keyword argument is not followed by a valid value: float or double.')
					else
						dtype = tmparg;
					end
				otherwise
					if((isstr(varargin{i}) || ischar(varargin{i}))  && ~ strcmp(tmparg, 'cpu') && ~ startsWith(tmparg, 'gpu'))
						error([ tmparg ' unrecognized argument'])
					end
			end
		end
	end
	if (strcmp(tradeoff, 'memory'))
		poly_meth = 1;
	else
		poly_meth = 2; % tradeoff is time
	end
	if(~ exist('group_coeffs'))
		if(strcmp(tradeoff, 'memory'))
			group_coeffs = false;
		else
			group_coeffs = true;
			if(poly_meth == 1)
				error('can''t use poly_meth == 2 and group_coeffs at the same time')
			end
		end
	end
	if (~ issparse(A))
		error('A must be a sparse matrix.')
	end
	if (size(A,1) ~= size(A,2))
		error('A must be symmetric positive definite')
	end
	phi = eigs(A, 1) / 2;
	T = matfaust.poly.basis(A/phi-speye(size(A)), K, 'chebyshev', 'dev', dev, 'dtype', dtype);
	if (~ ismatrix(t) || ~ isreal(t) || size(t, 1) ~= 1 && size(t, 2) ~= 1)
		error('t must be a real value or a real vector')
	end
	m = size(B, 1);
	n = size(B, 2);
	npts = numel(t);
	Y = zeros(m, n, npts);
	if (poly_meth == 2)
		TB = T*B;
	end
	if(group_coeffs)
		% poly_meth == 2
		coeffs = zeros(K+1, npts);
		for i=1:npts
			tau = t(i);
			if(tau >= 0)
				error('matfaust.poly.expm_multiply handles only negative time points')
			end
			coeffs(:, i) = calc_coeffs(tau, K, phi);
		end
		Y = reshape(matfaust.poly.poly(coeffs, TB, 'dev', dev), m, n, npts);
	else
		for i=1:npts
			tau = t(i);
			if(tau >= 0)
				error('matfaust.poly.expm_multiply handles only negative time points')
			end
			coeff = calc_coeffs(tau, K, phi);
			if(poly_meth == 2)
				Y(:, :, i) = matfaust.poly.poly(coeff, TB, 'dev', dev);
			elseif(poly_meth == 1)
				Y(:, :, i) = matfaust.poly.poly(coeff, T, 'X', B, 'dev', dev);
			else
				error('poly_meth must be 1 or 2')
			end
		end
	end
	C = Y;
end

function coeff = calc_coeffs(tau, K, phi)
	% Compute the K+1 Chebychev coefficients
	coeff = zeros(K+1, 1);
	coeff(end) = 2 * besseli(K, tau * phi, 1);
	coeff(end-1) = 2 * besseli(K-1, tau * phi, 1);
	for j=K-1:-1:1
		coeff(j) = coeff(j+2) - (2 * (j-1) + 2) / (-tau * phi) * coeff(j+1);
	end
	coeff(1) = coeff(1)*.5;
end
