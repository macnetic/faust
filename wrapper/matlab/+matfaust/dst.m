%=========================================
%> @brief Returns the Direct Sine Transform (Type II) Faust of order n.
%>
%> The aialytical formula of DST II used here is:
%> \f$2 \sum_{i=0}^{n-1} x_i sin \left( {\pi (k+1) (2i + 1)} \over {2n} \right)\f$
%>
%> @param n: the order of the DST (it must be a power of two).
%> @param 'dev', str: 'gpu' or 'cpu' to create the Faust on CPU or GPU ('cpu' is the default).
%> @param 'normed',bool: true (by default) to normalize the returned Faust as if Faust.normalize() was called, false otherwise.
%>
%> @b Example
%> @code
%> % in a matlab terminal
%> >> import matfaust.*
%> >> F = dst(8);
%> >> x = 1:8;
%> >> % apply the DST to x
%> >> F*x'
%>
%> ans =
%>
%>    72.0000
%>   -25.7693
%>          0
%>    -2.6938
%>          0
%>    -0.8036
%>          0
%>    -0.2028
%>
%> >> % check the density with a larger DST Faust of size 1024
%> >> density(dst(1024))
%>
%> ans =
%>
%>     0.1581
%> >> % it is smaller than 1
%> @endcode
%>
%>@b See also matfaust.dft, matfaust.dct, matfaust.Faust.density
%=========================================
function D = dst(n, varargin)
    import matfaust.Faust
	% check n (must be integer > 0)
	if(~ isreal(n) || n < 0 || abs(n-floor(n)) > 0)
		error('n must be an integer greater than zero')
	end
	log2n = floor(log2(n));
	if(2^log2n < n)
		error('n must be a power of 2')
	end
	normed = true; % normalization by default
	dev = 'cpu';
	argc = length(varargin);
	class = 'double';
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
				case 'class'
					if(argc == i || ~ strcmp(tmparg, 'double') && ~ startsWith(tmparg, 'single'))
						error('class keyword argument is not followed by a valid value: single or double.')
					else
						class = tmparg;
					end
				otherwise
					if((isstr(varargin{i}) || ischar(varargin{i}))  && ~ strcmp(tmparg, 'cpu') && ~ startsWith(tmparg, 'gpu'))
						error([ tmparg ' unrecognized argument'])
					end
			end
		end
	end
    % bit reversal permutation
    MDFT = mod_fft(n, 'dev', dev);
    d1 = zeros(1, n);
    for k=1:n
        d1(k) = -1j * exp(-1j * pi / 2 / n * k);
    end
    % D1 = -2 * sparse(diag(d1));
    D1 = -2 * sparse(1:n, 1:n, d1, n, n, n);
    d2 = zeros(1, n);
    for k=1:n
        d2(k) = exp(-1j * pi / n * k);
    end
    % D2 = sparse(diag(d2));
    D2 = sparse(1:n, 1:n, d2, n, n, n);
    %    P_ = zeros(n*2, n);
    %    for i=1:n/2
    %        P_(i, 2*(i-1)+1) = 1;
    %        P_(i+n, 2*i) = 1;
    %    end
    % P_ as sparse matrix
    P_row_inds = [ 1:n/2, 1+n:n/2+n ];
    P_col_inds = [ 1:2:n-1, 2:2:n ];
    P_ = sparse(P_row_inds, P_col_inds, ones(n, 1), n*2, n, n);
    F_even = Faust(D1, 'dev', dev) * MDFT;
    F_odd = Faust(D1, 'dev', dev) * Faust(D2, 'dev', dev) * MDFT;
    F = [F_even, F_odd];
    F = F * Faust(P_);
    D = real(F);
	if(normed)
		D = normalize(D);
	end
	if(strcmp(dev, 'gpu'))
		D = clone(D, 'dev', 'gpu');
	end
	if(strcmp(class, 'single'))
		D = single(D);
	end
end

function O = omega(N)
    %% List of n-th root of unit raised to the power of -(k+1) (instead of k as
    %in FFT, the purpose is to write the DST).
    o = exp(pi*1j/N);
    vo = zeros(1, N);
    for k=1:N
        vo(k) = o^-k;
    end
    % O = sparse(diag(vo));
    O = sparse(1:N, 1:N, vo, N, N, N);
end

function  B = butterfly_(N)
    %% Butterfly factor of order N.
    I_N2 = speye(N/2);
    O_N2 = omega(N/2);
    B = [[ I_N2, O_N2]; [I_N2, - O_N2]];
end

function F = mod_fft(N, varargin)
	import matfaust.bitrev_perm
    %% Modified FFT for DST computation (it would have been the standard FFT
    %% if omega was raised to the power of -k instead of -(k+1)).
    N_ = N;
    Bs = cell(1, log2(N));
    i = 1;
    while N_ ~= 1
        B = butterfly_(N_);
        Bs{i} = kron(speye(N/N_), B);
        i = i + 1;
        N_ = N_ / 2;
    end
    F = matfaust.Faust([Bs, {complex(bitrev_perm(N))}], varargin{:});
end
