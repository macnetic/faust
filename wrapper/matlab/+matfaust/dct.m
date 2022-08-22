%=========================================
%> @brief Returns the Direct Cosine Transform (Type II) Faust of order n.
%>
%> The analytical formula of DCT II used here is:
%> \f$2 \sum_{i=0}^{n-1} x_i cos \left( {\pi k (2i + 1)} \over {2n} \right)\f$
%>
%> @param n: the order of the DCT (it must be a power of two).
%> @param 'dev', str: 'gpu' or 'cpu' to create the Faust on CPU or GPU ('cpu' is the default).
%> @param 'normed',bool: true (by default) to normalize the returned Faust as if Faust.normalize() was called, false otherwise.
%> @param 'class', str: 'single' or 'double'.
%>
%>
%> @b Example
%> @code
%> % in a matlab terminal
%> >> import matfaust.*
%> >> F = dct(8);
%> >> x = ones(8, 1);
%> >> % apply the DCT to x
%> >> F*x
%>
%> ans =
%>
%>   10.2517
%>    0.0000
%>    3.5999
%>    0.0000
%>    2.4054
%>    0.0000
%>    2.0392
%>    0.0000
%>
%> >> % check the density with a larger DCT Faust of size 1024
%> >> dct(1024).density()
%> ans =
%>
%>     0.0610
%> >> % it is smaller than 1
%> @endcode
%>
%>@b See also matfaust.dft, matfaust.dst, matfaust.Faust.density
%=========================================
function D = dct(n, varargin)
	import matfaust.Faust
	DFT = matfaust.dft(n, 'normed', false);
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
					if((isstr(varargin{i}) || ischar(varargin{i}))  && ~ any(strcmp(tmparg, {'cpu'})) && ~ startsWith(tmparg, 'gpu'))
						error([ tmparg ' unrecognized argument'])
					end
			end
		end
	end
	%	P_ = zeros(n);
	%	for i=1:n/2
	%		P_(i, (i-1) * 2 + 1) = 1;
	%		P_(i + n/2, n - (i-1) * 2) = 1;
	%	end
	%	P_ as sparse matrix
	P_row_inds = 1:n;
	P_col_inds = [ 1:2:n, n:-2:1 ];
	P_ = sparse(P_row_inds, P_col_inds, ones(n, 1), n, n, n);
	d = zeros(n, 1);
	for k=1:n
		d(k) = 2*exp(-1j*pi*(k-1)/2/n);
	end
	%	E = diag(d);
	% E as sparse matrix
	E = sparse(1:n, 1:n, d, n, n, n);
	f0 = E * factors(DFT, 1);
	f_end = factors(DFT, numfactors(DFT)) * P_;
	F_mid = factors(DFT, 2:numfactors(DFT)-1);
	if ~ matfaust.isFaust(F_mid)
		F_mid = Faust(F_mid);
	end
	D = Faust(f0) * F_mid * Faust(f_end);
	D = real(D);
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
