%=========================================
%> @brief Returns the Direct Cosine Transform (Type II) Faust of order n.
%>
%> @param n: the order of the DCT (it must be a power of two).
%> @param 'dev', str: 'gpu' or 'cpu' to create the Faust on CPU or GPU ('cpu' is the default).
%>
%>
%> @b Example
%> @code
%> % in a matlab terminal
%> >> import matfaust.*
%> >> F = dst(8);
%> >> x = 1:8;
%> >> % apply the DST to x
%> >> real(F*x')
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
%> @endcode
%>
%>@b See also matfaust.dft, matfaust.dst
%=========================================
function D = dct(n, varargin)
	import matfaust.Faust
	DFT = matfaust.dft(n, varargin{:}, 'normed', false);
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
end
