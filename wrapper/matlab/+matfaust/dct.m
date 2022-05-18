%=========================================
%> Returns the Direct Cosine Transform (Type II) Faust of order n.
%=========================================
function D = dct(n, varargin)
	import matfaust.Faust
	DFT = matfaust.dft(n, varargin{:}, 'normed', false);
	P_ = zeros(n);
	for i=1:n/2
		P_(i, (i-1) * 2 + 1) = 1;
		P_(i + n/2, n - (i-1) * 2) = 1;
	end
	d = zeros(n, 1);
	for k=1:n
		d(k) = 2*exp(-1j*pi*(k-1)/2/n);
	end
	E = diag(d);
	f0 = sparse(E * factors(DFT, 1));
	f_end = sparse(factors(DFT, numfactors(DFT)) * P_);
	F_mid = factors(DFT, 2:numfactors(DFT)-1);
	if ~ matfaust.isFaust(F_mid)
		F_mid = Faust(F_mid);
	end
	D = Faust(f0) * F_mid * Faust(f_end);
end
