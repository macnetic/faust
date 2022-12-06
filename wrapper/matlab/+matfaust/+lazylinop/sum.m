%=============================================================
%> @brief Sums (lazily) all linear operators in varargin.
%>
%> @param varargin: the operators to add up.
%>
%> @retval LS the LazyLinearOp for the sum.
%>
%> @b Example:
%> @code
%> >> import matfaust.dft
%> >> d = 32;
%> >> H = dft(32);
%> >> v = rand(d, 1);
%> >> D = spdiags(v, 0, d, d);
%> >> HDH = H * matfaust.Faust(D) * H;
%> >> terms = {}; for i=1:10; terms = [ terms {HDH} ];end
%> >> LS = matfaust.lazylinop.sum(terms{:});
%>
%> ans =
%>
%>     32    32
%>
%> >> M = rand(size(LS, 2), 13);
%> >> LS * M
%> @endcode
%=============================================================
function LS = sum(varargin)
	import matfaust.lazylinop.LazyLinearOperator
	% TODO: check varargin is composed of numerical arrays or compatible operators
	lAx = @(A, x) A * x;
	lAHx = @(A, x) A' * x;
	n = length(varargin);
	function P = matmat(x, lmul)
		P = lmul(varargin{1}, x);
		for i=2:n
			A = varargin{i};
			P = P + lmul(A, x);
		end
	end
	LS = LazyLinearOperator(size(varargin{1}), 'matmat', @(x) matmat(x, lAx), 'rmatmat', @(x) matmat(x, lAHx));
end
