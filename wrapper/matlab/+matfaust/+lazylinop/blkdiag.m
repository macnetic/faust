%=============================================================
%> @brief Returns the block diagonal LazyLinearOp formed of operators in varargin.
%>
%> @param varargin: cell array of operators defining the diagonal blocks.
%>
%> @retval LB the diagonal block LazyLinearOp.
%>
%> @b Example
%> @code
%> >> import matfaust.lazylinop.blkdiag
%> >> A = rand(10, 12);
%> >> B = rand(14, 18);
%> >> C = rand(14, 22);
%> >> LB = blkdiag(A, B, C)
%>
%> LB =
%>
%>   38x52 LazyLinearOp array with no properties.
%>
%> >> M = rand(size(LB, 2), 13);
%> >> LB * M;
%> @endcode
%>
%> <p>@b See @b also blkdiag matlab built-in.
%=============================================================
function LB = blkdiag(varargin)
	import matfaust.lazylinop.LazyLinearOperator
	% TODO: check varargin is composed of numerical arrays or compatible operators
	lAx = @(A, x) A * x;
	lAHx = @(A, x) A' * x;
	n = length(varargin);
	roffsets = zeros(n+1, 1);
	coffsets = zeros(n+1, 1);
	for i=2:n+1
		roffsets(i) = roffsets(i-1) + size(varargin{i-1}, 2);
		coffsets(i) = coffsets(i-1) + size(varargin{i-1}, 1);
	end
	function P = matmat(x, lmul, offsets)
		P = [];
		for i=1:n
			A = varargin{i};
			P = [P ; lmul(A, x(offsets(i)+1: offsets(i+1), :))];
		end
	end
	LB = LazyLinearOperator([coffsets(n+1), roffsets(n+1)], 'matmat', @(x) matmat(x, lAx, roffsets), 'rmatmat', @(x) matmat(x, lAHx, coffsets));
end
