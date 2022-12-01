%=============================================================
%> @brief Constructs a diagonal LazyLinearOp.
%>
%> @param v a vector for the diagonal.
%> @param k (int, optional) diagonal to place the vector on. Default is 0 (main diagonal). Negative integer for a diagonal below the main diagonal, strictly positive integer for a diagonal above.
%>
%> @retval the diagonal LazyLinearOp.
%>
%=============================================================
function DL = diag(v, varargin)

	import matfaust.lazylinop.*
    p = inputParser;
	validK = @(k) isscalar(k) && 0 == k - floor(k);
	addOptional(p, 'k', 0, validK)

	parse(p, varargin{:})
	k = p.Results.k;

    m = numel(v) + abs(k);
    v = reshape(v, numel(v), 1);

    function y = matmat(x, v, k)
        if isLazyLinearOp(x)
            DL = diag(v, k) * x;
        else
            if k > 0
                y = v .* x(k+1:k+numel(v), :);
                y = [y ; zeros(k, size(x, 2))];
            elseif k < 0
                y = v .* x(1:numel(v), :);
                y = [zeros(abs(k), size(x, 2)) ; y];
            else % k == 0
                y = v .* x(1:numel(v), :);
            end
        end
    end
	DL = LazyLinearOperator([m, m], 'matmat', @(x) matmat(x, v, k), 'rmatmat', @(x) matmat(x, v, -k));
end
