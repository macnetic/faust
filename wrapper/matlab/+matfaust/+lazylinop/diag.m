%=============================================================
%> @brief Constructs a diagonal LazyLinearOp.
%>
%> @param v a vector for the diagonal.
%> @param k (int, optional) diagonal to place the vector on. Default is 0 (main diagonal). Negative integer for a diagonal below the main diagonal, strictly positive integer for a diagonal above.
%>
%> @retval DL the diagonal LazyLinearOp.
%>
%> @b Example
%> @code
%> >> import matfaust.lazylinop.diag
%> >> v = rand(1, 5)
%>
%> v =
%>
%>     0.6678    0.8444    0.3445    0.7805    0.6753
%>
%> >> ld1 = diag(v)
%>
%> ld1 =
%>
%>   5x5 LazyLinearOp array with no properties.
%>
%> >> full(ld1)
%> ans =
%>
%>     0.6678         0         0         0         0
%>          0    0.8444         0         0         0
%>          0         0    0.3445         0         0
%>          0         0         0    0.7805         0
%>          0         0         0         0    0.6753
%>
%> >> ld2 = diag(v, -2)
%>
%> ld2 =
%>
%>   7x7 LazyLinearOp array with no properties.
%>
%> >> full(ld2)
%>
%> ans =
%>
%>          0         0         0         0         0         0         0
%>          0         0         0         0         0         0         0
%>     0.6678         0         0         0         0         0         0
%>          0    0.8444         0         0         0         0         0
%>          0         0    0.3445         0         0         0         0
%>          0         0         0    0.7805         0         0         0
%>          0         0         0         0    0.6753         0         0
%>
%> >> ld3 = diag(v, 2)
%>
%> ld3 =
%>
%>   7x7 LazyLinearOp array with no properties.
%>
%> >> full(ld3)
%>
%> ans =
%>
%>          0         0    0.6678         0         0         0         0
%>          0         0         0    0.8444         0         0         0
%>          0         0         0         0    0.3445         0         0
%>          0         0         0         0         0    0.7805         0
%>          0         0         0         0         0         0    0.6753
%>          0         0         0         0         0         0         0
%>          0         0         0         0         0         0         0
%>
%> @endcode
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
