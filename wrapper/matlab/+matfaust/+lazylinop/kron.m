%=============================================================
%> @brief Returns the LazyLinearOp for the Kronecker product A x B.
%>
%>
%> @note this specialization is particularly optimized for applying the operator self to a vector.
%>
%> @param A: LinearOperator (scaling factor),
%> @param B: LinearOperator (block factor).
%>
%> @b Example:
%> @code
%> >> A = rand(100, 100);
%> >> B = rand(100, 100);
%> >> AxB = kron(A, B);
%> >> lAxB = matfaust.lazylinop.kron(A, B)
%>
%> lAxB =
%>
%>   10000x10000 LazyLinearOpKron array with no properties.
%>
%> >> x = rand(size(AxB, 2), 1);
%> >> timeit(@() AxB * x)
%>
%> ans =
%>
%>     0.0400
%>
%> >> timeit(@() lAxB * x)
%>
%> ans =
%>
%>    7.8475e-04
%>
%> >> norm(AxB * x - lAxB * x) < 1e-9
%>
%> ans =
%>
%>   logical
%>
%>    1
%> @endcode
%>
%> @b See @b also:  kron matlab built-in.
%=============================================================
function LK = kron(A, B)
    import matfaust.lazylinop.*
    hasmeth = @(x, meth) any(ismember(methods(x), meth));
    function LmK = kron_(A, B, shape, op)
        op_is_scalar = all(size(op) == [1, 1]);
        if ~ op_is_scalar && ~ shape(2) == size(op, 1)
            error('Dimensions must agree')
        end
        if op_is_scalar
            new_size = shape;
        else
            new_size = [shape(1), size(op, 2)];
        end
        if ismatrix(op) && isnumeric(op) && ~ op_is_scalar && hasmeth(op, 'reshape') && hasmeth(op, 'mtimes') %&& hasmeth(op, 'subsref')
            % op is not a LazyLinearOp because isnumeric(LazyLinearOp) is always false
            % op is a dense matrix that is not limited to one element
            LmK = zeros(new_size);
            for j=1:new_size(2)
                op_mat = reshape(op(:, j), [size(B, 2), size(A, 2)]);
                LmK(:, j) = reshape(B * op_mat * transpose(A), [new_size(1), 1]);
            end
        else
            error(['op must possess reshape and mtimes methods to'
            ' be multiplied by a Kronecker LazyLinearOp (use full() on the'
            ' latter to multily by the former)'])
        end
    end
    shape = [size(A, 1) * size(B, 1), size(A, 2) * size(B, 2)];
    LK = LazyLinearOperator(shape, 'matmat', @(x) kron_(A, B, shape, x), ...
        'rmatmat', @(x) kron_(A', B', [shape(2), shape(1)], x));
end
