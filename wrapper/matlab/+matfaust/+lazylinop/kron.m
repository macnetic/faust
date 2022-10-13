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
%>>> A = rand(100, 100);
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
%> @b See @b also: LazyLinearOpKron, kron matlab built-in.
%=============================================================
function KL = kron(A, B)
    import matfaust.lazylinop.LazyLinearOpKron
    KL = LazyLinearOpKron(A, B);
end
