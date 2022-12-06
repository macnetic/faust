%=============================================================
%> @brief Creates a LazyLinearOp based on the object obj which must be of a linear operator compatible type.
%>
%> @note obj must support operations and attributes defined in the LazyLinearOp class.
%> Any operation not supported would raise an exception at the evaluation time.
%> @param obj: the root object on which the LazyLinearOp is based (it could
%> be a dense matrix, a sparse matrix, a Faust object or almost any
%> object that supports the same kind of functions).
%>
%> @retval L: a LazyLinearOp instance based on obj.
%>
%> @b Example
%>
%> @code
%> >> import matfaust.lazylinop.aslazylinearoperator
%> >> M = rand(10, 12);
%> >> lM = aslazylinearoperator(M);
%> >> twolM = lM + lM
%>
%> twolM =
%>
%>   10x12 LazyLinearOp array with no properties.
%>
%> >> F = matfaust.rand(10, 12);
%> >> lF = aslazylinearoperator(F);
%> >> twolF = lF + lF
%> twolF =
%>
%>   10x12 LazyLinearOp array with no properties.
%> @endcode
%>
%> @b See @b also: matfaust.rand.
%=============================================================
function L = aslazylinearoperator(obj)
	import matfaust.lazylinop.LazyLinearOp
	L = LazyLinearOp.create_from_op(obj);
end
