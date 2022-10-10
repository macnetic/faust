%=============================================================
%> Returns true if obj is a LazyLinearOp, false otherwise.
%=============================================================
function B = isLazyLinearOp(obj)
	import matfaust.lazylinop.LazyLinearOp
	B = LazyLinearOp.isLazyLinearOp(obj);
end

