function L = aslazylinearoperator(obj)
	import matfaust.lazylinop.LazyLinearOp
	L = LazyLinearOp.create_from_op(obj);
end
