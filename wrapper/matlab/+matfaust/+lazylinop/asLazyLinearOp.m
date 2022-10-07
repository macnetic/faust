function L = asLazyLinearOp(obj)
	import matfaust.lazylinop.LazyLinearOp
	L = LazyLinearOp.create(obj);
end
