function KL = kron(A, B)
    import matfaust.lazylinop.LazyLinearOpKron
    KL = LazyLinearOpKron(A, B); 
end
