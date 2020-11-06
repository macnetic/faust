#ifdef USE_GPU_MOD
template<typename FPP>
  FaustCoreCppGPU<FPP>* FaustCoreCppGPU<FPP>::randFaustGPU(unsigned int t,
            unsigned int min_num_factors, unsigned int max_num_factors,
            unsigned int min_dim_size, unsigned int max_dim_size, float density, bool per_row)
{
	Faust::TransformHelper<FPP,GPU2>* th = Faust::TransformHelper<FPP,GPU2>::randFaust(Faust::RandFaustType(t), min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row);
	if(!th) return NULL;
	FaustCoreCppGPU<FPP>* core = new FaustCoreCppGPU<FPP>(th);
	return core;
}

template<typename FPP>
FaustCoreCppGPU<FPP>::FaustCoreCppGPU(Faust::TransformHelper<FPP,GPU2> *th)
{
    this->transform = th;
}


template<typename FPP>
     FaustCoreCppGPU<FPP>* FaustCoreCppGPU<FPP>::hadamardFaustGPU(unsigned int n, const bool norma)
{
	Faust::TransformHelper<FPP,GPU2>* th = Faust::TransformHelper<FPP,GPU2>::hadamardFaust(n, norma);
	if(!th) return NULL;
	FaustCoreCppGPU<FPP>* core = new FaustCoreCppGPU<FPP>(th);
	return core;
}

template<typename FPP>
   FaustCoreCppGPU<FPP>* FaustCoreCppGPU<FPP>::fourierFaustGPU(unsigned int n, const bool norma)
{
	Faust::TransformHelper<FPP,GPU2>* th = Faust::TransformHelper<FPP,GPU2>::fourierFaust(n, norma);
	if(!th) return NULL;
	FaustCoreCppGPU<FPP>* core = new FaustCoreCppGPU<FPP>(th);
	return core;
}

template<typename FPP>
    FaustCoreCppGPU<FPP>* FaustCoreCppGPU<FPP>::eyeFaustGPU(unsigned int n, unsigned int m)
{
	Faust::TransformHelper<FPP,GPU2>* th = Faust::TransformHelper<FPP,GPU2>::eyeFaust(n, m);
	if(!th) return NULL;
	FaustCoreCppGPU<FPP>* core = new FaustCoreCppGPU<FPP>(th);
	return core;
}

template<typename FPP>
void FaustCoreCppGPU<FPP>::get_product(FPP* y_data, int y_nrows, int y_ncols)
{
    Faust::MatDense<FPP, Cpu> Y = this->transform->get_product().tocpu();
    memcpy(y_data, Y.getData(), sizeof(FPP)*y_ncols*y_nrows);
}

template<typename FPP>
FaustCoreCppGPU<FPP>* FaustCoreCppGPU<FPP>::mul_faust_gpu(FaustCoreCppGPU<FPP>* right)
{
    Faust::TransformHelper<FPP,GPU2>* th = this->transform->multiply(right->transform);
    FaustCoreCppGPU<FPP>* core = new FaustCoreCppGPU<FPP>(th);
    return core;
}

template<typename FPP>
FaustCoreCppGPU<FPP>* FaustCoreCppGPU<FPP>::mul_scal_gpu(const FPP& scal)
{
	return (FaustCoreCppGPU<FPP>*) FaustCoreCpp<FPP, GPU2>::mul_scal(scal);
}

template<typename FPP>
FaustCoreCppGPU<FPP>* FaustCoreCppGPU<FPP>::normalize_gpu(int ord) const
{
	return (FaustCoreCppGPU<FPP>*) FaustCoreCpp<FPP, GPU2>::normalize(ord);
}

template<typename FPP>
FaustCoreCppGPU<FPP>* FaustCoreCppGPU<FPP>::left_gpu(const faust_unsigned_int id) const
{
	return (FaustCoreCppGPU<FPP>*) FaustCoreCpp<FPP, GPU2>::left(id);
}

template<typename FPP>
FaustCoreCppGPU<FPP>* FaustCoreCppGPU<FPP>::right_gpu(const faust_unsigned_int id) const
{
	return (FaustCoreCppGPU<FPP>*) FaustCoreCpp<FPP, GPU2>::right(id);
}

template<typename FPP>
FaustCoreCppGPU<FPP>* FaustCoreCppGPU<FPP>::transpose_gpu() const
{
	return (FaustCoreCppGPU<FPP>*) FaustCoreCpp<FPP, GPU2>::transpose();
}

template<typename FPP>
FaustCoreCppGPU<FPP>* FaustCoreCppGPU<FPP>::conjugate_gpu() const
{
	return (FaustCoreCppGPU<FPP>*) FaustCoreCpp<FPP, GPU2>::conjugate();
}

template<typename FPP>
FaustCoreCppGPU<FPP>* FaustCoreCppGPU<FPP>::adjoint_gpu() const
{
	return (FaustCoreCppGPU<FPP>*) FaustCoreCpp<FPP, GPU2>::adjoint();
}

    template<typename FPP>
FaustCoreCppGPU<FPP>* FaustCoreCppGPU<FPP>::zpruneout_gpu(const int nnz_tres, const int npasses, const bool only_forward)
{
	return (FaustCoreCppGPU<FPP>*) FaustCoreCpp<FPP, GPU2>::zpruneout(nnz_tres, npasses, only_forward);
}

    template<typename FPP>
FaustCoreCppGPU<FPP>* FaustCoreCppGPU<FPP>::clone_gpu(int dev_id/*=-1*/) const
{
    Faust::TransformHelper<FPP,GPU2>* th = this->transform->clone(dev_id);
    FaustCoreCppGPU<FPP>* core = new FaustCoreCppGPU<FPP>(th);
    return core;
}

template<typename FPP>
FaustCoreCpp<FPP,Cpu>* FaustCoreCppGPU<FPP>::clone_cpu() const
{
    Faust::TransformHelper<FPP,Cpu>* th = this->transform->tocpu();
    auto core = new FaustCoreCpp<FPP, Cpu>(th);
    return core;
}


template<typename FPP>
FaustCoreCppGPU<FPP>* FaustCoreCppGPU<FPP>::horzcat_gpu(FaustCoreCppGPU<FPP>* right) const
{
	return (FaustCoreCppGPU<FPP>*) FaustCoreCpp<FPP, GPU2>::horzcat(right);
}

template<typename FPP>
FaustCoreCppGPU<FPP>* FaustCoreCppGPU<FPP>::vertcat_gpu(FaustCoreCppGPU<FPP>* right) const
{
	return (FaustCoreCppGPU<FPP>*) FaustCoreCpp<FPP, GPU2>::vertcat(right);
}

template<typename FPP>
void  FaustCoreCppGPU<FPP>::device_gpu(char* dev) const
{
	//TODO: get the device id
	std::strncpy(dev, "gpu", 4);
}

template<typename FPP>
FaustCoreCppGPU<FPP>* FaustCoreCppGPU<FPP>::slice_gpu(unsigned int start_row_id, unsigned int end_row_id, unsigned int start_col_id, unsigned int end_col_id) const
{
	return (FaustCoreCppGPU<FPP>*) FaustCoreCpp<FPP, GPU2>::slice(start_row_id, end_row_id, start_col_id, end_col_id);
}

template<typename FPP>
FaustCoreCppGPU<FPP>* FaustCoreCppGPU<FPP>::fancy_idx_gpu(unsigned long int* row_ids, unsigned long int
        num_rows, unsigned long int* col_ids,
        unsigned long int num_cols) const
{
	return (FaustCoreCppGPU<FPP>*) FaustCoreCpp<FPP, GPU2>::fancy_idx(row_ids, num_rows, col_ids, num_cols);
}

template<typename FPP>
void FaustCoreCppGPU<FPP>::multiply_gpu(FPP* y_data, int y_nrows, int y_ncols, FPP* x_data, int* x_row_ptr, int* x_id_col, int x_nnz, int x_nrows, int x_ncols)
{
    Faust::MatSparse<FPP, GPU2> X(x_nrows, x_ncols, x_nnz, x_data, x_row_ptr, x_id_col);
    Faust::MatDense<FPP, GPU2> Y;
    Y = this->transform->multiply(X);
	Faust::MatDense<FPP, Cpu> Y_cpu;
	Y.tocpu(Y_cpu);
    memcpy(y_data, Y_cpu.getData(), sizeof(FPP)*y_ncols*y_nrows);
}
#endif
