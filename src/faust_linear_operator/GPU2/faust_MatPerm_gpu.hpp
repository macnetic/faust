namespace Faust
{

	template<typename FPP>
		MatPerm<FPP, GPU2>::MatPerm(const MatSparse<FPP, Cpu> &factor) : MatPerm()
	{
		MatPerm<FPP, Cpu> cpu_bmat(factor);
		auto cpu_d = cpu_bmat.getD();
		d = Vect<FPP, GPU2>(cpu_d.rows(), cpu_d.diagonal().data());
		auto sd_ids_vec = cpu_bmat.get_perm_ids();
		perm_ids = new int[sd_ids_vec.size()];
		memcpy(perm_ids, sd_ids_vec.data(), sizeof(int) * sd_ids_vec.size());
		dt.resize(0);
	}

	template<typename FPP>
		void MatPerm<FPP, GPU2>::setZeros()
		{
			throw std::runtime_error("setZeros is not available on a MatPerm matrix.");
		}


	template<typename FPP>
		void MatPerm<FPP, GPU2>::multiply(MatDense<FPP, GPU2> &other, const char op_this)
		{
			if(op_this != 'N' && op_this != 'T')
				throw std::runtime_error("MatButtermfly::multiply only handle 'N' and 'T' for op_this");
			bool use_dt = is_transp ^ op_this == 'T';
			other.eltwise_mul(use_dt?dt:d, perm_ids);
		}

	template<typename FPP>
		size_t MatPerm<FPP, GPU2>::getNBytes() const
		{
			return (d.size() + (is_transp?dt.size():0)) * sizeof(FPP) + d.size() * sizeof(int);
		}


	template<typename FPP>
		void MatPerm<FPP, GPU2>::multiply(MatSparse<FPP, GPU2>& M, const char opThis)
		{
			MatDense<FPP, GPU2> Y(M);
			this->multiply(Y, opThis);
			M = Y;
		}


	template<typename FPP>
		MatSparse<FPP, GPU2> MatPerm<FPP, GPU2>::toMatSparse() const
		{
			MatSparse<FPP, GPU2> sp(this->getNbRow(), this->getNbCol());
			sp.setEyes();
			const_cast<MatPerm<FPP, GPU2>*>(this)->multiply(sp, 'N');
			return sp;
		}


	template<typename FPP>
		MatSparse<FPP,GPU2>* MatPerm<FPP, GPU2>::get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const
		{
			return toMatSparse().get_cols(col_id_start, num_cols);
		}

	template<typename FPP>
		MatSparse<FPP,GPU2>* MatPerm<FPP, GPU2>::get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const
		{
			return toMatSparse().get_rows(row_id_start, num_rows);
		}

	template<typename FPP>
		MatSparse<FPP,GPU2>* MatPerm<FPP, GPU2>::get_cols(faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const
		{
			return toMatSparse().get_cols(col_ids, num_cols);
		}

	template<typename FPP>
		MatSparse<FPP,GPU2>* MatPerm<FPP, GPU2>::get_rows(faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const
		{
			return toMatSparse().get_rows(row_ids, num_rows);
		}


	template<typename FPP>
		Real<FPP> MatPerm<FPP, GPU2>::norm() const
		{
			return toMatSparse().norm();
		}

	template<typename FPP>
		MatPerm<FPP,GPU2>* MatPerm<FPP, GPU2>::clone(const int32_t dev_id/*=-1*/, const void* stream/*=nullptr*/) const
		{
			return new MatPerm<FPP, GPU2>(*this);
		}

	template<typename FPP>
		MatPerm<FPP,GPU2>* MatPerm<FPP, GPU2>::Clone(const bool isOptimize/*=false*/) const
		{
			if (isOptimize)
			{
				throw std::runtime_error("MatPerm doesn't handle isOptimize flag");
			} else
			{
				return clone();
			}
		}


	template<typename FPP>
		faust_unsigned_int MatPerm<FPP, GPU2>::getNonZeros() const
		{
			return d.getNonZeros();
		}


	template<typename FPP>
		void MatPerm<FPP, GPU2>::transpose()
		{
			init_transpose(); // free cost if already called once
			is_transp = ! is_transp;
		}


	template<typename FPP>
		void MatPerm<FPP, GPU2>::init_transpose()
		{
			//TODO: simplify in case of symmetric matrix (it happens for the FFT)
			if(dt.size() == 0)
			{
				dt.resize(d.size());
				perm_ids_T = new int[d.size()];
				Vect<FPP, Cpu> d_cpu = d.tocpu();
				Vect<FPP, Cpu> dt_cpu(d.size());
				auto dt_ptr = dt_cpu.getData();
				auto d_ptr = d_cpu.getData();
				for(int i=0;i < d.size(); i++)
				{
					auto cid = perm_ids[i];
					perm_ids_T[cid] = i;
					dt_ptr[cid] = d_ptr[i];
				}
				dt = dt_cpu;
			}
		}

	template<typename FPP>
		void MatPerm<FPP, GPU2>::conjugate()
		{
			d.conjugate();
			if(dt.size() > 0)
				dt.conjugate();
		}


	template<typename FPP>
		void MatPerm<FPP, GPU2>::adjoint()
		{
			transpose();
			conjugate();
		}


	template<typename FPP>
		void MatPerm<FPP, GPU2>::Display() const
		{
			//TODO: adjust consistently with MatGeneric Display (using to_string)
			std::cout << "MatPerm on GPU: ";
			std::cout << "D: ";
			d.tocpu().Display();
			if(dt.size() > 0)
			{
				std::cout << "DT: ";
				dt.tocpu().Display();
			}
			std::cout << "perm_ids: ";
			for(int i=0;i < d.size();i++)
				std::cout << perm_ids[i] << " ";
			std::cout << std::endl;
		}

	template<typename FPP>
		MatDense<FPP, GPU2> MatPerm<FPP, GPU2>::multiply(const FPP* X, int X_ncols)
		{
			MatDense<FPP, GPU2> gpu_X(d.size(), X_ncols, X);
			return multiply(gpu_X);
		}

	template<typename FPP>
		MatDense<FPP, GPU2> MatPerm<FPP, GPU2>::multiply(const FPP* x)
		{
			return multiply(x, 1);
		}


	template<typename FPP>
		void MatPerm<FPP, GPU2>::multiply(MatDense<FPP, GPU2> &gpu_X, MatDense<FPP, Cpu> &cpu_out)
		{
			multiply(gpu_X);
			gpu_X.tocpu(cpu_out);
		}

	template<typename FPP>
		MatDense<FPP, GPU2> MatPerm<FPP, GPU2>::multiply(MatDense<FPP, GPU2> &gpu_X)
		{
			//TODO: do we really need this method?
			if(is_transp)
				gpu_X.eltwise_mul(dt, perm_ids);
			else
				gpu_X.eltwise_mul(d, perm_ids);
			return gpu_X;
		}

	template<typename FPP>
		MatPerm<FPP, GPU2>::~MatPerm()
		{
			if(perm_ids)
				delete perm_ids;
			if(dt.size() > 0 && perm_ids_T)
				delete perm_ids_T;
		}
}
