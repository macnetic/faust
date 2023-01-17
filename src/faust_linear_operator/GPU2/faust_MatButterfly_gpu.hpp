namespace Faust
{


	template<typename FPP>
		MatButterfly<FPP, GPU2>::MatButterfly(const MatSparse<FPP, Cpu> &factor, int level) : level(level), is_transp(false)
	{
		MatButterfly<FPP, Cpu> cpu_bmat(factor, level);
		auto cpu_d1 = cpu_bmat.getD1();
		auto cpu_d2 = cpu_bmat.getD2();
		d1 = Vect<FPP, GPU2>(cpu_d1.rows(), cpu_d1.diagonal().data());
		d2 = Vect<FPP, GPU2>(cpu_d2.rows(), cpu_d2.diagonal().data());
		auto sd_ids_vec = cpu_bmat.get_subdiag_ids();
		subdiag_ids = new int[sd_ids_vec.size()];
		memcpy(subdiag_ids, sd_ids_vec.data(), sizeof(int) * sd_ids_vec.size());
		d2t.resize(0);
	}

	template<typename FPP>
		void MatButterfly<FPP, GPU2>::setZeros()
		{
			throw std::runtime_error("setZeros is not available on a MatButterfly matrix.");
		}


	template<typename FPP>
		void MatButterfly<FPP, GPU2>::multiply(MatDense<FPP, GPU2> &other, const char op_this) const
		{
			if(op_this != 'N' && op_this != 'T')
				throw std::runtime_error("MatButtermfly::multiply only handle 'N' and 'T' for op_this");
			bool use_d2t = is_transp ^ op_this == 'T';
			butterfly_diag_prod(other, d1, use_d2t?d2t:d2, subdiag_ids);
		}

	template<typename FPP>
		size_t MatButterfly<FPP, GPU2>::getNBytes() const
		{
			return (d1.size() + d2.size() + (is_transp?d2.size():0)) * sizeof(FPP) + d2.size() * sizeof(int);
		}


	template<typename FPP>
		void MatButterfly<FPP, GPU2>::multiply(MatSparse<FPP, GPU2>& M, const char opThis)
		{
			MatDense<FPP, GPU2> Y(M);
			this->multiply(Y, opThis);
			M = Y;
		}


	template<typename FPP>
		MatSparse<FPP, GPU2> MatButterfly<FPP, GPU2>::toMatSparse() const
		{
			MatSparse<FPP, GPU2> sp(this->getNbRow(), this->getNbCol());
			sp.setEyes();
			const_cast<MatButterfly<FPP, GPU2>*>(this)->multiply(sp, 'N');
			return sp;
		}


	template<typename FPP>
		MatSparse<FPP,GPU2>* MatButterfly<FPP, GPU2>::get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const
		{
			return toMatSparse().get_cols(col_id_start, num_cols);
		}

	template<typename FPP>
		MatSparse<FPP,GPU2>* MatButterfly<FPP, GPU2>::get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const
		{
			return toMatSparse().get_rows(row_id_start, num_rows);
		}

	template<typename FPP>
		MatSparse<FPP,GPU2>* MatButterfly<FPP, GPU2>::get_cols(faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const
		{
			return toMatSparse().get_cols(col_ids, num_cols);
		}

	template<typename FPP>
		MatSparse<FPP,GPU2>* MatButterfly<FPP, GPU2>::get_rows(faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const
		{
			return toMatSparse().get_rows(row_ids, num_rows);
		}


	template<typename FPP>
		Real<FPP> MatButterfly<FPP, GPU2>::norm() const
		{
			return toMatSparse().norm();
		}

	template<typename FPP>
		MatButterfly<FPP,GPU2>* MatButterfly<FPP, GPU2>::clone(const int32_t dev_id/*=-1*/, const void* stream/*=nullptr*/) const
		{
			return new MatButterfly<FPP, GPU2>(*this);
		}

	template<typename FPP>
		MatButterfly<FPP,GPU2>* MatButterfly<FPP, GPU2>::Clone(const bool isOptimize/*=false*/) const
		{

			if (isOptimize)
			{
				throw std::runtime_error("MatButterfly doesn't handle isOptimize flag");
			} else
			{
				return clone();
			}
		}


	template<typename FPP>
		faust_unsigned_int MatButterfly<FPP, GPU2>::getNonZeros() const
		{
			return d1.getNonZeros() + d2.getNonZeros();
		}


	template<typename FPP>
		void MatButterfly<FPP, GPU2>::transpose()
		{
			init_transpose(); // free cost if already called once
			is_transp = ! is_transp;
		}


	template<typename FPP>
		void MatButterfly<FPP, GPU2>::init_transpose()
		{
			//TODO: simplify in case of symmetric matrix (it happens for the FFT)
			if(d2t.size() == 0)
			{
				//TODO: do it all in GPU memory
				auto size = d2.size();
				FPP *d2_ptr, *d2t_ptr;
				auto cpu_d2 = d2.tocpu();
				d2_ptr = cpu_d2.getData();
				d2t.resize(size);
				Vect<FPP, Cpu> cpu_d2t(size);
				d2t_ptr = cpu_d2t.getData();

				auto d_offset = size >> (level+1);
				// D1 doesn't change
				// swap every pair of D2 contiguous blocks to form D2T
				for(int i = 0;i < size; i += d_offset * 2)
				{
					// swap two next blocks of size d_offset into d2t_ptr
					std::copy(d2_ptr + i, d2_ptr + i + d_offset, d2t_ptr + i + d_offset);
					std::copy(d2_ptr + i + d_offset, d2_ptr + i + 2 * d_offset, d2t_ptr + i);
				}
				d2t =  cpu_d2t;
			}
		}

	template<typename FPP>
		void MatButterfly<FPP, GPU2>::conjugate()
		{
			d1.conjugate();
			d2.conjugate();
			if(d2t.size() > 0)
				d2t.conjugate();
		}


	template<typename FPP>
		void MatButterfly<FPP, GPU2>::adjoint()
		{
			transpose();
			conjugate();
		}


	template<typename FPP>
		void MatButterfly<FPP, GPU2>::Display() const
		{
			//TODO: adjust consistently with MatGeneric Display (using to_string)
			std::cout << "MatButterfly on GPU: ";
			std::cout << "D1: ";
			d1.tocpu().Display();
			std::cout << "D2: ";
			d2.tocpu().Display();
			if(d2t.size() > 0)
			{
				std::cout << "D2T: ";
				d2t.tocpu().Display();
			}
			std::cout << "subdiag_ids: ";
			for(int i=0;i < d2.size();i++)
				std::cout << subdiag_ids[i] << " ";
			std::cout << std::endl;
		}

	template<typename FPP>
		MatDense<FPP, GPU2> MatButterfly<FPP, GPU2>::multiply(const FPP* X, int X_ncols)
		{
			MatDense<FPP, GPU2> gpu_X(d1.size(), X_ncols, X);
			return multiply(gpu_X);
		}

	template<typename FPP>
		MatDense<FPP, GPU2> MatButterfly<FPP, GPU2>::multiply(const FPP* x)
		{
			return multiply(x, 1);
		}


	template<typename FPP>
		void MatButterfly<FPP, GPU2>::multiply(MatDense<FPP, GPU2> &gpu_X, MatDense<FPP, Cpu> &cpu_out)
		{
			multiply(gpu_X);
			gpu_X.tocpu(cpu_out);
		}

	template<typename FPP>
		MatDense<FPP, GPU2> MatButterfly<FPP, GPU2>::multiply(MatDense<FPP, GPU2> &gpu_X)
		{
			//TODO: do we really need this method?
			butterfly_diag_prod(gpu_X, d1, d2, subdiag_ids);
			return gpu_X;
		}

	template<typename FPP>
		MatButterfly<FPP, GPU2>::~MatButterfly()
		{
			if(subdiag_ids)
				delete[] subdiag_ids;
		}
}
