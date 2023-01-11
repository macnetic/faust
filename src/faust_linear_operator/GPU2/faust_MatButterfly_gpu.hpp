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
	}

	template<typename FPP>
		void MatButterfly<FPP, GPU2>::setZeros()
		{
			throw std::runtime_error("setZeros is not available on a MatButterfly matrix.");
		}


	template<typename FPP>
		void MatButterfly<FPP, GPU2>::multiply(MatDense<FPP, GPU2> &other, const char op_this)
		{
			butterfly_diag_prod(other, d1, d2, subdiag_ids);
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
		MatSparse<FPP, GPU2> MatButterfly<FPP, GPU2>::toMatSparse()
		{
			MatSparse<FPP, GPU2> sp(this->getNbRow(), this->getNbCol());
			sp.setEyes();
			multiply(sp, 'N');
			return sp;
		}
}
