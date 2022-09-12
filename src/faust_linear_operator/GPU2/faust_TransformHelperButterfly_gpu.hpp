#include "faust_TransformHelperButterfly.h"
namespace Faust
{

	template<typename FPP>
		TransformHelperButterfly<FPP, GPU2>::TransformHelperButterfly(const std::vector<MatGeneric<FPP,Cpu> *>& facts, const FPP lambda_ /*= (FPP)1.0*/, const bool optimizedCopy/*=false*/, const bool cloning_fact /*= true*/, const bool internal_call/*=false*/)
		{
			int i = 0;
			auto size = this->getNbRow();
			//		for(auto csr_fac: facts)
			// use rather recorded factors in the Faust::Transform because one might have been multiplied with lambda_
			auto log2nf = 1 << (this->size() - 1);
			has_permutation = (log2nf - this->getNbRow()) == 0;
			auto end_it = has_permutation?this->end()-1:this->end();
			for(auto csr_fac_it = this->begin(); csr_fac_it != end_it; csr_fac_it++)
			{
				auto csr_fac = *csr_fac_it;
				opt_factors.insert(opt_factors.begin(), ButterflyMat<FPP, GPU2>(*dynamic_cast<const MatSparse<FPP, Cpu>*>(csr_fac), i++));
				this->push_back(csr_fac);
			}
			if(has_permutation)
			{
				// set the permutation factor
				auto csr_fac = dynamic_cast<const MatSparse<FPP, Cpu>*>(*(this->end()-1));
				this->push_back(csr_fac);
				d_perm.resize(size);
				// only ones should be enough because this is a permutation matrix but it could be normalized
				d_perm = Vect<FPP, GPU2>(size, csr_fac->getValuePtr());
				perm_ids = new int[size];
				copy(csr_fac->getColInd(), csr_fac->getColInd()+size, perm_ids);
			}
		}

	template<typename FPP>
		void TransformHelperButterfly<FPP, GPU2>::multiply(const FPP* A, int A_ncols, FPP* C)
	{

			MatDense<FPP, GPU2> gpu_X(d_perm.size(), A_ncols, A);

			if(has_permutation)
				gpu_X.eltwise_mul(d_perm, perm_ids);

			for(auto gpu_bmat: opt_factors)
				gpu_bmat.multiply(gpu_X);

			gpu_X.tocpu(C);
	}


	template<typename FPP>
			Vect<FPP, Cpu> TransformHelperButterfly<FPP, GPU2>::multiply(const Vect<FPP, Cpu>& x)
			{
				Vect<FPP, Cpu> y;
				y.resize(d_perm.size());
				multiply(x.getData(), y.getData());
				return y;
			}


	template<typename FPP>
		void TransformHelperButterfly<FPP, GPU2>::multiply(const FPP* x, FPP* y)
		{
			multiply(x, 1, y);
		}

	template<typename FPP>
		Vect<FPP, Cpu> TransformHelperButterfly<FPP, GPU2>::multiply(const FPP* x)
		{
			Vect<FPP, Cpu> y;
			y.resize(d_perm.size());
			multiply(x, 1, y.getData());
			return y;
		}


	template<typename FPP>
		MatDense<FPP, Cpu> TransformHelperButterfly<FPP, GPU2>::multiply(const MatDense<FPP,Cpu> &A)
		{
			MatDense<FPP, Cpu> out;
			out.resize(d_perm.size(), A.getNbCol());
			multiply(A.getData(), A.getNbCol(), out.getData());
			return out;
		}

	template<typename FPP>
		MatDense<FPP, Cpu> TransformHelperButterfly<FPP,GPU2>::multiply(const MatSparse<FPP,Cpu> &X)
		{
			return multiply(MatDense<FPP, GPU2>(X));
		}

	template<typename FPP>
	TransformHelper<FPP,GPU2>* TransformHelperButterfly<FPP,GPU2>::fourierFaust(unsigned int n, const bool norma/*=true*/)
	{
		std::vector<MatGeneric<FPP,Cpu>*> factors(n+1);
		TransformHelper<FPP, GPU2>* fourierFaust = nullptr;
		try
		{
			fft_factors(n, factors);
			FPP alpha = norma?FPP(1/sqrt((double)(1 << n))):FPP(1.0);
			fourierFaust = new TransformHelperButterfly<FPP, GPU2>(factors, alpha, false, false, /* internal call */ true);
		}
		catch(std::bad_alloc e)
		{
			//nothing to do, out of memory, return nullptr
		}
		return fourierFaust;
	}

	template<typename FPP>
		ButterflyMat<FPP, GPU2>::ButterflyMat(const MatSparse<FPP, Cpu> &factor, int level) : level(level)
	{
		ButterflyMat<FPP, Cpu> cpu_bmat(factor, level);
		auto cpu_d1 = cpu_bmat.getD1();
		auto cpu_d2 = cpu_bmat.getD2();
		d1 = Vect<FPP, GPU2>(cpu_d1.size(), cpu_d1.diagonal().data());
		d2 = Vect<FPP, GPU2>(cpu_d2.size(), cpu_d2.diagonal().data());

		auto sd_ids_vec = cpu_bmat.get_subdiag_ids();
		subdiag_ids = new int[sd_ids_vec.size()];
		memcpy(subdiag_ids, sd_ids_vec.data(), sizeof(int) * sd_ids_vec.size());
	}


	template<typename FPP>
	void ButterflyMat<FPP, GPU2>::Display() const
	{
		std::cout << "ButterflyMat on GPU: ";
		std::cout << "D1: ";
		d1.Display();
		std::cout << "D2: ";
		d1.Display();
		cout << "subdiag_ids: ";
		for(int i=0;i < d1.size();i++)
			cout << subdiag_ids[i] << " ";
		cout << std::endl;
	}

	template<typename FPP>
		MatDense<FPP, GPU2> ButterflyMat<FPP, GPU2>::multiply(const FPP* X, int X_ncols)
		{
			MatDense<FPP, GPU2> gpu_X(d1.size(), X_ncols, X);
			return multiply(gpu_X);
		}

	template<typename FPP>
		MatDense<FPP, GPU2> ButterflyMat<FPP, GPU2>::multiply(const FPP* x)
		{
			return multiply(x, 1);
		}

	template<typename FPP>
		MatDense<FPP, GPU2> ButterflyMat<FPP, GPU2>::multiply(MatDense<FPP, GPU2> &gpu_X)
		{
			MatDense<FPP, GPU2> gpu_X2(gpu_X);
			gpu_X.eltwise_mul(d2, subdiag_ids);
			gpu_X2.eltwise_mul(d1);
			gpu_X += gpu_X2;
			return gpu_X;
		}

	template<typename FPP>
		void ButterflyMat<FPP, GPU2>::multiply(MatDense<FPP, GPU2> &gpu_X, MatDense<FPP, Cpu> &cpu_out)
		{
			multiply(gpu_X);
			gpu_X.tocpu(cpu_out);
		}
}
