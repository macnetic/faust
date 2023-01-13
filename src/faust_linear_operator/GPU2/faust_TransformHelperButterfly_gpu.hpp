namespace Faust
{

	template<typename FPP>
		TransformHelperButterfly<FPP, GPU2>::TransformHelperButterfly(const std::vector<MatGeneric<FPP,Cpu> *>& facts, const FPP lambda_ /*= (FPP)1.0*/, const bool optimizedCopy/*=false*/, const bool cloning_fact /*= true*/, const bool internal_call/*=false*/)
		{
			int i = 0;
			auto size = facts[0]->getNbRow();
			//		for(auto csr_fac: facts)
			// use rather recorded factors in the Faust::Transform because one might have been multiplied with lambda_
			auto log2nf = 1 << (facts.size() - 1);
			has_permutation = (log2nf - size) == 0;
			auto end_it = has_permutation?this->end()-1:this->end();
			for(auto gen_fac: facts)
			{
				auto csr_fac = dynamic_cast<const MatSparse<FPP, Cpu>*>(gen_fac);
				if(csr_fac == nullptr)
					throw std::runtime_error("TransformHelperButterfly can receive only MatSparse CSR matrices");
				if(i < facts.size()-1 || ! has_permutation)
				{
					if( i == 0)
					{

						auto mul_csr = new MatSparse<FPP, Cpu>(*csr_fac);
						*mul_csr *= lambda_;
						opt_factors.insert(opt_factors.begin(),
								MatButterfly<FPP, GPU2>(*mul_csr, i++));
						this->push_back(mul_csr);
					}
					else
					{
						opt_factors.insert(opt_factors.begin(),
								MatButterfly<FPP, GPU2>(*csr_fac, i++));
						this->push_back(csr_fac);
					}
				}

			}
			if(has_permutation)
			{
				// set the permutation factor
				auto csr_fac = dynamic_cast<const MatSparse<FPP, Cpu>*>(*(facts.end()-1));
				P = MatPerm<FPP, GPU2>(*csr_fac);
			}
		}

	template<typename FPP>
		TransformHelperButterfly<FPP, GPU2>::TransformHelperButterfly(const TransformHelper<FPP, Cpu> & cputh)
		{
			MatButterfly<FPP, Cpu>* mbf;
			MatPerm<FPP, Cpu>* mp;
			MatSparse<FPP, Cpu>* ms;
			MatSparse<FPP, Cpu> sp;
			has_permutation = false;
			for(auto gen_fac: cputh)
			{
				if(mbf = dynamic_cast<MatButterfly<FPP, Cpu>*>(gen_fac))
				{
					opt_factors.insert(opt_factors.begin(), MatButterfly<FPP, GPU2>(*mbf));	
					sp = mbf->toMatSparse();
				}
				else if(mp = dynamic_cast<MatPerm<FPP, Cpu>*>(gen_fac))
				{
					P = MatPerm<FPP, GPU2>(*mp);	
					sp = mp->toMatSparse();
					has_permutation = true;
				}
				else if((ms = dynamic_cast<MatSparse<FPP, Cpu>*>(gen_fac)) && MatPerm<FPP, Cpu>::isPerm(*ms, false))
				{
					sp = *ms;
					P = MatPerm<FPP, GPU2>(sp);	
					has_permutation = true;
				}
				else
					throw std::runtime_error("Cannot convert CPU TransformHelper to GPU TransformHelperButterfly if it contains other matrix type than MatButterfly and MatPerm or MatSparse permutation");
				this->push_back(new MatSparse<FPP, GPU2>(sp));
			}
		}

	template<typename FPP>
		void TransformHelperButterfly<FPP, GPU2>::multiply(const FPP* A, int A_ncols, FPP* C)
	{

			MatDense<FPP, GPU2> gpu_X(this->getNbRow(), A_ncols, A);

			int i = 0;
			if(has_permutation)
				P.multiply(gpu_X, 'N');

			for(auto gpu_bmat: opt_factors)
				gpu_bmat.multiply(gpu_X);

			gpu_X.tocpu(C, nullptr);
	}


	template<typename FPP>
			Vect<FPP, Cpu> TransformHelperButterfly<FPP, GPU2>::multiply(const Vect<FPP, Cpu>& x)
			{
				Vect<FPP, Cpu> y;
				y.resize(this->getNbRow());
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
			y.resize(this->getNbRow());
			multiply(x, 1, y.getData());
			return y;
		}


	template<typename FPP>
		MatDense<FPP, Cpu> TransformHelperButterfly<FPP, GPU2>::multiply(const MatDense<FPP,Cpu> &A)
		{
			MatDense<FPP, Cpu> out;
			out.resize(this->getNbRow(), A.getNbCol());
			multiply(A.getData(), A.getNbCol(), out.getData());
			return out;
		}

	template<typename FPP>
		MatDense<FPP, Cpu> TransformHelperButterfly<FPP,GPU2>::multiply(const MatSparse<FPP,Cpu> &X)
		{
			return multiply(MatDense<FPP, Cpu>(X));
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




}
