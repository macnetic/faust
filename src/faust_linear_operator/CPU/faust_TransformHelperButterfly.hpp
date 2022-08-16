namespace Faust
{

	template<typename FPP>
		TransformHelperButterfly<FPP,Cpu>::TransformHelperButterfly(const std::vector<MatGeneric<FPP,Cpu> *>& facts, const FPP lambda_ /*= (FPP)1.0*/, const bool optimizedCopy/*=false*/, const bool cloning_fact /*= true*/, const bool internal_call/*=false*/) : TransformHelper<FPP, Cpu>(facts, lambda_, optimizedCopy, cloning_fact, internal_call)
	{
		int i = 0;
		for(auto csr_fac: facts)
			if(i < facts.size()-1)
				opt_factors.insert(opt_factors.begin(), ButterflyMat<FPP>(*dynamic_cast<const MatSparse<FPP, Cpu>*>(csr_fac), i++));
		// set the permutation factor
		auto csr_fac = *(facts.end()-1);
		auto size = csr_fac->getNbRow();
		perm_d = Vect<FPP, Cpu>(size);
		perm_d_ptr = perm_d.getData();
		//TODO: only a setOnes should be enough, this is a permutation matrix
		memcpy(perm_d_ptr, dynamic_cast<const MatSparse<FPP, Cpu>*>(csr_fac)->getValuePtr(), size*sizeof(FPP));
		auto bitrev_perm_ids = new unsigned int[size];
		iota(bitrev_perm_ids, bitrev_perm_ids+size, 0);
		bit_rev_permu(facts.size()-1, bitrev_perm_ids);
		bitrev_perm.resize(size);
		copy(bitrev_perm_ids, bitrev_perm_ids+size, bitrev_perm.begin());
		delete[] bitrev_perm_ids;
	}

	template<typename FPP>
		TransformHelperButterfly<FPP,Cpu>* TransformHelperButterfly<FPP,Cpu>::fourierFaust(unsigned int n, const bool norma)
		{

			std::vector<MatGeneric<FPP,Cpu>*> factors(n+1);
			TransformHelperButterfly<FPP,Cpu>* fourierFaust = nullptr;
			try
			{
				fft_factors(n, factors);
				FPP alpha = norma?FPP(1/sqrt((double)(1 << n))):FPP(1.0);
				fourierFaust = new TransformHelperButterfly<FPP, Cpu>(factors, alpha, false, false, /* internal call */ true);
			}
			catch(std::bad_alloc e)
			{
				//nothing to do, out of memory, return nullptr
			}
			return fourierFaust;
		}


	template<typename FPP>
		Vect<FPP, Cpu> TransformHelperButterfly<FPP,Cpu>::multiply(const Vect<FPP, Cpu>& x)
		{
			Vect<FPP, Cpu> y(this->getNbRow());
//			for(int i=0;i < this->getNbRow(); i++)
//				y.getData()[i] = perm_d_ptr[i] * x.getData()[bitrev_perm[i]];
//
//			int i = 0;
//			for(auto fac: opt_factors)
//				y = fac.multiply(y);
			multiply(x.getData(), y.getData());
			return y;
		}

	template<typename FPP>
		void TransformHelperButterfly<FPP,Cpu>::multiply(const FPP* x, FPP* y)
		{
			auto size = this->getNbRow();
			if(x == y)
			{
				// an intermediate vector is needed to index x
				auto x_ = new FPP[size];
				for(int i=0;i < size; i++)
					x_[i] = x[bitrev_perm[i]];
				for(int i=0;i < size; i++)
					y[i] = perm_d_ptr[i] * x_[i];
				delete[] x_;
			}
			else
				for(int i=0;i < this->getNbRow(); i++)
					y[i] = perm_d.getData()[i] * x[bitrev_perm[i]];

			Vect<FPP, Cpu> z(size);
			int i = 0;
			for(auto fac: opt_factors)
			{
				if(i & 1)
					fac.multiply(z.getData(), y, this->getNbRow());
				else
					fac.multiply(y, z.getData(), this->getNbRow());
				i++;
			}
			if(i & 1)
				memcpy(y, z.getData(), size*sizeof(FPP));
		}

	template<typename FPP>
		Vect<FPP, Cpu> TransformHelperButterfly<FPP,Cpu>::multiply(const FPP* x)
		{
			Vect<FPP, Cpu> y(this->getNbRow());
			multiply(x, y.getData());
			return y;
		}
}


namespace Faust
{
	template<typename FPP>
	ButterflyMat<FPP>::ButterflyMat(const MatSparse<FPP, Cpu> &factor, int level)
	{
		// build a d1, d2 pair from the butterfly factor
		auto size = factor.getNbRow();
		d1 = Vect<FPP, Cpu>(size);
		d2 = Vect<FPP, Cpu>(size);
		FPP *d1_ptr, *d2_ptr;
		d1_ptr = d1.getData();
		d2_ptr = d2.getData();
		auto d_offset = size >> (level+1);
		auto data = factor.getValuePtr();
		auto rowptr = factor.getRowPtr();
		for(int i=0;i < size; i++)
		{
			if((i / d_offset) % 2)
			{
				// d2 coeff is the first elt of row i
				d2_ptr[i] = data[rowptr[i]];
				d1_ptr[i] = data[rowptr[i]+1]; // diag elt is just after
			}
			else
			{
				// d2 coeff is the last elt of row i
				d2_ptr[i] = data[rowptr[i+1]-1];
				d1_ptr[i] = data[rowptr[i+1]-2]; // diag elt is just before
			}
		}
		std::vector<int> seq(size);
		iota(seq.begin(), seq.end(), 0);
		subdiag_ids.resize(size);
		for(int i = 0;i < size; i += d_offset * 2)
		{
			copy(seq.begin()+i+d_offset, seq.begin()+i+2*d_offset, subdiag_ids.begin()+i);
			copy(seq.begin()+i, seq.begin()+i+d_offset, subdiag_ids.begin()+i+d_offset);
		}
		this->level = level;
	}

	template<typename FPP>
	void ButterflyMat<FPP>::Display() const
	{
		std::cout << "d1: ";
		d1.Display();
		std::cout << "d2: ";
		d2.Display();
		cout << "subdiag_ids: ";
		for(int i=0;i < subdiag_ids.size();i++)
			cout << subdiag_ids[i] << " ";
		cout << std::endl;
	}

	template<typename FPP>
	Vect<FPP, Cpu> ButterflyMat<FPP>::multiply(const Vect<FPP, Cpu>& x) const
	{
		Vect<FPP, Cpu> z(x.size());
//		const FPP *x_ptr = x.getData(), *d1_ptr = d1.getData(), *d2_ptr = d2.getData();
//		for(int i=0;i < x.size(); i++)
//		{
//			z[i] = d1_ptr[i] * x[i] + d2_ptr[i] * x_ptr[subdiag_ids[i]];
//		}
		multiply(x.getData(), z.getData(), x.size());
		return z;
	}

	template<typename FPP>
	void ButterflyMat<FPP>::multiply(const FPP* x, FPP* y, size_t size) const
	{
		const FPP *d1_ptr = d1.getData(), *d2_ptr = d2.getData();
		for(int i=0;i < size; i++)
			y[i] = d1_ptr[i] * x[i] + d2_ptr[i] * x[subdiag_ids[i]];
	}
}
