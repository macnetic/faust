namespace Faust
{
	template<typename FPP,FDevice DEVICE> class Transform;

	template<typename FPP>
		TransformHelper<FPP,GPU2>::TransformHelper() : TransformHelperGen<FPP,GPU2>()
	{
	}

	template<typename FPP>
		TransformHelper<FPP,GPU2>::TransformHelper(const TransformHelper<FPP,GPU2>& th, bool transpose, bool conjugate) : TransformHelperGen<FPP,GPU2>(&th, transpose, conjugate)
	{
	}

	template<typename FPP>
		TransformHelper<FPP,GPU2>::TransformHelper(const std::vector<MatGeneric<FPP,GPU2> *>& facts, const FPP lambda_/*= (FPP)1.0*/, const bool optimizedCopy/*=false*/, const bool cloning_fact /*= true*/, const bool internal_call/*=false*/) : TransformHelper<FPP,GPU2>()
	{

		//if lambda is not 1.0 a factor will be multiplied and so it needs to be copied to preserve the original that could be used elsewhere
		// in an optimization purpose, the smallest factor is copied
		int min_size_id = 0;
		if(lambda_ != FPP(1.0))
		{
			std::vector<int> fact_ids(facts.size());
			int i = -1;
			std::generate(fact_ids.begin(), fact_ids.end(), [&i](){return ++i;});
			std::vector<int>::iterator result = std::min_element(fact_ids.begin(), fact_ids.end(), [&facts](const int &a, const int &b){return facts[a]->getNBytes() < facts[b]->getNBytes();});
			min_size_id = std::distance(fact_ids.begin(), result);
		}
		for(int i=0; i < facts.size(); i++)
		{
			if(i == min_size_id)
				this->push_back(facts[min_size_id], false, cloning_fact || lambda_ != (FPP) 1.0);
			else
				this->push_back(facts[i], false, cloning_fact);
		}
		this->transform->multiply(lambda_, min_size_id);
	}

	template<typename FPP>
		TransformHelper<FPP,GPU2>::TransformHelper(const TransformHelper<FPP,Cpu>& cpu_t, int32_t dev_id/*=-1*/, void* stream/*=nullptr*/) : TransformHelper<FPP,GPU2>()
		{
			for(auto cpu_fact: cpu_t)
				this->push_back(cpu_fact, false, dev_id, stream);
            this->is_transposed = cpu_t.is_transposed;
            this->is_conjugate = cpu_t.is_conjugate;
            //TODO: slice etc.
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>::TransformHelper(TransformHelper<FPP,GPU2>* th, Slice s[2]): TransformHelper<FPP,GPU2>()
	{
		this->init_sliced_transform(th, s);
	}

	template<typename FPP>
		TransformHelper<FPP,GPU2>::TransformHelper(TransformHelper<FPP,GPU2>* th, faust_unsigned_int* row_ids, faust_unsigned_int num_rows, faust_unsigned_int* col_ids, faust_unsigned_int num_cols): TransformHelper<FPP,GPU2>()
	{
		this->init_fancy_idx_transform(th, row_ids, num_rows, col_ids, num_cols);
	}

#ifndef IGNORE_TRANSFORM_HELPER_VARIADIC_TPL
	template<typename FPP>
		template<typename ... GList>
		TransformHelper<FPP,GPU2>::TransformHelper(GList& ... t): TransformHelper<FPP,GPU2>()
		{
			this->push_back_(t...);
		}
#endif



	template<typename FPP>
		void TransformHelper<FPP,GPU2>::push_back(const MatGeneric<FPP,GPU2>* M, const bool optimizedCopy/*=false*/, const bool copying/*=true*/, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
			//optimizedCopy is ignored because not handled yet by Transform<FPP,GPU2> // TODO ? (it's not used by wrappers anyway)
			this->transform->push_back(M, copying, transpose, conjugate);
		}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::push_back(const MatGeneric<FPP,Cpu>* M, const bool optimizedCopy/*=false*/, const int32_t dev_id/*=-1*/, const void* stream/*=nullptr*/)
		{

			MatGeneric<FPP,GPU2>* gpu_M = nullptr;
			const MatDense<FPP,Cpu>* cpu_dM = nullptr;
			const MatSparse<FPP,Cpu>* cpu_sM = nullptr;
			if(nullptr != (cpu_dM = dynamic_cast<const MatDense<FPP,Cpu>*>(M)))
			{
				auto gpu_dM = new MatDense<FPP,GPU2>(M->getNbRow(), M->getNbCol(), cpu_dM->getData());
				gpu_M = gpu_dM;
			}
			else if(nullptr != (cpu_sM = dynamic_cast<const MatSparse<FPP,Cpu>*>(M)))
			{
				auto gpu_sM = new MatSparse<FPP,GPU2>(M->getNbRow(), M->getNbCol(), cpu_sM->getNonZeros(), cpu_sM->getValuePtr(), cpu_sM->getRowPtr(), cpu_sM->getColInd(), dev_id);
				gpu_M = gpu_sM;
			}
			this->transform->push_back(gpu_M, false);
		}


	template<typename FPP>
		void TransformHelper<FPP,GPU2>::push_back(const FPP* data, const int* row_ptr, const int* id_col, const int nnz, const int nrows, const int ncols, const bool optimizedCopy/*=false*/, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
			auto sparse_mat = new MatSparse<FPP,GPU2>(nrows, ncols, nnz, data, row_ptr, id_col);
			this->push_back(sparse_mat, optimizedCopy, false, transpose, conjugate);
//			if(optimizedCopy) delete sparse_mat; // optimizedCopy not supported on GPU2
		}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::push_first(const MatGeneric<FPP,GPU2>* M, const bool optimizedCopy/*=false*/, const bool copying/*=true*/)
		{
			return this->transform->push_first(M, copying);
		}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::display() const
		{
			this->transform->Display(this->is_transposed);
		}

	template<typename FPP>
		std::string TransformHelper<FPP,GPU2>::to_string() const
		{
			return this->transform->to_string(this->is_transposed);
		}

	template<typename FPP>
		template<typename Head, typename ... Tail>
		void TransformHelper<FPP,GPU2>::push_back_(Head& h, Tail&... t)
		{
			for(auto it=h.begin(); it < h.end(); it++)
			{
				auto f = *it;
				this->push_back(f, false, false);
			}
			this->push_back_(t...);
		}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::push_back_()
		{
			// do nothing, here just for empty tail of above function
		}

	template<typename FPP>
		MatDense<FPP,GPU2> TransformHelper<FPP,GPU2>::get_product(int prod_mod/*=-1*/)
		{
			return this->transform->get_product(this->isTransposed2char(), this->is_conjugate);
		}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::get_product(MatDense<FPP,GPU2>& M, int prod_mod/*=-1*/)
		{
			this->transform->get_product(M, this->isTransposed2char(), this->is_conjugate);
		}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::get_product(MatDense<FPP,Cpu>& M, int prod_mod/*=-1*/)
		{
			MatDense<FPP,GPU2> gpuM;
			this->get_product(gpuM, prod_mod);
			M = gpuM.tocpu();
		}

	template<typename FPP>
		Real<FPP> TransformHelper<FPP,GPU2>::normFro() const
		{
			return this->transform->get_product().norm();
		}

	template<typename FPP>
		Real<FPP> TransformHelper<FPP,GPU2>::normL1() const
		{
			return this->transform->get_product().normL1();
		}

	template<typename FPP>
		Real<FPP> TransformHelper<FPP,GPU2>::normInf() const
		{
			return this->transform->normL1(!this->is_transposed);
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,GPU2>::size() const
		{
			return this->transform->size();
		}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::update_total_nnz() const
		{
			this->transform->update_total_nnz();
		}

	template<typename FPP>
		Real<FPP> TransformHelper<FPP,GPU2>::spectralNorm(int32_t nb_iter_max, float threshold, int& flag)
		{
			return this->transform->spectralNorm(nb_iter_max, threshold, flag);
		}

	template<typename FPP>
		const MatGeneric<FPP,GPU2>* TransformHelper<FPP,GPU2>::get_gen_fact(const faust_unsigned_int id) const
		{
			return this->transform->get_fact(id, false);
		}

	template<typename FPP>
		MatGeneric<FPP,GPU2>* TransformHelper<FPP,GPU2>::get_gen_fact_nonconst(const faust_unsigned_int id) const
		{
			return this->transform->get_fact(id, false);
		}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::update(const MatGeneric<FPP, GPU2>& M,const faust_unsigned_int id)
		{
			return this->transform->update(M, id);
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::multiply(const TransformHelper<FPP,GPU2>* right)
		{
			// The goal is to minimize the number of factors copied (but maybe the criterion should be the sum of the size of these factors rather than their number)
//			std::cout << "===" << this->is_transposed << std::endl;
//			this->display();
//			std::cout << "first fact:" << std::endl;
//			this->get_gen_fact(0)->Display();
//			std::cout << "===" << right->is_transposed << std::endl;
//			right->display();
//			std::cout << "===" << std::endl;
			bool copying_this = false;
			bool copying_right = false;
			bool conj_diff = this->is_conjugate != right->is_conjugate;
			bool trans_diff = this->is_transposed != right->is_transposed;
			bool transconj_diff = conj_diff || trans_diff;
			bool out_transp = false, out_conj = false;
			bool transp_this = false, transp_right = false, conj_this = false, conj_right = false;
			if(transconj_diff)
				if(this->size() < right->size())
				{
					copying_this = true;
					out_transp = trans_diff && right->is_transposed;
					out_conj = conj_diff && right->is_conjugate;
					transp_this = trans_diff;
					conj_this = conj_diff;
				}
				else
				{
					copying_right = true;
					out_transp = trans_diff && this->is_transposed;
					out_conj = conj_diff && this->is_conjugate;
					transp_right = trans_diff;
					conj_right = conj_diff;
				}
			auto mul_faust = new TransformHelper<FPP,GPU2>();
//			std::cout << "transp_this: " << transp_this << " conj_this: " << conj_this << std::endl;
//			std::cout << "transp_right: " << transp_right << " conj_right: " << conj_right << std::endl;
//			std::cout << "out_trans: " << out_transp << " out_conj: " << out_conj << std::endl;
			std::function<void()> push_right_facts = [&out_transp, &transp_right, &mul_faust, &right, &copying_right, &conj_right]()
			{
				if(transp_right)
					for(int i=right->size()-1; i>= 0; i--)
						mul_faust->push_back(right->get_gen_fact(i), /*OptimizedCopy*/ false, copying_right, /*transp_right*/ true, conj_right);
				else
					for(auto f: *right)
						mul_faust->push_back(f, /*OptimizedCopy*/ false, copying_right, /*transp_right*/ false, conj_right);
			};
			std::function<void()> push_this_facts = [&transp_this, &mul_faust, this, &copying_this, &conj_this]()
			{
				if(transp_this)
					for(int i=size()-1; i>= 0; i--)
						mul_faust->push_back(get_gen_fact(i), /*OptimizedCopy*/ false, copying_this, /*transp_this*/ true, conj_this);
				else
					for(auto f: *this)
						mul_faust->push_back(f, /*OptimizedCopy*/ false, copying_this, /*transp_this*/ false, conj_this);
			};
			if(out_transp)
			{
				push_right_facts();
				push_this_facts();
			}
			else
			{
				push_this_facts();
				push_right_facts();
			}
			mul_faust->is_transposed = out_transp;
			mul_faust->is_conjugate = out_conj;
			return mul_faust;
		}

	template<typename FPP>
		MatDense<FPP,GPU2> TransformHelper<FPP,GPU2>::multiply(const Faust::MatDense<FPP,GPU2> &A, const bool transpose /* deft to false */, const bool conjugate/*=false*/)
		{
			this->is_transposed ^= transpose;
			this->is_conjugate ^= conjugate;
			MatDense<FPP,GPU2> M = this->transform->multiply(A, this->isTransposed2char());
			this->is_transposed ^= transpose;
			this->is_conjugate ^= conjugate;
			return M;
		}

	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,GPU2>::multiply(const Faust::MatDense<FPP,Cpu> &A, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
			this->is_transposed ^= transpose;
			this->is_conjugate ^= conjugate;
			MatDense<FPP,GPU2> M = this->multiply(MatDense<FPP,GPU2>(A), transpose, conjugate);
			this->is_conjugate ^= conjugate;
			this->is_transposed ^= transpose;
			return M.tocpu();
		}

	template<typename FPP>
		Vect<FPP,Cpu> TransformHelper<FPP,GPU2>::multiply(const Faust::Vect<FPP,Cpu> &A, const bool transpose /* deft to false */, const bool conjugate/*=false*/)
		{
			this->is_transposed ^= transpose;
			this->is_conjugate ^= conjugate;
			Vect<FPP,GPU2> gpu_A(A.size(), A.getData());
			Vect<FPP,GPU2> v = this->multiply(gpu_A , transpose, conjugate); //TODO: handle transpose and conjugate
			this->is_transposed ^= transpose;
			this->is_conjugate ^= conjugate;
			return v.tocpu();
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::multiply(const FPP& a)
		{
			const vector<MatGeneric<FPP,GPU2>*>& vec = this->transform->data; //TransformHelper is a friend class of Transform // we can access private attribute data
			TransformHelper<FPP,GPU2>* th = new TransformHelper<FPP,GPU2>(vec, a, false, false, true);
			th->copy_transconj_state(*this);
			th->copy_slice_state(*this);
			return th;
		}

	template<typename FPP>
		Vect<FPP,GPU2> TransformHelper<FPP,GPU2>::multiply(const Faust::Vect<FPP,GPU2>& a, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
			this->is_transposed ^= transpose;
			this->is_conjugate ^= conjugate;
			Vect<FPP,GPU2> v = this->transform->multiply(a, this->isTransposed2char());
			this->is_transposed ^= transpose;
			this->is_conjugate ^= conjugate;
			return v;
		}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::set_FM_mul_mode() const
		{
			throw std::runtime_error("set_FM_mul_mode is yet to implement in Faust C++ core for GPU.");
		}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::set_Fv_mul_mode() const
		{
			throw std::runtime_error("set_Fv_mul_mode is yet to implement in Faust C++ core for GPU.");
		}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::pop_front()
		{
			return this->transform->pop_front();
		}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::pop_back()
		{
			return this->transform->pop_back();
		}

	template<typename FPP>
		void Faust::TransformHelper<FPP,GPU2>::pack_factors(const int mul_order_opt_mode/*=DEFAULT*/)
		{
			TransformHelperGen<FPP,Cpu>::pack_factors(mul_order_opt_mode);
		}

	template<typename FPP>
		void Faust::TransformHelper<FPP,GPU2>::pack_factors(faust_unsigned_int start_id, faust_unsigned_int end_id,const int mul_order_opt_mode/*=DEFAULT*/)
		{
			if(start_id < 0 || start_id >= size())
				throw out_of_range("start_id is out of range.");
			if(end_id < 0 || end_id >= size())
				throw out_of_range("end_id is out of range.");
			Faust::MatDense<FPP,GPU2> * packed_fac = nullptr;
			if(end_id == start_id)
			{
				//nothing to do except converting to MatDense if start_id
				//factor is a MatSparse
				packed_fac = dynamic_cast<Faust::MatDense<FPP,GPU2>*>(*(begin()+start_id));
				if(packed_fac == nullptr)
				{// factor start_id is not at MatDense, convert it
					packed_fac = new MatDense<FPP,GPU2>(*dynamic_cast<Faust::MatSparse<FPP,GPU2>*>(*(begin()+start_id)));
				}
				else
				{
					return; //no change
				}
			}
			else
			{
				// we have to multiply factors from start_id to end_id into one matrix
				// simple way to do, 1) create a overhead-free TransformHelper with these factors
				// 2) call get_product() to override the start_id factors with the result on the end
				// 3) erase factors from start_id to end_id and insert packed factor too to replace them (that's Transform object responsibility).
				// 1)
				std::vector<Faust::MatGeneric<FPP,GPU2>*> topack_factors;
				for(int i=start_id;i <= end_id; i++)
					topack_factors.push_back(get_gen_fact_nonconst(i));
				//				std::vector<Faust::MatGeneric<FPP,GPU2>*> topack_factors(begin()+start_id, begin()+end_id+1);
				Faust::TransformHelper<FPP,GPU2> t(topack_factors, 1.0, false, false, false);
				//TODO: not yet implemented for GPU2
//				t.set_FM_mul_mode(mul_order_opt_mode);
				// 2)
				packed_fac = new MatDense<FPP,GPU2>(t.get_product());
			}
			// 3)
			faust_unsigned_int i = end_id;
			while(i>=start_id)
			{
				this->transform->erase(i);
				if(i == 0) break;
				i--;
			}
			this->transform->insert(start_id, packed_fac, false);
		}

	template<typename FPP>
		typename Transform<FPP,GPU2>::iterator TransformHelper<FPP,GPU2>::begin() const
		{
			return this->transform->begin();
		}

	template<typename FPP>
		typename Transform<FPP,GPU2>::iterator TransformHelper<FPP,GPU2>::end() const
		{
			return this->transform->end();
		}


	template<typename FPP>
		void TransformHelper<FPP,GPU2>::operator=(TransformHelper<FPP,GPU2>& th)
		{
			copy_state(th);
		}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::save_mat_file(const char* filepath) const
		{
			this->transform->save_mat_file(filepath, this->is_transposed, this->is_conjugate);
		}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::tocpu(TransformHelper<FPP,Cpu>& cpu_transf) const
		{
            //TODO: tocpu support of arguments transpose and conjugate
			auto t = this->transform->tocpu();
			for(auto fac: t)
            {
				cpu_transf.push_back(fac, false, false);
            }
            cpu_transf.is_transposed = this->is_transposed;
            cpu_transf.is_conjugate = this->is_conjugate;
            //TODO: slice etc.
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,GPU2>::tocpu() const
		{
			auto cpu_t = new TransformHelper<FPP,Cpu>();
			tocpu(*cpu_t);
			return cpu_t;
		}


	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,GPU2>::get_total_nnz() const
		{
			return this->transform->get_total_nnz();
		}


	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::normalize(const int meth /* 1 for 1-norm, 2 for 2-norm (2-norm), -1 for inf-norm */) const
		{
			TransformHelper<FPP,Cpu> th;
			this->tocpu(th);
			auto thn = th.normalize(meth);
			auto gpu_thn = new TransformHelper<FPP,GPU2>(*thn, -1, nullptr);
			delete thn;
			return gpu_thn;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::transpose()
		{

			auto t = new TransformHelper<FPP,GPU2>(*this, true, false);
            return t;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::conjugate()
		{
			auto t = new TransformHelper<FPP,GPU2>(*this, false, true);
            return t;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::adjoint()
		{
			auto t = new TransformHelper<FPP,GPU2>(*this, true, true);
            return t;
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,GPU2>::getNBytes() const
		{
			faust_unsigned_int nbytes = 0;
			for(auto fac : *this)
			{
				if(dynamic_cast<Faust::MatDense<FPP, GPU2>*>(fac))
					nbytes += fac->getNbCol() * fac->getNbRow() * sizeof(FPP);
				else if (dynamic_cast<Faust::MatSparse<FPP, GPU2>*>(fac))
					nbytes += fac->getNonZeros() * (sizeof(FPP) + sizeof(int)) + (fac->getNbRow() + 1) * sizeof(int); // by default storage index is int
				else
					throw runtime_error("Unknown matrix type.");
			}
			return nbytes;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::optimize_time(const bool transp/*=false*/, const bool inplace/*=false*/, const int nsamples/*=1*/)
		{
			throw std::runtime_error("optimize_time is yet to implement in Faust C++ core for GPU.");
			return nullptr;
//			TransformHelper<FPP,Cpu> th;
//			this->tocpu(th);
//			auto thn = th.optimize_time(transp, /*inplace*/ true, nsamples);
//			auto gpu_thn = new TransformHelper<FPP,GPU2>(*thn, -1, nullptr);
//			delete thn;
//			return gpu_thn;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::optimize(const bool transp/*=false*/)
		{
			throw std::runtime_error("optimize is yet to implement in Faust C++ core for GPU.");
			return nullptr;
//			TransformHelper<FPP,Cpu> th;
//			this->tocpu(th);
//			auto thn = th.optimize(transp);
//			auto gpu_thn = new TransformHelper<FPP,GPU2>(*thn, -1, nullptr);
//			delete thn;
//			return gpu_thn;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::fourierFaust(unsigned int n, const bool norma/*=true*/)
		{
			auto cpu_faust = TransformHelper<FPP,Cpu>::fourierFaust(n, norma);
			TransformHelper<FPP,GPU2>* gpu_faust = new TransformHelper<FPP,GPU2>(*cpu_faust, -1, nullptr /*TODO: dev_id and stream ?*/);
			delete cpu_faust;
			return gpu_faust;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::pruneout(const int nnz_tres, const int npasses/*=-1*/, const bool only_forward/*=false*/)
		{
			TransformHelper<FPP,Cpu> th;
			this->tocpu(th);
			auto thn = th.pruneout(nnz_tres, npasses, only_forward);
			auto gpu_thn = new TransformHelper<FPP,GPU2>(*thn, -1, nullptr);
			delete thn;
			return gpu_thn;
		}


	template<typename FPP>
	TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::clone(int32_t dev_id/*=-1*/, void* stream/*=nullptr*/)
	{
		//TODO: take the requested device into account
		return this->clone();
	}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::randFaust(int faust_nrows, int faust_ncols, RandFaustType t, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int min_dim_size, unsigned int max_dim_size, float density, bool per_row)
		{
			auto cpu_faust = TransformHelper<FPP,Cpu>::randFaust(faust_nrows, faust_ncols, t, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row);
//	TransformHelper<FPP,GPU2>::TransformHelper(const std::vector<MatGeneric<FPP,GPU2> *>& facts, const FPP lambda_/*= (FPP)1.0*/, const bool optimizedCopy/*=false*/, const bool cloning_fact /*= true*/, const bool internal_call/*=false*/)
			TransformHelper<FPP,GPU2>* gpu_faust = new TransformHelper<FPP,GPU2>(*cpu_faust/*TODO: dev_id and stream ?*/);
//		void TransformHelper<FPP,GPU2>::push_back(const MatGeneric<FPP,Cpu>* M, const bool optimizedCopy/*=false*/, const int32_t dev_id/*=-1*/, const void* stream/*=nullptr*/)
//			for(auto cpu_fact: *cpu_faust)
//				gpu_faust->push_back(cpu_fact, false, -1/*TODO: replace dev_id by an arg passed to the func */, nullptr);

			delete cpu_faust;
			return gpu_faust;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::randFaust(RandFaustType t, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int min_dim_size, unsigned int max_dim_size, float density, bool per_row)
		{
			return TransformHelper<FPP,GPU2>::randFaust(-1, -1, t, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row);
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::hadamardFaust(unsigned int n, const bool norma/*=true*/)
		{
//			std::cerr << "Warning: GPU2 hadamardFaust is implemented by copying the Faust on CPU RAM and copying them back." << std::endl;
			auto cpu_faust = TransformHelper<FPP,Cpu>::hadamardFaust(n, norma);
			TransformHelper<FPP,GPU2>* gpu_faust = new TransformHelper<FPP,GPU2>(*cpu_faust, -1, nullptr /*TODO: dev_id and stream ?*/);
			delete cpu_faust;
			return gpu_faust;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::eyeFaust(unsigned int n, unsigned int m)
		{
//			std::cerr << "Warning: GPU2 eyeFaust is implemented by copying the Faust on CPU RAM and copying them back." << std::endl;
			auto cpu_faust = TransformHelper<FPP,Cpu>::eyeFaust(n, m);
			TransformHelper<FPP,GPU2>* gpu_faust = new TransformHelper<FPP,GPU2>(*cpu_faust, -1, nullptr /*TODO: dev_id and stream ?*/);
			delete cpu_faust;
			return gpu_faust;
		}

	template<typename FPP>
	void TransformHelper<FPP,GPU2>::get_fact(const faust_unsigned_int id,
			int* rowptr,
			int* col_ids,
			FPP* elts,
			faust_unsigned_int* nnz,
			faust_unsigned_int* num_rows,
			faust_unsigned_int* num_cols,
			const bool transpose /* = false*/) const
	{
		this->transform->get_fact(this->is_transposed?this->size()-id-1:id, rowptr, col_ids, elts, nnz, num_rows, num_cols, this->is_transposed ^ transpose);
		if(this->is_conjugate)
			Faust::conjugate(elts, *nnz);
	}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::get_fact(const faust_unsigned_int id,
				FPP* elts,
				faust_unsigned_int* num_rows,
				faust_unsigned_int* num_cols,
				const bool transpose /* = false*/) const
		{
			TransformHelperGen<FPP,GPU2>::get_fact(id, elts, num_rows, num_cols, transpose);
		}
}

#include "faust_TransformHelper_cat_gpu.hpp"
