namespace Faust
{
	template<typename FPP,FDevice DEVICE> class Transform;

	template<typename FPP>
		TransformHelper<FPP,GPU2>::TransformHelper() : TransformHelperGen<FPP,GPU2>()
	{
	}

	template<typename FPP>
		TransformHelper<FPP,GPU2>::TransformHelper(const std::vector<MatGeneric<FPP,GPU2> *>& facts, const FPP lambda_/*= (FPP)1.0*/, const bool optimizedCopy/*=false*/, const bool cloning_fact /*= true*/, const bool internal_call/*=false*/) : TransformHelper<FPP,GPU2>()
	{
		for(auto f: facts)
		{
			this->push_back(f, false, cloning_fact);
		}
		this->multiply(lambda_);
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
		void TransformHelper<FPP,GPU2>::push_back(const MatGeneric<FPP,GPU2>* M, const bool optimizedCopy/*=false*/, const bool copying/*=true*/)
		{
			//optimizedCopy is ignored because not handled yet by Transform<FPP,GPU2> // TODO ? (it's not used by wrappers anyway)
			this->transform->push_back(M, copying);
		}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::push_back(const MatGeneric<FPP,Cpu>* M, const bool optimizedCopy/*=false*/, const int32_t dev_id/*=-1*/, const void* stream/*=nullptr*/)
		{

			MatGeneric<FPP,GPU2>* gpu_M = nullptr;
			MatDense<FPP,Cpu>* cpu_dM = nullptr;
			MatSparse<FPP,Cpu>* cpu_sM = nullptr;
			if(nullptr != (cpu_dM = dynamic_cast<MatDense<FPP,Cpu>*>(M)))
			{
				auto gpu_dM = new MatDense<FPP,GPU2>(M->getNbRow(), M->getNbCol(), cpu_dM->getData());
				gpu_M = gpu_dM;
			}
			else if(nullptr != (cpu_sM = dynamic_cast<MatSparse<FPP,Cpu>*>(M)))
			{
				auto gpu_sM = new MatSparse<FPP,GPU2>(M.getNbRow(), M.getNbCol(), cpu_sM.getNonZeros(), cpu_sM.getValuePtr(), cpu_sM.getRowPtr(), cpu_sM.getColInd(), dev_id);
				gpu_M = gpu_sM;
			}
			this->transform->push_back(gpu_M, false);
		}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::push_first(const MatGeneric<FPP,GPU2>* M, const bool optimizedCopy/*=false*/, const bool copying/*=true*/)
		{
			return this->transform->push_first(M, copying);
		}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::display() const
		{
			this->transform->Display();
		}

	template<typename FPP>
		std::string TransformHelper<FPP,GPU2>::to_string() const
		{
			return this->transform->to_string(/*this->is_transposed*/);
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
		MatDense<FPP,GPU2> TransformHelper<FPP,GPU2>::get_product()
		{
			return this->transform->get_product();
		}
	template<typename FPP>
		void TransformHelper<FPP,GPU2>::get_product(MatDense<FPP,GPU2>& M)
		{
			return this->transform->get_product(M);
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
		bool TransformHelper<FPP,GPU2>::is_fact_sparse(int id) const
		{
			return this->transform->is_fact_sparse(id);
		}

	template<typename FPP>
		bool TransformHelper<FPP,GPU2>::is_fact_dense(int id) const
		{
			return this->transform->is_fact_dense(id);
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
		MatDense<FPP,Cpu> TransformHelper<FPP,GPU2>::multiply(const Faust::MatDense<FPP,Cpu> &A, const bool transpose /* deft to false */, const bool conjugate/*=false*/)
		{
			MatDense<FPP,GPU2> M = this->multiply(MatDense<FPP,GPU2>(A), transpose, conjugate);
			return M.tocpu();
		}

	template<typename FPP>
		Vect<FPP,Cpu> TransformHelper<FPP,GPU2>::multiply(const Faust::Vect<FPP,Cpu> &A, const bool transpose /* deft to false */, const bool conjugate/*=false*/)
		{
//			Vect<FSFG,GPU2>::Vect(const faust_unsigned_int size,
//					const FSFG* cpu_data,
//					const bool no_alloc,
//					const int32_t dev_id,
//					const void* stream): MatDense<FSFG,GPU2>(size, 1, cpu_data, no_alloc, dev_id, stream)

			Vect<FPP,GPU2> gpu_A(A.size(), A.getData());
			Vect<FPP,GPU2> v = this->multiply(gpu_A/*, transpose, conjugate*/); //TODO: handle transpose and conjugate
			return v.tocpu();
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::multiply(const FPP& a)
		{
			this->transform->multiply(a);
			return this;
		}

	template<typename FPP>
		Vect<FPP,GPU2> TransformHelper<FPP,GPU2>::multiply(const Faust::Vect<FPP,GPU2>& a)
		{
			throw std::runtime_error("multiply is yet to implement in Faust C++ core for GPU.");
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
		void Faust::TransformHelper<FPP,GPU2>::pack_factors()
		{
			TransformHelperGen<FPP,Cpu>::pack_factors();
		}

	template<typename FPP>
		void Faust::TransformHelper<FPP,GPU2>::pack_factors(faust_unsigned_int start_id, faust_unsigned_int end_id)
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
			//TODO: factorize with TranformHelper<FPP,Cpu> into parent class TransformHelperGen.
			this->transform = th.transform;
			this->is_transposed = th.is_transposed;
			this->is_conjugate = th.is_conjugate;
			this->is_sliced = th.is_sliced;
			if(th.is_sliced)
				copy_slices(&th);
			this->mul_order_opt_mode = th.mul_order_opt_mode;
			this->Fv_mul_mode = th.Fv_mul_mode;
		}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::save_mat_file(const char* filename) const
		{
			this->transform->save_mat_file(filename, this->is_transposed, this->is_conjugate);
		}

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::tocpu(TransformHelper<FPP,Cpu>& cpu_transf)
		{
			auto t = this->transform->tocpu();
			for(auto fac: t)
				cpu_transf.push_back(fac, false, false);
			cpu_transf.display();
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::vertcat(const TransformHelper<FPP,GPU2>* G)
		{
			//TODO
			throw std::runtime_error("vertcat is yet to implement in Faust C++ core for GPU.");
			return nullptr;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::horzcat(const TransformHelper<FPP,GPU2>* G)
		{
			//TODO
			throw std::runtime_error("horzcat is yet to implement in Faust C++ core for GPU.");
			return nullptr;
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,GPU2>::get_total_nnz() const
		{
			return this->transform.get_total_nnz();
		}


	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::normalize(const int meth /* 1 for 1-norm, 2 for 2-norm (2-norm), -1 for inf-norm */) const
		{
			throw std::runtime_error("normalize is yet to implement in Faust C++ core for GPU.");
			return nullptr;
		}

	template<typename FPP>
		unsigned int TransformHelper<FPP,GPU2>::get_fact_nb_rows(const faust_unsigned_int id) const
		{
			throw std::runtime_error("get_fact_nb_rows is yet to implement in Faust C++ core for GPU.");
			return 0;
		}

	template<typename FPP>
		unsigned int TransformHelper<FPP,GPU2>::get_fact_nb_cols(const faust_unsigned_int id) const
		{
			throw std::runtime_error("get_fact_nb_cols is yet to implement in Faust C++ core for GPU.");
			return 0;
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,GPU2>::get_fact_nnz(const faust_unsigned_int id) const
		{
			throw std::runtime_error("get_fact_nnz is yet to implement in Faust C++ core for GPU.");
			return 0;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::transpose()
		{
			throw std::runtime_error("transpose is yet to implement in Faust C++ core for GPU.");
			return nullptr;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::conjugate()
		{
			throw std::runtime_error("conjugate is yet to implement in Faust C++ core for GPU.");
			return nullptr;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::adjoint()
		{
			throw std::runtime_error("adjoint is yet to implement in Faust C++ core for GPU.");
			return nullptr;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::right(const faust_unsigned_int id, const bool copy/*=false*/) const
		{
			throw std::runtime_error("right is yet to implement in Faust C++ core for GPU.");
			return nullptr;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::left(const faust_unsigned_int id, const bool copy/*=false*/) const
		{
			throw std::runtime_error("left is yet to implement in Faust C++ core for GPU.");
			return nullptr;
		}

	template<typename FPP>
		TransformHelper<FPP, GPU2>* TransformHelper<FPP, GPU2>::slice(faust_unsigned_int start_row_id, faust_unsigned_int end_row_id,
				faust_unsigned_int start_col_id, faust_unsigned_int end_col_id)
		{
			throw std::runtime_error("slice is yet to implement in Faust C++ core for GPU.");
			return nullptr;
		}

	template<typename FPP>
		TransformHelper<FPP, GPU2>* TransformHelper<FPP,GPU2>::fancy_index(faust_unsigned_int* row_ids, faust_unsigned_int num_rows, faust_unsigned_int* col_ids, faust_unsigned_int num_cols)
		{
			throw std::runtime_error("fancy_index is yet to implement in Faust C++ core for GPU.");
			return nullptr;
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,GPU2>::getNBytes() const
		{
			throw std::runtime_error("getNBytes is yet to implement in Faust C++ core for GPU.");
			return 0;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::optimize_time(const bool transp/*=false*/, const bool inplace/*=false*/, const int nsamples/*=1*/)
		{
			throw std::runtime_error("optimize_time is yet to implement in Faust C++ core for GPU.");
			return nullptr;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::optimize(const bool transp/*=false*/)
		{
			throw std::runtime_error("optimize is yet to implement in Faust C++ core for GPU.");
			return nullptr;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::optimize_storage(const bool time/*=true*/)
		{
			throw std::runtime_error("optimize_storage is yet to implement in Faust C++ core for GPU.");
			return nullptr;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::hadamardFaust(unsigned int n, const bool norma/*=true*/)
		{
			throw std::runtime_error("hadamardFaust is yet to implement in Faust C++ core for GPU.");
			return nullptr;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::fourierFaust(unsigned int n, const bool norma/*=true*/)
		{
			throw std::runtime_error("fourierFaust is yet to implement in Faust C++ core for GPU.");
			return nullptr;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::eyeFaust(unsigned int n, unsigned int m)
		{
			throw std::runtime_error("eyeFaust is yet to implement in Faust C++ core for GPU.");
			return nullptr;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::pruneout(const int nnz_tres, const int npasses/*=-1*/, const bool only_forward/*=false*/)
		{
			throw std::runtime_error("pruneout is yet to implement in Faust C++ core for GPU.");
			return nullptr;
		}

	template<typename FPP>
			void TransformHelper<FPP,GPU2>::get_fact(const faust_unsigned_int id,
					int* rowptr,
					int* col_ids,
					FPP* elts,
					faust_unsigned_int* nnz,
					faust_unsigned_int* num_rows,
					faust_unsigned_int* num_cols,
					const bool transpose/*=false*/) const
			{
				throw std::runtime_error("get_fact is yet to implement in Faust C++ core for GPU.");
			}
}
