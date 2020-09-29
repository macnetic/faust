namespace Faust
{
	template<typename FPP,FDevice DEVICE> class Transform;

	template<typename FPP>
	TransformHelper<FPP,GPU2>::TransformHelper() : TransformHelperGen<FPP,GPU2>()
	{
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
		void TransformHelper<FPP,GPU2>::Display() const
		{
			this->transform->Display();
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
		Real<FPP> TransformHelper<FPP,GPU2>::normFro() const
		{
			return this->transform->get_product().norm();
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
		MatDense<FPP,GPU2> TransformHelper<FPP,GPU2>::multiply(const Faust::MatDense<FPP,GPU2> &A, const bool transpose /* deft to false */, const bool conjugate)
		{
			this->is_transposed ^= transpose;
			this->is_conjugate ^= conjugate;
			MatDense<FPP,GPU2> M = this->transform->multiply(A, this->isTransposed2char());
			this->is_transposed ^= transpose;
			this->is_conjugate ^= conjugate;
			return M;
		}
}
