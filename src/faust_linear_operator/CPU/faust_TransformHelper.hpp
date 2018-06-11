
namespace Faust {

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(const std::vector<MatGeneric<FPP,Cpu> *>& facts, const FPP lambda_, const bool optimizedCopy) : is_transposed(false)
	{
		Transform<FPP,Cpu>* transform = new Transform<FPP,Cpu>(facts, lambda_, optimizedCopy);
		this->transform = make_shared<Transform<FPP,Cpu>>(facts, lambda_, optimizedCopy);
	}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper() : is_transposed(false)
	{
		Transform<FPP,Cpu>* transform = new Transform<FPP,Cpu>();
		this->transform = make_shared<Transform<FPP,Cpu>>();
	}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(TransformHelper<FPP,Cpu>* th, bool transpose, bool conjugate)
		{
			this->transform = th->transform;
			this->is_transposed = transpose?!th->is_transposed:th->is_transposed;
			this->is_conjugate = conjugate?!th->is_conjugate:th->is_conjugate;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(TransformHelper<FPP,Cpu>* th)
		{
			this->transform = th->transform;
			this->is_transposed = th->is_transposed;
			this->is_conjugate = th->is_conjugate;
		}

	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const MatDense<FPP,Cpu> A) const
		{
			return this->transform->multiply(A, isTransposed2char());
		}

	template<typename FPP>
		Vect<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const Vect<FPP,Cpu> x) const
		{
			return this->transform->multiply(x, isTransposed2char());
		}

	template<typename FPP>
		Vect<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const Vect<FPP,Cpu> x, const bool transpose)
		{
			is_transposed ^= transpose;
			Vect<FPP,Cpu> v = this->transform->multiply(x, isTransposed2char());
			is_transposed ^= transpose;
			return v;
		}


	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const MatDense<FPP,Cpu> A, const bool transpose)
		{
			is_transposed ^= transpose;
			MatDense<FPP,Cpu> M = this->transform->multiply(A, isTransposed2char());
			is_transposed ^= transpose;
			return M;
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::push_back(const MatGeneric<FPP,Cpu>* M)
		{
			this->transform->push_back(M);
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,Cpu>::getNbRow() const
		{
			return is_transposed?transform->getNbCol():transform->getNbRow();
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,Cpu>::getNbCol() const
		{
			return is_transposed?transform->getNbRow():transform->getNbCol();
		}

	template<typename FPP>
		bool TransformHelper<FPP,Cpu>::isReal() const
		{
			return this->transform->isReal();
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,Cpu>::get_total_nnz() const
		{
			return this->transform->get_total_nnz();
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,Cpu>::size() const
		{
			return this->transform->size();
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::display() const
		{
			this->transform->Display();
		}

	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::get_fact(faust_unsigned_int id) const
		{
			MatGeneric<FPP,Cpu>* const factor_generic = this->transform->get_fact(is_transposed?size()-id-1:id);
			MatDense<FPP,Cpu> dense_factor;

			switch (factor_generic->getType())
			{
				case Dense :
					{
						MatDense<FPP,Cpu>* factor_dense_ptr = dynamic_cast<MatDense<FPP,Cpu>* > (factor_generic);
						dense_factor = (*factor_dense_ptr); //copy
					}
					break;

				case Sparse :
					{
						MatSparse<FPP,Cpu>* factor_sparse_ptr = dynamic_cast<MatSparse<FPP,Cpu>* > (factor_generic);
						dense_factor = (*factor_sparse_ptr); //copy
					}
					break;

				default:
					handleError("Faust::TransformHelper", "get_fact : unknown type of the factor matrix");
			}
			delete factor_generic;
			if(is_transposed) dense_factor.transpose();
			return dense_factor;
		}

	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::get_product() const {
			return this->transform->get_product(isTransposed2char());
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::save_mat_file(const char* filename) const
		{
			this->transform->save_mat_file(filename, is_transposed);
		}

	template<typename FPP>
		double TransformHelper<FPP,Cpu>::spectralNorm(const int nbr_iter_max, double threshold, int &flag) const
		{
			return this->transform->spectralNorm(nbr_iter_max, threshold, flag);
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::transpose()
		{
			return new TransformHelper<FPP,Cpu>(this, true, false);
		}

	template<typename FPP>
		bool TransformHelper<FPP,Cpu>::isTransposed() const
		{
			return is_transposed;
		}

	template<typename FPP>
		const char TransformHelper<FPP,Cpu>::isTransposed2char() const
		{
			return is_transposed?'T':'N';
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::~TransformHelper() {
			// transform is deleted auto. when no TransformHelper uses it (no more weak refs left)
		}

}
