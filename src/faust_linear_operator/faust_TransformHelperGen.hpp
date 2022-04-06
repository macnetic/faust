#include "faust_linear_algebra.h"
namespace Faust
{
	template<typename FPP, FDevice DEV>
		TransformHelperGen<FPP,DEV>::TransformHelperGen() : is_transposed(false), is_conjugate(false), is_sliced(false), is_fancy_indexed(false), transform(std::make_shared<Transform<FPP,DEV>>()), mul_order_opt_mode(0)

	{
	}

	template<typename FPP, FDevice DEV>
		TransformHelperGen<FPP,DEV>::TransformHelperGen(const TransformHelperGen<FPP,DEV>* th, bool transpose, bool conjugate) : TransformHelperGen<FPP,DEV>()
	{
		this->transform = th->transform;
		this->is_transposed = transpose?!th->is_transposed:th->is_transposed;
		this->is_conjugate = conjugate?!th->is_conjugate:th->is_conjugate;
    }


	template<typename FPP, FDevice DEV>
	void TransformHelperGen<FPP,DEV>::init_fancy_idx_transform(TransformHelper<FPP,DEV>* th, faust_unsigned_int* row_ids, faust_unsigned_int num_rows, faust_unsigned_int* col_ids, faust_unsigned_int num_cols)
	{
		this->transform = th->transform; //do not remove this line, necessary for eval*()
		this->copy_transconj_state(*th);
		this->is_sliced = false;
		//TODO: check indices
		//				handleError("Faust::TransformHelper::TransformHelper(TransformHelper,Slice)", "Fancy indexing overflows a Faust dimension.");
		unsigned int id0=0, id1=1;
		this->fancy_num_cols = num_cols;
		this->fancy_num_rows = num_rows;
		if(this->is_transposed)
		{
			id0 = 1;
			id1 = 0;
			this->fancy_num_cols = num_rows;
			this->fancy_num_rows = num_cols;
		}
		this->fancy_indices[id0] = new faust_unsigned_int[num_rows];
		this->fancy_indices[id1] = new faust_unsigned_int[num_cols];
		this->is_fancy_indexed = true;
		memcpy(this->fancy_indices[id0], row_ids, num_rows*sizeof(faust_unsigned_int));
		memcpy(this->fancy_indices[id1], col_ids, num_cols*sizeof(faust_unsigned_int));
		this->eval_fancy_idx_Transform();
		delete[] this->fancy_indices[0];
		delete[] this->fancy_indices[1];
		copy_mul_mode_state(*th);
	}

	template<typename FPP, FDevice DEV>
		void TransformHelperGen<FPP,DEV>::init_sliced_transform(TransformHelper<FPP,DEV>* th, Slice s[2])

	{
		this->transform = th->transform; //do not remove this line, necessary for eval_sliced_transform()
		this->copy_transconj_state(*th);
		if(! (s[0].belong_to(0, th->getNbRow()) || s[1].belong_to(0, th->getNbCol())))
			handleError("Faust::TransformHelper::TransformHelper(TransformHelper,Slice)", "Slice overflows a Faust dimension.");
		this->slices[0] = s[0];
		this->slices[1] = s[1];
		this->is_sliced = true;
		this->eval_sliced_Transform();
		this->copy_mul_mode_state(*th);
	}

	template<typename FPP, FDevice DEV>
		const char TransformHelperGen<FPP,DEV>::isTransposed2char() const
		{
            return this->transposed2char(this->is_transposed, this->is_conjugate);
		}

	template<typename FPP, FDevice DEV>
		const char TransformHelperGen<FPP,DEV>::transposed2char(bool is_transposed, bool is_conjugate) const
		{
			return is_transposed?(is_conjugate?'H':'T'):'N';
		}

	template<typename FPP, FDevice DEV>
		bool TransformHelperGen<FPP,DEV>::isTransposed() const
		{
			return this->is_transposed;
		}

	template <typename FPP, FDevice DEV>
		void TransformHelperGen<FPP,DEV>::pack_factors(const int mul_order_opt_mode/*=DEFAULT_L2R*/)
		{
			//pack all factors in one
			this->pack_factors(0, this->size()-1, mul_order_opt_mode);
		}

	template<typename FPP, FDevice DEV>
		unsigned int TransformHelperGen<FPP,DEV>::get_fact_nb_rows(const faust_unsigned_int id) const
		{
			return get_fact_dim_size(id,0);
		}

	template<typename FPP, FDevice DEV>
		unsigned int TransformHelperGen<FPP,DEV>::get_fact_nb_cols(const faust_unsigned_int id) const
		{
			return get_fact_dim_size(id,1);
		}

	template<typename FPP, FDevice DEV>
		unsigned int TransformHelperGen<FPP,DEV>::get_fact_dim_size(const faust_unsigned_int id, unsigned short dim) const
		{
			//dim == 0 to get num of cols, >0 otherwise
			faust_unsigned_int rid; //real id
			if(this->is_transposed) {
				rid = size()-id-1;
				dim = !dim;
			}
			else {
				rid = id;
				//dim = dim;
			}
			MatGeneric<FPP,DEV>* mat;
			if(rid == 0 || rid == this->size()-1)
				mat = this->transform->get_fact(rid, false);
			else
				mat = this->transform->get_fact(rid, false);
			if(dim)
				return mat->getNbCol();
			else
				return mat->getNbRow();
		}

	template<typename FPP, FDevice DEV>
		faust_unsigned_int TransformHelperGen<FPP,DEV>::get_fact_nnz(const faust_unsigned_int id) const
		{
			if(id == 0 || id == this->size()-1)
				return this->transform->get_fact_nnz(this->is_transposed?size()-id-1:id);
			else
				return this->transform->get_fact_nnz(this->is_transposed?size()-id-1:id);
		}

	template<typename FPP, FDevice DEV>
		bool TransformHelperGen<FPP,DEV>::is_fact_sparse(const faust_unsigned_int id) const
		{
			return this->transform->is_fact_sparse(this->is_transposed?size()-id-1:id);
		}

	template<typename FPP, FDevice DEV>
		bool TransformHelperGen<FPP,DEV>::is_fact_bsr(const faust_unsigned_int id) const
		{
			return this->transform->is_fact_bsr(this->is_transposed?size()-id-1:id);
		}

	template<typename FPP, FDevice DEV>
		MatType TransformHelperGen<FPP,DEV>::get_fact_type(const faust_unsigned_int id) const
		{
			return this->transform->get_fact(this->is_transposed?size()-id-1:id, false)->getType();
		}


	template<typename FPP, FDevice DEV>
		bool TransformHelperGen<FPP,DEV>::is_fact_dense(const faust_unsigned_int id) const
		{
			return this->transform->is_fact_dense(this->is_transposed?size()-id-1:id);
		}

	template<typename FPP, FDevice DEV>
		void TransformHelperGen<FPP,DEV>::get_fact(const faust_unsigned_int &id,
				FPP* elts,
				faust_unsigned_int* num_rows,
				faust_unsigned_int* num_cols,
				const bool transpose /* default to false */) const
		{
			this->transform->get_fact(this->is_transposed?this->size()-id-1:id, elts, num_rows, num_cols, this->is_transposed ^ transpose);
			if(this->is_conjugate)
				conjugate(elts,*num_cols*(*num_rows));
		}

	template<typename FPP, FDevice DEV>
		void TransformHelperGen<FPP,DEV>::get_fact(const faust_unsigned_int id,
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

	template<typename FPP, FDevice DEV>
		TransformHelper<FPP,DEV>* TransformHelperGen<FPP,DEV>::left(const faust_unsigned_int id, const bool copy /* default to false */) const
		{
			if(id >= this->size()) throw out_of_range("factor id is lower than zero or greater or equal to the size of Transform.");
			std::vector<Faust::MatGeneric<FPP,DEV>*> left_factors;
			for(int i=0; (faust_unsigned_int)i <= id; i++)
				left_factors.push_back(const_cast<Faust::MatGeneric<FPP,DEV>*>(this->get_gen_fact(i)));
			return new TransformHelper<FPP,DEV>(left_factors, FPP(1.0), false, copy, true);
		}

	template<typename FPP, FDevice DEV>
		TransformHelper<FPP,DEV>* TransformHelperGen<FPP,DEV>::right(const faust_unsigned_int id, const bool copy /* default to false */) const
		{
			if(id >= this->size()) throw out_of_range("factor id is lower than zero or greater or equal to the size of Transform.");
			std::vector<Faust::MatGeneric<FPP,DEV>*> right_factors;
			for(int i=id; (faust_unsigned_int)i < size(); i++)
				right_factors.push_back(const_cast<Faust::MatGeneric<FPP,DEV>*>(this->get_gen_fact(i)));
			return new TransformHelper<FPP,DEV>(right_factors, FPP(1.0), false, copy, true);
		}

	template<typename FPP, FDevice DEV>
		faust_unsigned_int TransformHelperGen<FPP,DEV>::getNbRow() const
		{
			if(this->is_sliced){
				faust_unsigned_int id = this->is_transposed?1:0;
				return	this->slices[id].end_id-this->slices[id].start_id;
			}
			else
				return this->is_transposed?this->transform->getNbCol():this->transform->getNbRow();
		}

	template<typename FPP, FDevice DEV>
		faust_unsigned_int TransformHelperGen<FPP,DEV>::getNbCol() const
		{
			if(this->is_sliced)
			{
				faust_unsigned_int id = this->is_transposed?0:1;
				return this->slices[id].end_id-this->slices[id].start_id;
			}
			else
				return this->is_transposed?this->transform->getNbRow():this->transform->getNbCol();
		}

	template<typename FPP, FDevice DEV>
		void TransformHelperGen<FPP, DEV>::copy_slices(const TransformHelper<FPP, DEV> *th, const bool transpose /*=false*/)
		{
			if(transpose)
			{
				this->slices[1].copy(th->slices[0]);
				this->slices[0].copy(th->slices[1]);
			}
			else {
				this->slices[0].copy(th->slices[0]);
				this->slices[1].copy(th->slices[1]);
			}
		}

	template<typename FPP, FDevice DEV>
		void TransformHelperGen<FPP,DEV>::copy_slice_state(const TransformHelper<FPP,DEV>& th)
		{
			this->is_sliced = th.is_sliced;
			if(th.is_sliced)
				copy_slices(&th);
		}

	template<typename FPP, FDevice DEV>
		void TransformHelperGen<FPP,DEV>::copy_transconj_state(const TransformHelper<FPP,DEV>& th)
		{
			this->is_transposed = th.is_transposed;
			this->is_conjugate = th.is_conjugate;
		}

	template<typename FPP, FDevice DEV>
		void TransformHelperGen<FPP,DEV>::copy_state(const TransformHelper<FPP,DEV>& th)
		{
			transform = th.transform;
			copy_transconj_state(th);
			copy_slice_state(th);
			copy_mul_mode_state(th);
		}

	template<typename FPP, FDevice DEV>
		void TransformHelperGen<FPP,DEV>::copy_mul_mode_state(const TransformHelper<FPP,DEV>& th)
		{
			this->mul_order_opt_mode = th.mul_order_opt_mode;
		}

	template<typename FPP, FDevice DEV>
		int TransformHelperGen<FPP,DEV>::get_mul_order_opt_mode() const
		{
			return this->mul_order_opt_mode;
		}

	template<typename FPP, FDevice DEV>
		void TransformHelperGen<FPP, DEV>::eval_sliced_Transform()
		{
			bool cloning_fact = true;
			std::vector<MatGeneric<FPP,DEV>*> factors((size_t) this->size());
			faust_unsigned_int size = this->size();
			MatGeneric<FPP,DEV>* gen_fac, *first_sub_fac, *last_sub_fac;
			gen_fac = this->transform->get_fact(0, cloning_fact);
			first_sub_fac = gen_fac->get_rows(this->slices[0].start_id, this->slices[0].end_id-this->slices[0].start_id);
			//		first_sub_fac->Display();
			//
			if(cloning_fact)
				delete gen_fac;
			if(size > 1)
			{
				gen_fac = this->transform->get_fact(size-1, cloning_fact);
				last_sub_fac = gen_fac->get_cols(this->slices[1].start_id, this->slices[1].end_id-this->slices[1].start_id);
				//		std::cout << "---" << std::endl;
				//		last_sub_fac->Display();
				if(cloning_fact)
					delete gen_fac;
				factors.reserve(size);
				factors.insert(factors.begin(), first_sub_fac);
				if(size > 2)
				{
					for(faust_unsigned_int i = 1; i < size-1; i++)
					{
						gen_fac = this->transform->get_fact(i, cloning_fact);
						factors[i] = gen_fac;
					}

				}
				factors.insert(factors.begin()+(size-1), last_sub_fac);
				factors.resize(size);
			}
			else
			{ //only one factor
				last_sub_fac = first_sub_fac->get_cols(this->slices[1].start_id, this->slices[1].end_id-this->slices[1].start_id);
				delete first_sub_fac;
				factors[0] = last_sub_fac;
				factors.resize(1);
			}
			this->transform = make_shared<Transform<FPP,DEV>>(factors, 1.0, false, cloning_fact);
			if(cloning_fact)
			{
				for(faust_unsigned_int i = 0; i < size; i++)
					delete factors[i];
			}
		}

	template<typename FPP, FDevice DEV>
		TransformHelper<FPP, DEV>* TransformHelperGen<FPP, DEV>::slice(faust_unsigned_int start_row_id, faust_unsigned_int end_row_id,
				faust_unsigned_int start_col_id, faust_unsigned_int end_col_id)
		{
			Slice sr(start_row_id, end_row_id);
			Slice sc(start_col_id, end_col_id);
			Slice s[2];
			if(this->is_transposed)
			{
				s[0] = sc;
				s[1] = sr;
			}
			else
			{
				s[0] = sr;
				s[1] = sc;
			}
			return new TransformHelper<FPP, DEV>(dynamic_cast<TransformHelper<FPP, DEV>*>(this), s);
		}

	template<typename FPP, FDevice DEV>
		TransformHelper<FPP, DEV>* TransformHelperGen<FPP,DEV>::fancy_index(faust_unsigned_int* row_ids, faust_unsigned_int num_rows, faust_unsigned_int* col_ids, faust_unsigned_int num_cols)
		{
			return new TransformHelper<FPP,DEV>(dynamic_cast<TransformHelper<FPP, DEV>*>(this), row_ids, num_rows, col_ids, num_cols);
		}

	template<typename FPP, FDevice DEV>
		void TransformHelperGen<FPP, DEV>::eval_fancy_idx_Transform()
		{
			bool cloning_fact = false;
			faust_unsigned_int size = this->size();
			std::vector<MatGeneric<FPP,DEV>*> factors((size_t) size);
			MatGeneric<FPP,DEV>* gen_fac, *first_sub_fac, *last_sub_fac;
			gen_fac = this->transform->get_fact(0, cloning_fact);
			//				first_sub_fac = gen_fac->get_rows(this->slices[0].start_id, this->slices[0].end_id-this->slices[0].start_id);
			//		first_sub_fac->Display();
			first_sub_fac = gen_fac->get_rows(this->fancy_indices[0], this->fancy_num_rows);
			if(cloning_fact)
				delete gen_fac;
			if(size > 1) {
				gen_fac = this->transform->get_fact(size-1, cloning_fact);
				//					last_sub_fac = gen_fac->get_cols(this->slices[1].start_id, this->slices[1].end_id-this->slices[1].start_id);
				last_sub_fac = gen_fac->get_cols(this->fancy_indices[1], this->fancy_num_cols);	//		std::cout << "---" << std::endl;
				//		last_sub_fac->Display();
				if(cloning_fact)
					delete gen_fac;
				factors.reserve(size);
				factors.insert(factors.begin(), first_sub_fac);
				if(size > 2)
				{
					auto it = factors.begin();
					for(faust_unsigned_int i = 1; i < size-1; i++)
					{
						gen_fac = this->transform->get_fact(i, cloning_fact);
						factors[i] = gen_fac;
					}
				}
				factors.insert(factors.begin()+(size-1), last_sub_fac);
				factors.resize(size);
			}
			else { //only one factor
				last_sub_fac = first_sub_fac->get_cols(this->fancy_indices[1], this->fancy_num_cols);
				delete first_sub_fac;
				factors[0] = last_sub_fac;
				factors.resize(1);
			}
			this->transform = make_shared<Transform<FPP, DEV>>(factors, 1.0, false, cloning_fact);
			if(cloning_fact) {
				for(faust_unsigned_int i = 0; i < size; i++)
					delete factors[i];
			}
		}

	template<typename FPP, FDevice DEV>
		TransformHelper<FPP,DEV>* TransformHelperGen<FPP,DEV>::optimize_storage(const bool time/*=false*/)
		{
//			throw std::runtime_error("optimize_storage is yet to implement in Faust C++ core for GPU.");
//			return nullptr;
			Faust::MatDense<FPP,DEV> *dfac = nullptr;
			Faust::MatSparse<FPP,DEV> *sfac = nullptr;
			faust_unsigned_int sparse_weight, dense_weight;
			std::vector<Faust::MatGeneric<FPP,DEV>*> opt_factors;

			for(auto fac : *(this->transform))
			{
				if((dfac = dynamic_cast<Faust::MatDense<FPP,DEV>*>(fac)))
					sfac = nullptr;
				else
					sfac = dynamic_cast<Faust::MatSparse<FPP,DEV>*>(fac);

				if(time)
				{
					if((dfac = dynamic_cast<Faust::MatDense<FPP,DEV>*>(fac)))
					{
						opt_factors.push_back(dfac->Clone(true));
					}
					else
					{
						opt_factors.push_back(sfac->Clone(true));
					}
				}
				else
				{ // storage size is the main criterion
					sparse_weight = fac->getNonZeros()*(sizeof(FPP)+sizeof(int))+(fac->getNbRow()+1)*sizeof(int);
					dense_weight = fac->getNbCol()*fac->getNbRow()*sizeof(FPP);
					if(sparse_weight < dense_weight)
					{
						// choose CSR format
						if(dfac)
							opt_factors.push_back(new Faust::MatSparse<FPP, DEV>(*dfac));
						else
							opt_factors.push_back(fac);
					}
					else
					{
						// choose dense format
						if(sfac)
							opt_factors.push_back(new Faust::MatDense<FPP, DEV>(*sfac));
						else
							opt_factors.push_back(fac);
					}
				}
			}
			TransformHelper<FPP,DEV> *pth = new TransformHelper<FPP,DEV>(opt_factors, FPP(1), false, false, true);
			//TODO: should use copy_transconj_state copy_slice_state
			pth->is_transposed = this->is_transposed;
			pth->is_conjugate = this->is_conjugate;
//			this->mul_order_opt_mode = th.mul_order_opt_mode; // do not copy the mul method because it might be inefficient for the transpose case
			return pth;
		}

	template<typename FPP, FDevice DEV>
		TransformHelper<FPP,DEV>* TransformHelperGen<FPP,DEV>::clone()
		{
			auto clone = new TransformHelper<FPP,DEV>(this->transform->getData(), /* lambda_ */(FPP)1.0, /*optimizedCopy*/false, /*cloning_fact*/ true, /*internal_call*/ true);
			auto th = dynamic_cast<TransformHelper<FPP,DEV>*>(this);
			clone->copy_transconj_state(*th);
			clone->copy_slice_state(*th);
			clone->copy_mul_mode_state(*th);
			return clone;
		}

	template<typename FPP, FDevice DEV>
		bool TransformHelperGen<FPP,DEV>::isConjugate() const
		{
			return this->is_conjugate;
		}

	template<typename FPP, FDevice DEV>
		bool TransformHelperGen<FPP,DEV>::isReal() const
		{
			return typeid(FPP) == typeid(double) || typeid(FPP) == typeid(float);
		}

	template<typename FPP, FDevice DEV>
		bool TransformHelperGen<FPP, DEV>::is_all_sparse(bool csr/*=true*/, bool bsr/*=false*/) const
		{
			if(this->size() == 0)
				return false;
			for(int i=0;i<this->size();i++)
			{
				if(! (csr && dynamic_cast<const MatSparse<FPP,DEV>*>(this->get_gen_fact(i)) || bsr && dynamic_cast<const MatBSR<FPP,DEV>*>(this->get_gen_fact(i))))
					return false;
			}
			return true;
		}

	template<typename FPP, FDevice DEV>
		bool TransformHelperGen<FPP, DEV>::is_all_dense() const
		{
			for(int i=0;i<this->size();i++)
			{
				if(! dynamic_cast<const MatDense<FPP,DEV>*>(this->get_gen_fact(i)))
					return false;
			}
			return true;
		}

	template<typename FPP, FDevice DEV>
		void TransformHelperGen<FPP, DEV>::get_item(faust_unsigned_int i, faust_unsigned_int j, MatDense<FPP, DEV> &out, faust_unsigned_int &out_id)
		{
			FPP item;
			TransformHelper<FPP, DEV> *th;
			if(this->getNbCol() < this->getNbRow())
			{
				th = this->slice(i, i+1, 0, this->getNbCol());
				th->transpose();
				out_id = j;
			}
			else
			{
				th = this->slice(0, this->getNbRow(), j, j+1);
				out_id = i;
			}
			out = th->get_product();
			delete th;
		}
}
