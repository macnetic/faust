#include "faust_linear_algebra.h"
#include <iostream>

namespace Faust
{
	template<typename FPP, FDevice DEV>
		TransformHelperGen<FPP,DEV>::TransformHelperGen() : is_transposed(false), is_conjugate(false), is_sliced(false), is_fancy_indexed(false), transform(std::make_shared<Transform<FPP,DEV>>()), mul_order_opt_mode(0)

	{
		fancy_indices[0] = fancy_indices[1] = nullptr; // fancy_indices({nullptr, nullptr}) in ctor attribute init list is a GNU extension
	}

	template<typename FPP, FDevice DEV>
		TransformHelperGen<FPP,DEV>::TransformHelperGen(const TransformHelperGen<FPP,DEV>* th, bool transpose, bool conjugate) : TransformHelperGen<FPP,DEV>()
	{
		this->transform = th->transform;
		this->is_transposed = transpose?!th->is_transposed:th->is_transposed;
		this->is_conjugate = conjugate?!th->is_conjugate:th->is_conjugate;
		this->copy_slice_state(*th, transpose);
		this->copy_fancy_idx_state(*th, transpose);
    }


	template<typename FPP, FDevice DEV>
	void TransformHelperGen<FPP,DEV>::init_fancy_idx_transform(TransformHelper<FPP,DEV>* th, faust_unsigned_int* row_ids, faust_unsigned_int num_rows, faust_unsigned_int* col_ids, faust_unsigned_int num_cols)
	{
		this->transform = th->transform; //do not remove this line, necessary for eval*()
		this->copy_transconj_state(*th);
		//				handleError("Faust::TransformHelper::TransformHelper(TransformHelper,Slice)", "Fancy indexing overflows a Faust dimension.");
		this->fancy_num_cols = num_cols;
		this->fancy_num_rows = num_rows;
		this->fancy_indices[0] = new faust_unsigned_int[num_rows]; // deleted after evaluation or in destructor
		this->fancy_indices[1] = new faust_unsigned_int[num_cols];
		memcpy(this->fancy_indices[0], row_ids, num_rows*sizeof(faust_unsigned_int));
		memcpy(this->fancy_indices[1], col_ids, num_cols*sizeof(faust_unsigned_int));
		// the Faust is sliced, shift the indices to do as if the Faust wasn't sliced
		if(th->is_sliced)
		{
			if(th->slices[0].start_id != 0)
				for(int i=0; i < num_rows; i++)
				{
					this->fancy_indices[0][i] += th->slices[0].start_id;
				}
			if(th->slices[1].start_id != 0)
				for(int i=1; i < num_cols; i++)
				{
					this->fancy_indices[1][i] += th->slices[1].start_id;
				}
		}
		else if(th->is_fancy_indexed) // a Faust can't be sliced and fancy_indexed at the same time
		{
			// the Faust is already indexed, convert the indices taking account of th indices
			for(int i=0;i < num_rows; i++)
			{
				auto id = this->fancy_indices[0][i];
				if(id < th->fancy_num_rows)
					this->fancy_indices[0][i] = th->fancy_indices[0][id];
				// else error delayed to the sanity check below
			}
			for(int j=0;j < num_cols; j++)
			{
				auto id = this->fancy_indices[1][j];
				if(id < th->fancy_num_cols)
					this->fancy_indices[1][j] = th->fancy_indices[1][id];
				// else error delayed to the sanity check below
			}
		}
		// prevent overflows
		for(int i=0;i < num_rows; i++)
		{
			auto size = is_transposed?transform->getNbCol():transform->getNbRow();
			if(this->fancy_indices[0][i] > size)
				throw std::runtime_error("Faust indexing error: row index is greater than Faust number of rows.");

		}
		for(int j=0;j < num_cols; j++)
		{
			auto size = is_transposed?transform->getNbRow():transform->getNbCol();
			if(this->fancy_indices[1][j] > size)
				throw std::runtime_error("Faust indexing error: column index is greater than Faust number of columns.");
		}
		is_sliced = false; // absolute indices were computed, consider the Faust is not sliced in any case (even if it was before)
		this->is_fancy_indexed = true;
//		this->eval_fancy_idx_Transform(); // lazy indexing
//		copy_mul_mode_state(*th); // the best product method can be totally different after indexing
	}

	template<typename FPP, FDevice DEV>
		void TransformHelperGen<FPP,DEV>::init_sliced_transform(TransformHelper<FPP,DEV>* th, Slice s[2])

	{
		th->eval_fancy_idx_Transform();
		this->transform = th->transform; //do not remove this line, necessary for eval_sliced_Transform()
		this->copy_transconj_state(*th);
		if(! (s[0].belong_to(0, th->getNbRow()) || s[1].belong_to(0, th->getNbCol())))
			handleError("Faust::TransformHelper::TransformHelper(TransformHelper,Slice)", "Slice overflows a Faust dimension.");
		if(th->is_sliced)
		{
			// the original Faust is already sliced, combine slices
			this->slices[0] = Slice(s[0].start_id+th->slices[0].start_id, s[0].end_id+th->slices[0].start_id);
			this->slices[1] = Slice(s[1].start_id+th->slices[1].start_id, s[1].end_id+th->slices[1].start_id);
		}
		else
		{
			this->slices[0] = s[0];
			this->slices[1] = s[1];
		}
		this->is_sliced = true;
//		this->eval_sliced_Transform(); // lazy slicing
//		this->copy_mul_mode_state(*th); // the best product method can be totally different after indexing
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
			if(is_sliced && id == 0 && dim == 0)
				return slices[0].size();
			else if (is_sliced && id == size()-1 && dim == 1)
				return slices[1].size();
			else
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
		}

	template<typename FPP, FDevice DEV>
		faust_unsigned_int TransformHelperGen<FPP,DEV>::get_fact_nnz(const faust_unsigned_int id) const
		{
			if(id == 0 || id == size()-1)
			{
				if(is_sliced)
					const_cast<Faust::TransformHelperGen<FPP, DEV>*>(this)->eval_sliced_Transform();
				const_cast<Faust::TransformHelperGen<FPP, DEV>*>(this)->eval_fancy_idx_Transform();
			}
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
			if(id == 0 || id == size()-1)
			{
				const_cast<Faust::TransformHelperGen<FPP, DEV>*>(this)->eval_sliced_Transform();
				const_cast<Faust::TransformHelperGen<FPP, DEV>*>(this)->eval_fancy_idx_Transform();
			}
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
			if(id == 0 || id == size()-1)
			{
				const_cast<Faust::TransformHelperGen<FPP, DEV>*>(this)->eval_sliced_Transform();
				const_cast<Faust::TransformHelperGen<FPP, DEV>*>(this)->eval_fancy_idx_Transform();
			}
			this->transform->get_fact(this->is_transposed?this->size()-id-1:id, rowptr, col_ids, elts, nnz, num_rows, num_cols, this->is_transposed ^ transpose);
			if(this->is_conjugate)
				Faust::conjugate(elts, *nnz);
		}

	template<typename FPP, FDevice DEV>
		TransformHelper<FPP,DEV>* TransformHelperGen<FPP,DEV>::left(const faust_unsigned_int id, const bool copy /* default to false */) const
		{
			if(id >= this->size()) throw out_of_range("factor id is lower than zero or greater or equal to the size of Transform.");
			const_cast<Faust::TransformHelperGen<FPP, DEV>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelperGen<FPP, DEV>*>(this)->eval_fancy_idx_Transform();
			std::vector<Faust::MatGeneric<FPP,DEV>*> left_factors;
			for(int i=0; (faust_unsigned_int)i <= id; i++)
				left_factors.push_back(const_cast<Faust::MatGeneric<FPP,DEV>*>(this->get_gen_fact(i)));
			return new TransformHelper<FPP,DEV>(left_factors, FPP(1.0), false, copy, true);
		}

	template<typename FPP, FDevice DEV>
		TransformHelper<FPP,DEV>* TransformHelperGen<FPP,DEV>::right(const faust_unsigned_int id, const bool copy /* default to false */) const
		{
			if(id >= this->size()) throw out_of_range("factor id is lower than zero or greater or equal to the size of Transform.");
			const_cast<Faust::TransformHelperGen<FPP, DEV>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelperGen<FPP, DEV>*>(this)->eval_fancy_idx_Transform();
			std::vector<Faust::MatGeneric<FPP,DEV>*> right_factors;
			for(int i=id; (faust_unsigned_int)i < size(); i++)
				right_factors.push_back(const_cast<Faust::MatGeneric<FPP,DEV>*>(this->get_gen_fact(i)));
			return new TransformHelper<FPP,DEV>(right_factors, FPP(1.0), false, copy, true);
		}

	template<typename FPP, FDevice DEV>
		faust_unsigned_int TransformHelperGen<FPP,DEV>::getNbRow() const
		{
			if(this->is_sliced)
				return	this->slices[0].size();
			else if(this->is_fancy_indexed)
				return this->fancy_num_rows;
			else
				return this->is_transposed?this->transform->getNbCol():this->transform->getNbRow();
		}

	template<typename FPP, FDevice DEV>
		faust_unsigned_int TransformHelperGen<FPP,DEV>::getNbCol() const
		{
			if(this->is_sliced)
				return this->slices[1].size();
			else if(this->is_fancy_indexed)
				return this->fancy_num_cols;
			else
				return this->is_transposed?this->transform->getNbRow():this->transform->getNbCol();
		}

	template<typename FPP, FDevice DEV>
		void TransformHelperGen<FPP, DEV>::copy_slices(const TransformHelperGen<FPP, DEV> *th, const bool transpose /*=false*/)
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
		void TransformHelperGen<FPP,DEV>::copy_slice_state(const TransformHelperGen<FPP,DEV>& th, bool transpose /*=false*/)
		{
			this->is_sliced = th.is_sliced;
			if(th.is_sliced)
				copy_slices(&th, transpose);
		}

	template<typename FPP, FDevice DEV>
		void TransformHelperGen<FPP,DEV>::copy_fancy_idx_state(const TransformHelperGen<FPP,DEV>& th, bool transpose /*=false*/)
		{
			this->is_fancy_indexed = th.is_fancy_indexed;
			if(th.is_fancy_indexed)
			{
				if(fancy_indices[0] != nullptr)
					delete[] fancy_indices[0];
				if(fancy_indices[1] != nullptr)
					delete[] fancy_indices[1];
				int id0, id1;
				if(transpose)
				{
					fancy_num_rows = th.fancy_num_cols;
					fancy_num_cols = th.fancy_num_rows;
					id0 = 1, id1 = 0;
				}
				else
				{
					fancy_num_rows = th.fancy_num_rows;
					fancy_num_cols = th.fancy_num_cols;
					id0 = 0, id1 = 1;
				}
				fancy_indices[0] = new faust_unsigned_int[fancy_num_rows];
				fancy_indices[1] = new faust_unsigned_int[fancy_num_cols];
				memcpy(fancy_indices[0], th.fancy_indices[id0], fancy_num_rows*sizeof(faust_unsigned_int));
				memcpy(fancy_indices[1], th.fancy_indices[id1], fancy_num_cols*sizeof(faust_unsigned_int));
			}
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
			copy_fancy_idx_state(th);
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
		void Faust::TransformHelperGen<FPP, DEV>::eval_sliced_Transform()
		{
			eval_fancy_idx_Transform();
			if(is_sliced)
			{
				// new transform object (to avoid modifying original Faust::Transform*) but the matrices are shared (not copied)
				vector<Faust::MatGeneric<FPP, DEV>*> factors;
				for(auto ite=transform->begin(); ite != transform->end(); ite++)
					factors.push_back(*ite);
				this->transform = make_shared<Transform<FPP,DEV>>(factors, 1.0, false, /* cloning_fact */ false);
				auto first_fact = factors[0];
				Slice first_fact_slice, last_fact_slice;
				if(this->is_transposed)
				{
					first_fact_slice = this->slices[1];
					last_fact_slice = this->slices[0];
				}
				else
				{
					first_fact_slice = this->slices[0];
					last_fact_slice = this->slices[1];
				}
				// Replace only the first and last factors if necessary
				if(first_fact_slice.start_id != 0 || first_fact_slice.end_id != first_fact->getNbRow())
					// the row slice is smaller than the first factor num of rows, replace the first factor
					this->transform->replace(first_fact->get_rows(first_fact_slice.start_id, first_fact_slice.size()), 0);
				auto last_fact = *(transform->end()-1); // if first_fact is also the last fact, it may have changed since first_fact was read
				if(last_fact_slice.start_id != 0 || last_fact_slice.end_id != last_fact->getNbCol())
					// the column slice is smaller than the last factor num of cols, replace the first factor
					this->transform->replace(last_fact->get_cols(last_fact_slice.start_id, last_fact_slice.size()), size()-1);
				// now that the sliced was evaluated, consider this Transform as a basic one
				this->is_sliced = false;
			}
		}

	template<typename FPP, FDevice DEV>
		TransformHelper<FPP, DEV>* TransformHelperGen<FPP, DEV>::slice(faust_unsigned_int start_row_id, faust_unsigned_int end_row_id,
				faust_unsigned_int start_col_id, faust_unsigned_int end_col_id)
		{
			Slice sr(start_row_id, end_row_id);
			Slice sc(start_col_id, end_col_id);
			Slice s[2] = {sr, sc};
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
			if(is_fancy_indexed)
			{
				bool cloning_fact = false;
				faust_unsigned_int size = this->size();
				std::vector<MatGeneric<FPP,DEV>*> factors((size_t) size);
				MatGeneric<FPP,DEV>* gen_fac, *first_sub_fac, *last_sub_fac;
				gen_fac = transform->get_fact(0, cloning_fact);
				int id_left = 0;
				int id_right = 1;
				size_t left_ind_size = fancy_num_rows;
				size_t right_ind_size = fancy_num_cols;
				if(is_transposed)
				{
					id_left = 1;
					id_right = 0;
					left_ind_size = fancy_num_cols;
					right_ind_size = fancy_num_rows;
				}
				first_sub_fac = gen_fac->get_rows(fancy_indices[id_left], left_ind_size);
				if(cloning_fact)
					delete gen_fac;
				if(size > 1)
				{
					gen_fac = transform->get_fact(size-1, cloning_fact);
					last_sub_fac = gen_fac->get_cols(fancy_indices[id_right], right_ind_size);	//		std::cout << "---" << std::endl;
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
							gen_fac = transform->get_fact(i, cloning_fact);
							factors[i] = gen_fac;
						}
					}
					factors.insert(factors.begin()+(size-1), last_sub_fac);
					factors.resize(size);
				}
				else
				{ //only one factor
					last_sub_fac = first_sub_fac->get_cols(fancy_indices[id_right], right_ind_size);
					delete first_sub_fac;
					factors[0] = last_sub_fac;
					factors.resize(1);
				}
				transform = make_shared<Transform<FPP, DEV>>(factors, 1.0, false, cloning_fact);
				if(cloning_fact) {
					for(faust_unsigned_int i = 0; i < size; i++)
						delete factors[i];
				}
				is_fancy_indexed = false;
				delete[] fancy_indices[id_left];
				delete[] fancy_indices[id_right];
				fancy_indices[id_left] = fancy_indices[id_right] = nullptr;
			}
		}

	template<typename FPP, FDevice DEV>
		TransformHelper<FPP,DEV>* TransformHelperGen<FPP,DEV>::optimize_storage(const bool time/*=false*/)
		{
			eval_sliced_Transform();
			eval_fancy_idx_Transform();
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
			//TODO: should use copy_transconj_state
			pth->is_transposed = this->is_transposed;
			pth->is_conjugate = this->is_conjugate;
//			this->mul_order_opt_mode = th.mul_order_opt_mode; // do not copy the mul method because it might be inefficient for the transpose case
			return pth;
		}

	template<typename FPP, FDevice DEV>
		TransformHelper<FPP,DEV>* TransformHelperGen<FPP,DEV>::clone()
		{
			eval_sliced_Transform();
			eval_fancy_idx_Transform();
			auto clone = new TransformHelper<FPP,DEV>(this->transform->getData(), /* lambda_ */(FPP)1.0, /*optimizedCopy*/false, /*cloning_fact*/ true, /*internal_call*/ true);
			auto th = dynamic_cast<TransformHelper<FPP,DEV>*>(this);
			clone->copy_transconj_state(*th);
			clone->copy_slice_state(*th);
			clone->copy_mul_mode_state(*th);
			clone->copy_fancy_idx_state(*th);
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
				// do not use get_gen_fact to avoid (slicing/indexing) evaluation
				auto fac = this->transform->data[this->is_transposed?size()-i-1:i];
				if(! (csr && dynamic_cast<const MatSparse<FPP,DEV>*>(fac) || bsr && dynamic_cast<const MatBSR<FPP,DEV>*>(fac)))
					return false;
			}
			return true;
		}

	template<typename FPP, FDevice DEV>
		bool TransformHelperGen<FPP, DEV>::is_all_dense() const
		{
			for(int i=0;i<this->size();i++)
			{
				// do not use get_gen_fact to avoid (slicing/indexing) evaluation
				if(! dynamic_cast<const MatDense<FPP,DEV>*>(this->transform->data[this->is_transposed?size()-i-1:i]))
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

	template<typename FPP, FDevice DEV>
		TransformHelperGen<FPP, DEV>::~TransformHelperGen()
		{
			if(fancy_indices[0] != nullptr)
				delete[] fancy_indices[0];
			if(fancy_indices[1] != nullptr)
				delete[] fancy_indices[1];
			fancy_indices[0] = nullptr;
			fancy_indices[1] = nullptr;
		}
}
