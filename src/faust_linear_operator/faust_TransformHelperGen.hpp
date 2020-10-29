#include "faust_linear_algebra.h"
namespace Faust
{
	template<typename FPP, FDevice DEV>
		TransformHelperGen<FPP,DEV>::TransformHelperGen() : is_transposed(false), is_conjugate(false), is_sliced(false), is_fancy_indexed(false), transform(std::make_shared<Transform<FPP,DEV>>())
	{
	}

	template<typename FPP, FDevice DEV>
		const char TransformHelperGen<FPP,DEV>::isTransposed2char() const
		{
			return this->is_transposed?(this->is_conjugate?'H':'T'):'N';
		}

	template<typename FPP, FDevice DEV>
		bool TransformHelperGen<FPP,DEV>::isTransposed() const
		{
			return this->is_transposed;
		}

	template <typename FPP, FDevice DEV>
		void TransformHelperGen<FPP,DEV>::pack_factors(const faust_unsigned_int id, const PackDir dir)
		{
			if(dir == PACK_RIGHT)
				this->pack_factors(id, size()-1);
			else // dir == PACK_LEFT
				this->pack_factors(0, id);
		}

	template <typename FPP, FDevice DEV>
		void TransformHelperGen<FPP,DEV>::pack_factors()
		{
			//pack all factors in one
			this->pack_factors(0, this->size()-1);
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
}
