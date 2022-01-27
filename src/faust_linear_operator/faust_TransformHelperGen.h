#ifndef __FAUST_TRANSFORM_HELPER_DEVICE__
#define __FAUST_TRANSFORM_HELPER_DEVICE__
#include "faust_Slice.h"
#include "faust_constant.h"
#include "faust_prod_opt.h"
#include <memory>

namespace Faust
{
	template<typename FPP,FDevice DEVICE> class Transform;
	template<typename FPP,FDevice DEVICE> class TransformHelper;
	template<typename FPP,FDevice DEVICE> class Vect;
	template<typename FPP,FDevice DEVICE> class MatDense;
	template<typename FPP,FDevice DEVICE> class MatGeneric;

	enum RandFaustType {
		DENSE,
		SPARSE,
		MIXED
	};

	template<typename FPP, FDevice DEV>
	class TransformHelperGen
	{
        friend TransformHelper<FPP,GPU2>; // needs to access is_Transposed of TransformHelper<FPP,Cpu>::tocpu
		public:
			TransformHelperGen();
            TransformHelperGen(const TransformHelperGen<FPP,DEV>* th, bool transpose, bool conjugate);
#ifndef IGNORE_TRANSFORM_HELPER_VARIADIC_TPL
			template<typename ...GList> TransformHelperGen(GList& ... t);
#endif


            virtual faust_unsigned_int getNbCol() const;
            virtual faust_unsigned_int getNbRow() const;
			virtual void push_back(const MatGeneric<FPP,DEV>* M, const bool optimizedCopy=false, const bool copying=true, const bool transpose=false, const bool conjugate=false)=0;

			const char isTransposed2char() const;
			bool isTransposed() const;
			bool isConjugate() const;
			bool isReal() const;
			bool is_all_sparse() const;
			bool is_all_dense() const;

			virtual faust_unsigned_int size() const=0;
			virtual void get_fact(const faust_unsigned_int &id,
					FPP* elts,
					faust_unsigned_int* num_rows,
					faust_unsigned_int* num_cols,
					const bool transpose = false) const;
			virtual void get_fact(const faust_unsigned_int id,
					int* rowptr,
					int* col_ids,
					FPP* elts,
					faust_unsigned_int* nnz,
					faust_unsigned_int* num_rows,
					faust_unsigned_int* num_cols,
					const bool transpose = false) const;
			virtual unsigned int get_fact_nb_rows(const faust_unsigned_int id) const;
			virtual unsigned int get_fact_nb_cols(const faust_unsigned_int id) const;
			virtual unsigned int get_fact_dim_size(const faust_unsigned_int id, unsigned short dim) const;
			virtual const MatGeneric<FPP,DEV>* get_gen_fact(const faust_unsigned_int id) const=0;
			virtual faust_unsigned_int get_fact_nnz(const faust_unsigned_int id) const;
			virtual bool is_fact_sparse(const faust_unsigned_int id) const;
			virtual bool is_fact_dense(const faust_unsigned_int id) const;
			virtual bool is_fact_bsr(const faust_unsigned_int id) const;
			virtual MatType get_fact_type(const faust_unsigned_int id) const;
			virtual void pack_factors(faust_unsigned_int start_id, faust_unsigned_int end_id, const int mul_order_opt_mode=DEFAULT_L2R)=0;

			virtual void pack_factors(const int mul_order_opt_mode=DEFAULT_L2R);
			/**
			  \brief Returns the left hand side factors of this from index 0 to id included (as a new TransformHelper obj).

*/
			virtual TransformHelper<FPP,DEV>* left(const faust_unsigned_int id, const bool copy=false) const;
			/**
			  \brief Returns the right hand side factors of this from index id to the size()-1 (as a new TransformHelper obj).

*/
			virtual TransformHelper<FPP,DEV>* right(const faust_unsigned_int id, const bool copy=false) const;

			void copy_slices(const TransformHelper<FPP,DEV>* th, const bool transpose = false);
			void copy_slice_state(const TransformHelper<FPP,DEV>& th);
			void copy_transconj_state(const TransformHelper<FPP,DEV>& th);
			virtual void copy_mul_mode_state(const TransformHelper<FPP,DEV>& th);
			void copy_state(const TransformHelper<FPP,DEV>& th);
			int get_mul_order_opt_mode() const;
			void eval_sliced_Transform();
			void eval_fancy_idx_Transform();
			virtual TransformHelper<FPP, DEV>* slice(faust_unsigned_int start_row_id, faust_unsigned_int end_row_id,
					faust_unsigned_int start_col_id, faust_unsigned_int end_col_id);
			TransformHelper<FPP, DEV>* fancy_index(faust_unsigned_int* row_ids, faust_unsigned_int num_rows, faust_unsigned_int* col_ids, faust_unsigned_int num_cols);
			virtual TransformHelper<FPP,DEV>* optimize_storage(const bool time=true);
			virtual TransformHelper<FPP,DEV>* clone();
			virtual void convertToSparse()=0;
			virtual void convertToDense()=0;
			/**
			 * \brief Returns Faust-matrix element without computing the Faust full array.
			 */
			virtual FPP get_item(faust_unsigned_int i, faust_unsigned_int j)=0;
			void get_item(faust_unsigned_int i, faust_unsigned_int j, MatDense<FPP, DEV>& out_mat, faust_unsigned_int &out_id);

		protected:
			void init_fancy_idx_transform(TransformHelper<FPP,DEV>* th, faust_unsigned_int* row_ids, faust_unsigned_int num_rows, faust_unsigned_int* col_ids, faust_unsigned_int num_cols);
			void init_sliced_transform(TransformHelper<FPP,DEV>* th, Slice s[2]);
			bool is_transposed;
			bool is_conjugate;
			bool is_sliced;
			Slice slices[2];
			bool is_fancy_indexed;
			faust_unsigned_int * fancy_indices[2];
			faust_unsigned_int fancy_num_rows;
			faust_unsigned_int fancy_num_cols;
			std::shared_ptr<Transform<FPP,DEV>> transform;
			int mul_order_opt_mode;
	};
}
#include "faust_TransformHelperGen.hpp"
#endif


