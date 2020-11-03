#ifndef __FAUST_TRANSFORM_HELPER_DEVICE__
#define __FAUST_TRANSFORM_HELPER_DEVICE__
#include "faust_Slice.h"
#include "faust_constant.h"
#include <memory>

namespace Faust
{
	template<typename FPP,FDevice DEVICE> class Transform;
	template<typename FPP,FDevice DEVICE> class TransformHelper;
	template<typename FPP,FDevice DEVICE> class Vect;
	template<typename FPP,FDevice DEVICE> class MatDense;
	template<typename FPP,FDevice DEVICE> class MatGeneric;

	enum PackDir {
		PACK_LEFT,
		PACK_RIGHT
	};

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


            faust_unsigned_int getNbCol() const;
            faust_unsigned_int getNbRow() const;
			virtual void push_back(const MatGeneric<FPP,DEV>* M, const bool optimizedCopy=false, const bool copying=true, const bool transpose=false, const bool conjugate=false)=0;

			const char isTransposed2char() const;
			bool isTransposed() const;
			void enable_gpu_meth_for_mul(){}; //TODO: remove later (it is only a special case of TransformHelper Cpu)

			virtual faust_unsigned_int size() const=0;
			void get_fact(const faust_unsigned_int &id,
					FPP* elts,
					faust_unsigned_int* num_rows,
					faust_unsigned_int* num_cols,
					const bool transpose = false) const;
			void get_fact(const faust_unsigned_int id,
					int* rowptr,
					int* col_ids,
					FPP* elts,
					faust_unsigned_int* nnz,
					faust_unsigned_int* num_rows,
					faust_unsigned_int* num_cols,
					const bool transpose = false) const;
			unsigned int get_fact_nb_rows(const faust_unsigned_int id) const;
			unsigned int get_fact_nb_cols(const faust_unsigned_int id) const;
			unsigned int get_fact_dim_size(const faust_unsigned_int id, unsigned short dim) const;
			virtual const MatGeneric<FPP,DEV>* get_gen_fact(const faust_unsigned_int id) const=0;
			faust_unsigned_int get_fact_nnz(const faust_unsigned_int id) const;
			bool is_fact_sparse(const faust_unsigned_int id) const;
			bool is_fact_dense(const faust_unsigned_int id) const;

			virtual void pack_factors(faust_unsigned_int start_id, faust_unsigned_int end_id)=0;

			void pack_factors(const faust_unsigned_int id, const PackDir dir);
			void pack_factors();
			/**
			  \brief Returns the left hand side factors of this from index 0 to id included (as a new TransformHelper obj).

*/
			TransformHelper<FPP,DEV>* left(const faust_unsigned_int id, const bool copy=false) const;
			/**
			  \brief Returns the right hand side factors of this from index id to the size()-1 (as a new TransformHelper obj).

*/
			TransformHelper<FPP,DEV>* right(const faust_unsigned_int id, const bool copy=false) const;

			void copy_slices(const TransformHelper<FPP,DEV>* th, const bool transpose = false);
			void copy_slice_state(const TransformHelper<FPP,DEV>& th);
			void copy_transconj_state(const TransformHelper<FPP,DEV>& th);
			virtual void copy_mul_mode_state(const TransformHelper<FPP,DEV>& th);
			void copy_state(const TransformHelper<FPP,DEV>& th);
			int get_mul_order_opt_mode() const;
			int get_Fv_mul_mode() const;

		protected:
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
			int Fv_mul_mode;

	};
}
#include "faust_TransformHelperGen.hpp"
#endif


