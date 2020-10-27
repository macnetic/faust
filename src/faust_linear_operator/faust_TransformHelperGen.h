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
		public:
		TransformHelperGen();
#ifndef IGNORE_TRANSFORM_HELPER_VARIADIC_TPL
		template<typename ...GList> TransformHelperGen(GList& ... t);
#endif

		virtual void push_back(const MatGeneric<FPP,DEV>* M, const bool optimizedCopy=false, const bool copying=true)=0;

		const char isTransposed2char() const;
		void enable_gpu_meth_for_mul(){}; //TODO: remove later (it is only a special case of TransformHelper Cpu)

		virtual faust_unsigned_int size() const=0;
		virtual void pack_factors(faust_unsigned_int start_id, faust_unsigned_int end_id)=0;

		void pack_factors(const faust_unsigned_int id, const PackDir dir);
		void pack_factors();

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
	};
}
#include "faust_TransformHelperGen.hpp"
#endif


