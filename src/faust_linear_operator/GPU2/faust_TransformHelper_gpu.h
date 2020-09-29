#ifndef __TRANSFORM_HELPER_GPU2__
#define __TRANSFORM_HELPER_GPU2__
#include "faust_constant.h"
#include "faust_Transform_gpu.h"
//#include "faust_Transform.h"
#include "faust_TransformHelperGen.h"
#include <memory>
namespace Faust
{
	template<typename FPP, FDevice DEVICE> class TransformHelper;
	template<typename FPP, FDevice DEVICE> class TransformHelperGen;
	template<typename FPP>
		class TransformHelper<FPP,GPU2> : public TransformHelperGen<FPP,GPU2>
		{
			public:
				TransformHelper();
#ifndef IGNORE_TRANSFORM_HELPER_VARIADIC_TPL
			template<typename ...GList> TransformHelper(GList& ... t);
#endif
				void push_back(const MatGeneric<FPP,GPU2>* M, const bool optimizedCopy=false, const bool copying=true);
				template<typename Head, typename ... Tail>
					void push_back_(Head& h, Tail&... t);
				void push_back_();
				void Display() const;
				MatDense<FPP,GPU2> get_product();
				MatDense<FPP,GPU2> multiply(const Faust::MatDense<FPP,GPU2> &A, const bool transpose=false, const bool conjugate=false);
				Real<FPP> normFro() const;
				faust_unsigned_int size() const;
				void update_total_nnz() const;
				Real<FPP> spectralNorm(int32_t nb_iter_max, float threshold, int& flag);
				bool is_fact_sparse(int id) const;
				bool is_fact_dense(int id) const;
				MatGeneric<FPP,GPU2>* get_gen_fact_nonconst(const faust_unsigned_int id) const;
				void update(const MatGeneric<FPP, GPU2>& M, const faust_unsigned_int id);
		};
}
#include "faust_TransformHelper_gpu.hpp"
#endif
