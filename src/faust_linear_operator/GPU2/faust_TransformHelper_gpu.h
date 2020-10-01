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
				TransformHelper(const std::vector<MatGeneric<FPP,GPU2> *>& facts, const FPP lambda_ = (FPP)1.0, const bool optimizedCopy=false, const bool cloning_fact = true, const bool internal_call=false);
#ifndef IGNORE_TRANSFORM_HELPER_VARIADIC_TPL
				template<typename ...GList> TransformHelper(GList& ... t);
#endif
				void push_back(const MatGeneric<FPP,GPU2>* M, const bool optimizedCopy=false, const bool copying=true);
				void pop_front();
				void pop_back();
				void push_first(const MatGeneric<FPP,GPU2>*, const bool optimizedCopy=false, const bool copying=true);
				template<typename Head, typename ... Tail>
					void push_back_(Head& h, Tail&... t);
				void push_back_();
				void display() const;
				MatDense<FPP,GPU2> get_product();
				void get_product(MatDense<FPP,GPU2>& M);
				MatDense<FPP,GPU2> multiply(const Faust::MatDense<FPP,GPU2> &A, const bool transpose=false, const bool conjugate=false);
				TransformHelper<FPP,GPU2>* multiply(const FPP& a);
				Real<FPP> normFro() const;
				faust_unsigned_int size() const;
				void update_total_nnz() const;
				Real<FPP> spectralNorm(int32_t nb_iter_max, float threshold, int& flag);
				bool is_fact_sparse(int id) const;
				bool is_fact_dense(int id) const;
				MatGeneric<FPP,GPU2>* get_gen_fact_nonconst(const faust_unsigned_int id) const;
				void pack_factors(faust_unsigned_int start_id, faust_unsigned_int end_id) {} //TODO
				void pack_factors() {} //TODO
				void pack_factors(const faust_unsigned_int id, const PackDir dir) {}//TODO
				void update(const MatGeneric<FPP, GPU2>& M, const faust_unsigned_int id);
				typename Transform<FPP,GPU2>::iterator begin() const;
				typename Transform<FPP,GPU2>::iterator end() const;
		};
}
#include "faust_TransformHelper_gpu.hpp"
#endif
