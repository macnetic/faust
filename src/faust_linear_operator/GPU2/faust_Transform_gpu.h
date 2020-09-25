#ifndef __FAUST_TRANSFORM_GPU2__
#define __FAUST_TRANSFORM_GPU2__
#ifdef USE_GPU_MOD
#include "faust_gpu_mod_utils.h"
#include "faust_constant.h"
#include "faust_MatGeneric_gpu.h"
#include "faust_MatDense_gpu.h"
#include <vector>

namespace Faust
{
	template<typename FPP, FDevice DEVICE> class Transform;

	template<typename FPP>
		class Transform<FPP,GPU2>
		{
			gm_MatArray_t gpu_mat_arr;
			public:
				Transform();
				Transform(const std::vector<MatGeneric<FPP,GPU2>*> &factors);
				Transform(const Transform<FPP,GPU2>& t);
				~Transform();
				void operator=(const Transform<FPP,GPU2>& t);
				void push_back(const MatGeneric<FPP,GPU2>*, bool copying=true);
				void push_first(const MatGeneric<FPP,GPU2>*, bool copying=true);
				void pop_front();
				void pop_back();
				void clear();
				MatGeneric<FPP,GPU2>* get_fact(int32_t id, bool cloning_fact=true) const;
				void get_facts(std::vector<MatGeneric<FPP,GPU2>*> &factors, bool cloning_facts=true) const;
				void transpose();
				int32_t getNbRow()const;
				int32_t getNbCol()const;
				void Display() const;
				int32_t size() const;
				faust_unsigned_int get_total_nnz() const;
				void update_total_nnz() const;
				void scalarMultiply(const FPP& alpha);
				MatDense<FPP,GPU2> get_product() const;
				void multiply(const Transform<FPP,GPU2> & A);
				void multiplyLeft(const Transform<FPP,GPU2> & A);
				Real<FPP> spectralNorm(int32_t nb_iter_max, float threshold, int& flag);
		};
}
#include "faust_Transform_gpu_double.hpp"
#endif
#endif
