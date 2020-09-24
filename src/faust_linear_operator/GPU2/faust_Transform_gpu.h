#ifndef __FAUST_TRANSFORM_GPU2__
#define __FAUST_TRANSFORM_GPU2__
#ifdef USE_GPU_MOD
#include "faust_gpu_mod_utils.h"
#include "faust_constant.h"
#include "faust_MatGeneric_gpu.h"
#include "faust_MatDense_gpu.h"
#include "faust_MatGeneric_gpu.h"
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
				~Transform();
				void push_back(const MatGeneric<FPP,GPU2>*, bool copying=true);
				int32_t getNbRow()const;
				int32_t getNbCol()const;
				void Display() const;
				int32_t size() const;
				MatDense<FPP,GPU2> get_product() const;
		};
}
#include "faust_Transform_gpu_double.hpp"
#endif
#endif
