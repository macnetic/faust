#ifndef __FAUST_MATGENERIC_GPU__
#define __FAUST_MATGENERIC_GPU__
#include "faust_constant.h"
namespace Faust
{
	template<typename FPP, FDevice DEVICE> class Transform;
	template<typename FPP, FDevice DEVICE> class MatGeneric;
	template<typename FPP>
		class Transform<FPP,GPU2>;
	// The interest of this class is mostly to make Transform capable of storing generic matrix
	// TODO: this class should extends MatGeneric<FPP,Device>
	template<typename FPP>
		class MatGeneric<FPP, GPU2>
		{
			friend Transform<FPP,GPU2>; // needs access to get_gpu_mat_ptr
			virtual void set_gpu_mat_ptr(void*)=0;
			protected:
				bool is_identity;
				bool is_zeros;
			public:
				virtual MatType getType() const=0;
				virtual int32_t getNbRow() const=0;
				virtual int32_t getNbCol() const=0;
				virtual MatGeneric<FPP,GPU2>* clone(const int32_t dev_id=-1, const void* stream=nullptr) const=0;
				virtual void* get_gpu_mat_ptr() const=0;
				virtual faust_unsigned_int getNonZeros() const=0;
				virtual void transpose()=0;
				virtual void conjugate()=0;
				virtual void adjoint()=0;
				MatGeneric();
				virtual ~MatGeneric();
		};

}
#include "faust_MatGeneric_gpu.hpp"
#endif
