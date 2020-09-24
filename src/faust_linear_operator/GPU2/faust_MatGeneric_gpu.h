#ifndef __FAUST_MATGENERIC_GPU__
#define __FAUST_MATGENERIC_GPU__
namespace Faust
{
	template<typename FPP, FDevice DEVICE> class Transform;
	template<typename FPP, FDevice DEVICE> class MatGeneric;
	template<typename FPP>
		class Transform<FPP,GPU2>;
	//TODO: this class is temporary, ideally MatSparse<FPP,GPU2> and MatDense<FPP,GPU2> should extend the MatGeneric<FPP, Device> class
	//TODO: keep this class until MatSparse<FPP,GPU2> and MatDense<FPP,GPU2> fully implement the MatGeneric<FPP, Device> methods
	// The interest of this class is mostly to make Transform capable of storing generic matrix
	template<typename FPP>
		class MatGeneric<FPP, GPU2>
		{
			friend Transform<FPP,GPU2>; // need to access to get_gpu_mat_ptr
			protected:
				bool is_identity;
				bool is_zeros;
			public:
				virtual int32_t getNbRow() const=0;
				virtual int32_t getNbCol() const=0;
				virtual MatGeneric<FPP,GPU2>* clone(const int32_t dev_id=-1, const void* stream=nullptr) const=0;
				virtual void* get_gpu_mat_ptr() const=0;
				MatGeneric();
				virtual ~MatGeneric();
		};

}
#include "faust_MatGeneric_gpu.hpp"
#endif
