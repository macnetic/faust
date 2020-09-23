#ifndef __FAUST_MATGENERIC_GPU__
#define __FAUST_MATGENERIC_GPU__
namespace Faust
{
	//TODO: this class is temporary, ideally MatSparse<FPP,GPU2> and MatDense<FPP,GPU2> should extend the MatGeneric<FPP, Device> class
	//TODO: keep this class until MatSparse<FPP,GPU2> and MatDense<FPP,GPU2> fully implement the MatGeneric<FPP, Device> methods
	// The interest of this class is mostly to make Transform capable of storing generic matrix
	template<typename FPP>
		class MatGeneric<FPP, GPU2>
		{
			protected:
				bool is_identity;
				bool is_zeros;
			public:
				virtual int32_t getNbRow() const=0;
				virtual int32_t getNbCol() const=0;
				MatGeneric();
				virtual ~MatGeneric();
		};

}
#include "faust_MatGeneric_gpu.hpp"
#endif
