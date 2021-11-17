#include "faust_TransformHelper.h"

namespace Faust
{
	template<>
		TransformHelper<float, Cpu>* TransformHelper<std::complex<float>, Cpu>::real()
		{
			throw std::runtime_error("real conversion from complex double to float is not yet supported.");
		}

}
