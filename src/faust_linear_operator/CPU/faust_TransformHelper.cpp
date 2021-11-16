#include "faust_TransformHelper.h"

namespace Faust
{
	//TODO: a generic cpp.in would refactor things


	template<>
		TransformHelper<double, Cpu>* TransformHelper<double, Cpu>::real()
		{
			return this->clone();
		}

	template<>
		TransformHelper<float, Cpu>* TransformHelper<float, Cpu>::real()
		{
			return this->clone();
		}

	template<>
		TransformHelper<float, Cpu>* TransformHelper<std::complex<float>, Cpu>::real()
		{
			throw std::runtime_error("real conversion from complex double to float is not yet supported.");
		}

}
