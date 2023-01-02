#include "faust_conj.h"
#include <complex>

namespace Faust
{
	template<>
		double conj<double>(const double& e)
		{
			return std::real(std::conj(e));
		}

	template<>
		float conj<float>(const float& e)
		{
			return std::real(std::conj(e));
		}

	template<>
		std::complex<double> conj<std::complex<double>>(const std::complex<double>& e)
		{
			return std::conj(e);
		}

	template<>
		std::complex<float> conj<std::complex<float>>(const std::complex<float>& e)
		{
			return std::conj(e);
		}
}
