#include "FaustCoreCpp.h"
#include "faust_TransformHelper.h"
using namespace Faust;
#ifdef USE_GPU_MOD
#include "faust_TransformHelper_gpu.h"
#endif

#ifdef BUILD_FLOAT_PYX
template<>
FaustCoreCpp<float, Cpu>* FaustCoreCpp2<std::complex<double>, float, Cpu>::to_float(FaustCoreCpp<std::complex<double>,Cpu>* core)
{
	throw std::runtime_error("Faust in complex single precision are not supported (yet).");
}

template<>
FaustCoreCpp<float, Cpu>* FaustCoreCpp2<double, float, Cpu>::to_float(FaustCoreCpp<double, Cpu>* core)
{
	auto th = core->transform->template cast<float>();
	auto ret_core = new FaustCoreCpp<float, Cpu>(th);
	return ret_core;
}

template<>
FaustCoreCpp<float, Cpu>* FaustCoreCpp2<float, float, Cpu>::to_float(FaustCoreCpp<float, Cpu>* core)
{
	auto th = new TransformHelper<float, Cpu>(*core->transform);
	auto ret_core = new FaustCoreCpp<float,Cpu>(th);
	return ret_core;
}

template<>
FaustCoreCpp<double, Cpu>* FaustCoreCpp2<float, double, Cpu>::to_double(FaustCoreCpp<float, Cpu>* core)
{
	auto th = core->transform->template cast<double>();
	auto ret_core = new FaustCoreCpp<double,Cpu>(th);
	return ret_core;
}
#endif

template<>
FaustCoreCpp<double, Cpu>* FaustCoreCpp2<std::complex<double>, double, Cpu>::to_double(FaustCoreCpp<std::complex<double>,Cpu>* core)
{
	throw std::runtime_error("No complex conversion to double, use Faust.real");
}



template<>
FaustCoreCpp<double, Cpu>* FaustCoreCpp2<double, double, Cpu>::to_double(FaustCoreCpp<double, Cpu>* core)
{
	auto th = new TransformHelper<double, Cpu>(*core->transform);
	auto ret_core = new FaustCoreCpp<double,Cpu>(th);
	return ret_core;
}

#ifdef USE_GPU_MOD
#ifdef BUILD_FLOAT_PYX
template<>
FaustCoreCpp<float, GPU2>* FaustCoreCpp2<std::complex<double>, float, GPU2>::to_float(FaustCoreCpp<std::complex<double>,GPU2>* core)
{
	throw std::runtime_error("Faust in complex single precision are not supported (yet).");
}

template<>
FaustCoreCpp<float, GPU2>* FaustCoreCpp2<double, float, GPU2>::to_float(FaustCoreCpp<double, GPU2>* core)
{
	auto th = core->transform->template cast<float>();
	auto ret_core = new FaustCoreCpp<float,GPU2>(th);
	return ret_core;
}

template<>
FaustCoreCpp<float, GPU2>* FaustCoreCpp2<float, float, GPU2>::to_float(FaustCoreCpp<float, GPU2>* core)
{
	auto th = new TransformHelper<float, GPU2>(*core->transform);
	auto ret_core = new FaustCoreCpp<float,GPU2>(th);
	return ret_core;
}

template<>
FaustCoreCpp<double, GPU2>* FaustCoreCpp2<float, double, GPU2>::to_double(FaustCoreCpp<float, GPU2>* core)
{
	auto th = core->transform->template cast<double>();
	auto ret_core = new FaustCoreCpp<double,GPU2>(th);
	return ret_core;
}

#endif

template<>
FaustCoreCpp<double, GPU2>* FaustCoreCpp2<std::complex<double>, double, GPU2>::to_double(FaustCoreCpp<std::complex<double>,GPU2>* core)
{
	throw std::runtime_error("No complex conversion to double, use Faust.real");
}
template<>
FaustCoreCpp<double, GPU2>* FaustCoreCpp2<double, double, GPU2>::to_double(FaustCoreCpp<double, GPU2>* core)
{
	auto th = new TransformHelper<double, GPU2>(*core->transform);
	auto ret_core = new FaustCoreCpp<double,GPU2>(th);
	return ret_core;
}
#endif
