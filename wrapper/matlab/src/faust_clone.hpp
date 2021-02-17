#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_clone_gpu2cpu(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV>>(prhs[1]);

	if (nlhs > 1)
		mexErrMsgTxt("clone_gpu2cpu: too many output arguments");
	else
	{
		Faust::TransformHelper<SCALAR,Cpu>* cpu_th = core_ptr->tocpu();
		plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,Cpu>>(cpu_th);
	}

}

template <typename SCALAR, FDevice DEV>
void faust_clone_cpu2gpu(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,Cpu>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,Cpu>>(prhs[1]);

	if (nlhs > 1)
		mexErrMsgTxt("clone_cpu2gpu: too many output arguments");
	else
	{
		TransformHelper<SCALAR,DEV>* gpu_th = new TransformHelper<SCALAR, DEV>(core_ptr);
		plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV>>(gpu_th);
	}

}
