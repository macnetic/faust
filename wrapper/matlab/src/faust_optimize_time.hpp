#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_optimize_time(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	if(nrhs != 5)
		mexErrMsgTxt("optimize_time mex error: invalid number of arguments.");
	if(nlhs != 1 && nlhs != 0)
		mexErrMsgTxt("optimize_time mex error: this function doesn't return more than one argument.");
	bool transp = (bool) mxGetScalar(prhs[2]);
	bool inplace = (bool) mxGetScalar(prhs[3]);
	int nsamples = (int) mxGetScalar(prhs[4]);
	Faust::TransformHelper<SCALAR,DEV>* th = core_ptr->optimize_time(transp, inplace, nsamples);
	if(inplace /*th == nullptr*/)
		return;
	plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(th);
}

template <typename SCALAR, FDevice DEV>
void faust_optimize_time_prod(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	if(nrhs < 6)
		mexErrMsgTxt("optimize_time_prod mex error: invalid number of arguments.");
	if(nlhs != 1 && nlhs != 0)
		mexErrMsgTxt("optimize_time_prod mex error: this function doesn't return more than one argument.");
	bool transp = (bool) mxGetScalar(prhs[2]);
	bool inplace = (bool) mxGetScalar(prhs[3]);
	int nsamples = (int) mxGetScalar(prhs[4]);
    const mxArray *mat = nullptr;
    const MatGeneric<SCALAR, Cpu>* matGen = nullptr;
    Faust::MatSparse<SCALAR,Cpu> sp_mat;
    Faust::MatDense<SCALAR,Cpu> ds_mat;
    mat = prhs[5];
	if(mxIsSparse(mat))
	{
		mxArray2FaustspMat<SCALAR>(mat, sp_mat);
        matGen = &sp_mat;
    }
    else
    {
		SCALAR* ptr_data = nullptr;
		mxArray2Ptr(mat, ptr_data);
		const size_t mat_nrows = mxGetM(mat);
		const size_t mat_ncols  = mxGetN(mat);
        ds_mat.resize(mat_nrows, mat_ncols);
        memcpy(ds_mat.getData(), ptr_data, mat_ncols*mat_ncols*sizeof(SCALAR));
        delete [] ptr_data;
        matGen = &ds_mat;
    }
	Faust::TransformHelper<SCALAR,DEV>* th = core_ptr->optimize_time_prod(matGen, transp, inplace, nsamples);
	if(inplace /*th == nullptr*/)
		return;
	plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(th);
}
