#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_subsref(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	int start_row_id = (int) mxGetScalar(prhs[2]);
	int end_row_id = (int) mxGetScalar(prhs[3]);
	int start_col_id = (int) mxGetScalar(prhs[4]);
	int end_col_id = (int) mxGetScalar(prhs[5]);
	Faust::TransformHelper<SCALAR, DEV>* th = core_ptr->slice(start_row_id-1, end_row_id, start_col_id-1, end_col_id);
	plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(th);
}

template <typename SCALAR, FDevice DEV>
void faust_subsref_byvec(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	unsigned long int *row_pr;
	row_pr = (unsigned long int*) mxGetData(prhs[2]);

	unsigned long int *col_pr;
	col_pr = (unsigned long int*) mxGetData(prhs[3]);

	Faust::TransformHelper<SCALAR, DEV>* th = core_ptr->fancy_index(row_pr, mxGetNumberOfElements(prhs[2]), col_pr, mxGetNumberOfElements(prhs[3]));
	plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(th);
}
