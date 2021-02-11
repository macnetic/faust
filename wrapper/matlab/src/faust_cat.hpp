#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_vertcat(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);

	//TODO: check number of args and results
	Faust::TransformHelper<SCALAR,DEV>* right_term = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[2]);
	if(right_term->getNbCol() != core_ptr->getNbCol()) mexErrMsgTxt("The dimensions of the two Fausts must agree.");
	Faust::TransformHelper<SCALAR,DEV>* th = core_ptr->vertcat(right_term);
	plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(th);
#ifdef FAUST_VERBOSE
	printf("new ptr for Faust concatened vertically TransformHelper obj=%p\n", th);
	printf("isConjugate()=%d\n",th->isConjugate());
	printf("isTransposed()=%d\n", th->isTransposed());
#endif
}

template <typename SCALAR, FDevice DEV>
void faust_horzcat(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	//TODO: check number of args and results
	Faust::TransformHelper<SCALAR,DEV>* right_term = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[2]);
	if(right_term->getNbRow() != core_ptr->getNbRow()) mexErrMsgTxt("The dimensions of the two Fausts must agree.");
	Faust::TransformHelper<SCALAR,DEV>* th = core_ptr->horzcat(right_term);
	plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(th);
#ifdef FAUST_VERBOSE
	printf("new ptr for Faust concatened horizontally TransformHelper obj=%p\n", th);
	printf("isConjugate()=%d\n",th->isConjugate());
	printf("isTransposed()=%d\n", th->isTransposed());
#endif
}

