#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_rand(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{

	if(nlhs!=1)
		mexErrMsgTxt("rand(): 1 variable result is expected.");
	if(nrhs < 6)
		mexErrMsgTxt("rand(): wrong number of arguments (must be 6 to 7)");

	int num_rows = mxGetScalar(prhs[1]), num_cols =mxGetScalar(prhs[2]);
	faust_unsigned_int t = (faust_unsigned_int) mxGetScalar(prhs[3]),
					   min_num_factors = (faust_unsigned_int) mxGetScalar(prhs[4]),
					   max_num_factors = (faust_unsigned_int) mxGetScalar(prhs[5]),
					   min_dim_size = (faust_unsigned_int) mxGetScalar(prhs[6]),
					   max_dim_size = (faust_unsigned_int) mxGetScalar(prhs[7]);
	float density = (float) mxGetScalar(prhs[8]);
	bool per_row = (bool) mxGetScalar(prhs[9]);

	Faust::TransformHelper<SCALAR,DEV>* F = Faust::TransformHelper<SCALAR,DEV>::randFaust(
			num_rows,
			num_cols,
			RandFaustType(t),
			min_num_factors,
			max_num_factors,
			min_dim_size,
			max_dim_size,
			density,
			per_row);

	if(F) //not NULL
		plhs[0]=convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(F);
	else {
		plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
		double* ptr_out = (double*) mxGetData(plhs[0]);
		ptr_out[0]=(double) 0;
	}
}
