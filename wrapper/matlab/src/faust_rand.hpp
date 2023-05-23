#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_rand(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{

	if(nlhs!=1)
		mexErrMsgTxt("rand(): 1 variable result is expected.");
	if(nrhs < 10)
		mexErrMsgTxt("rand(): wrong number of arguments (must be at least 10)");

	int num_rows = mxGetScalar(prhs[1]), num_cols =mxGetScalar(prhs[2]);
	faust_unsigned_int t = (faust_unsigned_int) mxGetScalar(prhs[3]),
					   min_num_factors = (faust_unsigned_int) mxGetScalar(prhs[4]),
					   max_num_factors = (faust_unsigned_int) mxGetScalar(prhs[5]),
					   min_dim_size = (faust_unsigned_int) mxGetScalar(prhs[6]),
					   max_dim_size = (faust_unsigned_int) mxGetScalar(prhs[7]);
	float density = (float) mxGetScalar(prhs[8]);
	bool per_row = (bool) mxGetScalar(prhs[9]);
	unsigned int seed = 0;

	if(nrhs >= 11)
		seed = (unsigned int) mxGetScalar(prhs[10]);

	Faust::TransformHelper<SCALAR,DEV>* F = Faust::TransformHelper<SCALAR,DEV>::randFaust(
			num_rows,
			num_cols,
			RandFaustType(t),
			min_num_factors,
			max_num_factors,
			min_dim_size,
			max_dim_size,
			density,
			per_row,
			seed);

	if(F) //not NULL
		plhs[0]=convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(F);
	else {
		plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
		double* ptr_out = (double*) mxGetData(plhs[0]);
		ptr_out[0]=(double) 0;
	}
}

template <typename SCALAR, FDevice DEV>
void faust_rand_bsr(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	if(nlhs!=1)
		mexErrMsgTxt("rand(): 1 variable result is expected.");
	if(nrhs < 8)
		mexErrMsgTxt("rand(): wrong number of arguments (must be 8)");

	unsigned int num_rows = mxGetScalar(prhs[1]), num_cols =mxGetScalar(prhs[2]);
	unsigned int min_num_factors = (unsigned int) mxGetScalar(prhs[3]),
						max_num_factors = (unsigned int) mxGetScalar(prhs[4]),
						bnrows = (unsigned int) mxGetScalar(prhs[5]),
						bncols = (faust_unsigned_int) mxGetScalar(prhs[6]);
	float density = (float) mxGetScalar(prhs[7]);

	Faust::TransformHelper<SCALAR,DEV>* F = Faust::TransformHelper<SCALAR,DEV>::randBSRFaust(
			num_rows,
			num_cols,
			min_num_factors,
			max_num_factors,
			bnrows,
			bncols,
			density);

	if(F) //not NULL
		plhs[0]=convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(F);
	else {
		plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
		double* ptr_out = (double*) mxGetData(plhs[0]);
		ptr_out[0]=(double) 0;
	}

}
