#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_multiply(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);
	if (nlhs > 1 ||  nrhs != 4)
		mexErrMsgTxt("Multiply: Unexcepted number of arguments.");

	mwSize nelem = mxGetNumberOfElements(prhs[3]);
	if (nelem != 1)
		mexErrMsgTxt("invalid char argument.");
	// boolean flag to know if the faust needs to be transposed
	bool transpose_flag = (bool) mxGetScalar(prhs[3]);

	// input matrix or vector from MATLAB
	const mxArray * inMatlabMatrix = prhs[2];

	if(mxIsSparse(inMatlabMatrix))
	{
		Faust::MatSparse<SCALAR,Cpu> spA;
		mxArray2FaustspMat<SCALAR>(inMatlabMatrix, spA);
		Faust::MatDense<SCALAR,Cpu> B;
		B = (*core_ptr).multiply(spA, transpose_flag);
		plhs[0] = FaustMat2mxArray(B);
	}
	else
	{
		const size_t nbRowA = mxGetM(inMatlabMatrix);
		const size_t nbColA = mxGetN(inMatlabMatrix);
		const size_t nbRowOp = core_ptr->getNbRow();
		const size_t nbColOp = core_ptr->getNbCol(); // getNbRow() and getNbCol is tranpose aware
		const size_t nbRowB = transpose_flag?nbColOp:nbRowOp;
		const size_t nbColB = nbColA;


		/** Check parameters **/

		//check dimension match
		if (mxGetNumberOfDimensions(inMatlabMatrix) != 2
				|| nbRowA != nbColOp && ! transpose_flag || nbRowA != nbRowOp && transpose_flag)
			mexErrMsgTxt("Multiply : Wrong number of dimensions for the input vector or matrix (third argument).");


		SCALAR* ptr_data = NULL;
		if (typeid(Faust::TransformHelper<double,DEV>) == typeid(Faust::TransformHelper<SCALAR,DEV>) && mxIsComplex(inMatlabMatrix) )
			mexErrMsgTxt("impossibile to multiply a real Faust with complex matrix");


#ifdef MX_HAS_INTERLEAVED_COMPLEX
		newMxGetData(ptr_data, inMatlabMatrix);
#else
		mxArray2Ptr(inMatlabMatrix, ptr_data);
#endif

		if(nbColA == 1)
		{
			// applying the Faust to a vector
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			SCALAR* ptr_out;
			mwSize out_dims[2] = {nbRowB, 1};
			plhs[0] = helperCreateNumMxArray<SCALAR>(out_dims);
			newMxGetData(ptr_out, plhs[0]);
			core_ptr->multiply(ptr_data, ptr_out, transpose_flag);
#else
			Faust::Vect<SCALAR,Cpu> A(nbRowA, ptr_data);
			Faust::Vect<SCALAR,Cpu> B(nbRowB);
			B = (*core_ptr).multiply(A, transpose_flag);
			plhs[0]=FaustVec2mxArray(B);
#endif


		}
		else
		{ // applying the Faust to a matrix
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			SCALAR* ptr_out;
			mwSize out_dims[2] = {nbRowB, nbColA};
			plhs[0] = helperCreateNumMxArray<SCALAR>(out_dims);
			newMxGetData(ptr_out, plhs[0]);
			core_ptr->multiply(ptr_data, nbColA, ptr_out);

#else
			Faust::MatDense<SCALAR,Cpu> A(ptr_data, nbRowA, nbColA);
			Faust::MatDense<SCALAR,Cpu> B(nbRowB, nbColA);
			B = (*core_ptr).multiply(A, transpose_flag);
			plhs[0]=FaustMat2mxArray(B);
#endif


		}

#ifndef MX_HAS_INTERLEAVED_COMPLEX
		if(ptr_data) {delete [] ptr_data ; ptr_data = NULL;}
#endif
	}

}

