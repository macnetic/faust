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
		faust_unsigned_int nbRowOp_,nbColOp_;
		const size_t nbRowOp = core_ptr->getNbRow();
		const size_t nbColOp = core_ptr->getNbCol();
		const size_t nbRowB = nbRowOp;
		const size_t nbColB = nbColA;


		/** Check parameters **/

		//check dimension match
		if (mxGetNumberOfDimensions(inMatlabMatrix) != 2
				|| nbRowA != nbColOp && (nbRowA != nbRowOp && transpose_flag))
			mexErrMsgTxt("Multiply : Wrong number of dimensions for the input vector or matrix (third argument).");


		SCALAR* ptr_data = NULL;
		if (typeid(Faust::TransformHelper<double,DEV>) == typeid(Faust::TransformHelper<SCALAR,DEV>) && mxIsComplex(inMatlabMatrix) )
			mexErrMsgTxt("impossibile to multiply a real Faust with complex matrix");

		mxArray2Ptr(inMatlabMatrix, ptr_data);

		if(nbColA == 1)
		{
			// applying the Faust to a vector
			Faust::Vect<SCALAR,Cpu> A(nbRowA, ptr_data);
			Faust::Vect<SCALAR,Cpu> B(nbRowB);
			B = (*core_ptr).multiply(A, transpose_flag);


			plhs[0]=FaustVec2mxArray(B);
		}
		else
		{ // applying the Faust to a matrix
			Faust::MatDense<SCALAR,Cpu> A(ptr_data, nbRowA, nbColA);
			Faust::MatDense<SCALAR,Cpu> B(nbRowB, nbColA);
			B = (*core_ptr).multiply(A, transpose_flag);

			plhs[0]=FaustMat2mxArray(B);

		}

		if(ptr_data) {delete [] ptr_data ; ptr_data = NULL;}
	}

}

