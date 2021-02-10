#include "class_handle.hpp"
#include "faust_TransformHelper.h"
template <typename SCALAR, FDevice DEV>
void faust_save(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);

	const size_t MAX_PATH_SIZE=512;
	char filepath[MAX_PATH_SIZE];
	// prhs[2] must be the filepath
    if(nrhs != 3){
            mexErrMsgTxt("The number of arguments for the save function is not valid. Must be one.");
            return;
        }
        if(mxGetString(prhs[2], filepath, sizeof(filepath)) || strlen(filepath) == 0){
            mexErrMsgTxt("The filepath is not valid.");
            return;
        }
        //	printf("save: filepath = %s, nprhs = %d\n", filepath, nprhs);
		try
		{
			core_ptr->save_mat_file(filepath/*, mxGetScalar(prhs[1])*/);
		}
		catch(exception& e)
		{
			mexErrMsgTxt("Failed to save the file.");
		}

}

