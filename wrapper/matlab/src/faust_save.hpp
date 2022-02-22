#include "class_handle.hpp"
#include "faust_TransformHelper.h"


template <typename SCALAR, FDevice DEV>
void extract_filepath(const mxArray* fp_arr, char* fp, size_t fp_max_len)
{
	if(mxGetString(fp_arr, fp, fp_max_len) || strlen(fp) == 0)
	{
		mexErrMsgTxt("The filepath is not valid.");
	}
}

template <typename SCALAR, FDevice DEV>
void faust_save(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	Faust::TransformHelper<SCALAR,DEV>* core_ptr = convertMat2Ptr<Faust::TransformHelper<SCALAR,DEV> >(prhs[1]);

	const size_t MAX_PATH_SIZE = 512;
	char filepath[MAX_PATH_SIZE];
	// prhs[2] must be the filepath
	if(nrhs != 3){
		mexErrMsgTxt("The number of arguments for the save function is not valid. Must be one.");
		return;
	}
	extract_filepath<SCALAR, DEV>(prhs[2], filepath, MAX_PATH_SIZE);
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

template <typename SCALAR, FDevice DEV>
void faust_restore(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	const size_t MAX_PATH_SIZE=512;
	char filepath[MAX_PATH_SIZE];
	if(nrhs != 2){
		mexErrMsgTxt("The number of arguments for the load function is not valid. Must be one.");
		return;
	}
	extract_filepath<SCALAR, DEV>(prhs[1], filepath, MAX_PATH_SIZE);
	auto th = new Faust::TransformHelper<SCALAR,DEV>();
	th->read_from_mat_file(filepath);
	plhs[0] = convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(th);
}

template <typename SCALAR, FDevice DEV>
void faust_get_mat_file_type(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	const size_t MAX_PATH_SIZE=512;
	char filepath[MAX_PATH_SIZE];
	if(nrhs != 2){
		mexErrMsgTxt("The number of arguments for the get_mat_file_type function is not valid. Must be one.");
		return;
	}
	extract_filepath<SCALAR, DEV>(prhs[1], filepath, MAX_PATH_SIZE);
	int file_type = Faust::TransformHelper<SCALAR,DEV>::get_mat_file_type(filepath);
	plhs[0] = mxCreateDoubleScalar(file_type);
}

