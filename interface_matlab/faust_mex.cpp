#include "mex.h"
#include "class_handle.hpp"
#include "faust_core.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
    // Get the command string
    char cmd[256];
	if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string less than 256 characters long.");
        
    // New
    if (!strcmp("new", cmd)) {
        // Check parameters
        if (nlhs != 1)
            mexErrMsgTxt("New: One output expected.");
        // Return a handle to a new C++ instance
        plhs[0] = convertPtr2Mat<faust_core>(new faust_core);
        return;
    }
    
    // Check there is a second input, which should be the class instance handle
    if (nrhs < 2)
		mexErrMsgTxt("Second input should be a class instance handle.");
    
    // Delete
    if (!strcmp("delete", cmd)) {
        // Destroy the C++ object
        destroyObject<faust_core>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    }
    
    // Get the class instance pointer from the second input
    faust_core* core_ptr = convertMat2Ptr<faust_core>(prhs[1]);
    
    if (!strcmp("multiply", cmd)) {
        // Check parameters
        if (nlhs != 1 || nrhs != 3)
            mexErrMsgTxt("Multiply: Unexpected arguments.");
        if (mxGetNumberOfDimensions(prhs[2])!= 2 || mxGetN(prhs[2])!=1)
            mexErrMsgTxt("Multiply: Wrong number of dimensions for the vector.");
        
        
        const size_t SIZE_V = mxGetM(prhs[2]);
        const size_t SIZE_W = core_ptr->getNbRow();
        
        void* ptr_data_void = mxGetData(prhs[2]);
        faust_real* ptr_data;
       
	const mxClassID V_CLASS_ID = mxGetClassID(prhs[2]);
        if(V_CLASS_ID == mxDOUBLE_CLASS) 
        {
            double* ptr_data_tmp = static_cast<double*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[SIZE_V];
            for (int i =0 ; i<SIZE_V ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
        else if(V_CLASS_ID == mxSINGLE_CLASS)
	{
            float* ptr_data_tmp = static_cast<float*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[SIZE_V];
            for (int i =0 ; i<SIZE_V ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
        else if(V_CLASS_ID == mxINT8_CLASS)
	{
            char* ptr_data_tmp = static_cast<char*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[SIZE_V];
            for (int i =0 ; i<SIZE_V ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
        else if(V_CLASS_ID == mxUINT8_CLASS)
	{
            unsigned char* ptr_data_tmp = static_cast<unsigned char*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[SIZE_V];
            for (int i =0 ; i<SIZE_V ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
        else if(V_CLASS_ID == mxINT16_CLASS)
	{
            short* ptr_data_tmp = static_cast<short*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[SIZE_V];
            for (int i =0 ; i<SIZE_V ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
        else if (V_CLASS_ID == mxUINT16_CLASS)
	{
            unsigned short* ptr_data_tmp = static_cast<unsigned short*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[SIZE_V];
            for (int i =0 ; i<SIZE_V ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
        else if (V_CLASS_ID == mxINT32_CLASS)
	{
            int* ptr_data_tmp = static_cast<int*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[SIZE_V];
            for (int i =0 ; i<SIZE_V ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
        else if (V_CLASS_ID == mxUINT32_CLASS)
	{
            unsigned int* ptr_data_tmp = static_cast<unsigned int*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[SIZE_V];
            for (int i =0 ; i<SIZE_V ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
        else if (V_CLASS_ID == mxINT64_CLASS)
	{
            long long* ptr_data_tmp = static_cast<long long*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[SIZE_V];
            for (int i =0 ; i<SIZE_V ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
        else if (V_CLASS_ID == mxUINT64_CLASS)
	{
            unsigned long long* ptr_data_tmp = static_cast<unsigned long long*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[SIZE_V];
            for (int i =0 ; i<SIZE_V ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
            
        faust_vec v(SIZE_V, ptr_data);
        faust_vec w(SIZE_W);
        
        delete [] ptr_data;
        ptr_data = NULL;
        
        w = (*core_ptr)*v;
        
        const mwSize dims[2]={SIZE_W,1};
        if(sizeof(faust_real)==sizeof(float))
            plhs[0] = mxCreateNumericArray(2, dims, 
         mxSINGLE_CLASS, mxREAL);
        else if(sizeof(faust_real)==sizeof(double))
            plhs[0] = mxCreateNumericArray(2, dims, 
         mxDOUBLE_CLASS, mxREAL);
        else
            mexErrMsgTxt("faust_real type is neither double nor float");
        
        faust_real* ptr_out = static_cast<faust_real*> (mxGetData(plhs[0]));
        memcpy(ptr_out, w.getData(), SIZE_W*sizeof(faust_real));
        
        return;
    }
    
    
//     // Call the various class methods
//     // Train    
//     if (!strcmp("train", cmd)) {
//         // Check parameters
//         if (nlhs < 0 || nrhs < 2)
//             mexErrMsgTxt("Train: Unexpected arguments.");
//         // Call the method
//         dummy_instance->train();
//         return;
//     }
//     // Test    
//     if (!strcmp("test", cmd)) {
//         // Check parameters
//         if (nlhs < 0 || nrhs < 2)
//             mexErrMsgTxt("Test: Unexpected arguments.");
//         // Call the method
//         dummy_instance->test();
//         return;
//     }
    
    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}
