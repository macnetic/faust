#include "mex.h"
#include "class_handle.hpp"
#include "faust_core.h"
#include "mexFaustMat.h"


// prhs[0] : name of command : 
//    "delete" to delete the faust_core object dynamically allocated previously
//    "multiply" to multiply the faust_core object by a vector or a matrix

// prhs[1] : address of the faust_core object dynamically allocated previously

// prhs[2] (only necessary if prhs[0] matches "multiply") : vector or matrix A to multiply by the faust_core object



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
	// Get the command string
	char cmd[256];
	if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string less than 256 characters long.");
	// Check there is a second input, which should be the class instance handle
	if (nrhs < 2)
		mexErrMsgTxt("Second input should be a class instance handle.");
 
	// New
	if (!strcmp("new", cmd)) 
	{
		if(nlhs!=1)
			mexErrMsgTxt("1 output is expected.");
		if(nrhs!=2)
			mexErrMsgTxt("2 inputs are expected.");


		int nbRow,nbCol;
		if(!mxIsCell(prhs[1]))
			mexErrMsgTxt("input must be a cell-array");

		std::vector<faust_spmat> vec_spmat;
		mwSize nb_element = mxGetNumberOfElements(prhs[1]);
		/*if (nb_element == 0)
			mexWarnMsgTxt("Empty cell array.");
	        else if (!mxIsSparse(mxGetCell(prhs[1],0)))
		{
			mexPrintf("Dense\n");	
			loadDenseFaust(prhs[1],vec_spmat);
		}
		else
		{	
			mexPrintf("Sparse\n");
			loadSpFaust(prhs[1],vec_spmat);
		}*/
	if (nb_element == 0)
		mexWarnMsgTxt("Empty cell array.");
    else 
	{
			mxArray * mxMat;	
			for (int i=0;i<nb_element;i++)
			{	
				mxMat=mxGetCell(prhs[1],i);
				addSpmat(mxMat,vec_spmat);
			}
	}
	
	
		faust_core* F = new faust_core(vec_spmat); 
		plhs[0]=convertPtr2Mat<faust_core>(F);
    
		return;
	}
    


    
	// Delete
	if (!strcmp("delete", cmd))
	{
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
       
        const size_t SIZE_A1 = mxGetM(prhs[2]);
        const size_t SIZE_A2 = mxGetN(prhs[2]);

        const size_t SIZE_B1 = core_ptr->getNbRow(); 
        const size_t SIZE_B2 = SIZE_A2; 

        
        // Check parameters
        if (nlhs != 1 || nrhs != 3)
            mexErrMsgTxt("Multiply: Unexpected arguments.");
        if (mxGetNumberOfDimensions(prhs[2]) != 2 
                || SIZE_A1 != core_ptr->getNbCol() )
            mexErrMsgTxt("Multiply: Wrong number of dimensions for the input vector or matrix (third argument).");
        
             
        void* ptr_data_void = mxGetData(prhs[2]);
        faust_real* ptr_data;
       
	const mxClassID V_CLASS_ID = mxGetClassID(prhs[2]);
	const size_t NB_ELEMENTS = mxGetNumberOfElements(prhs[2]);
        if(V_CLASS_ID == mxDOUBLE_CLASS) 
        {
            double* ptr_data_tmp = static_cast<double*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[NB_ELEMENTS];
            for (int i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
        else if(V_CLASS_ID == mxSINGLE_CLASS)
	{
            float* ptr_data_tmp = static_cast<float*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[NB_ELEMENTS];
            for (int i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
        else if(V_CLASS_ID == mxINT8_CLASS)
	{
            char* ptr_data_tmp = static_cast<char*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[NB_ELEMENTS];
            for (int i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
        else if(V_CLASS_ID == mxUINT8_CLASS)
	{
            unsigned char* ptr_data_tmp = static_cast<unsigned char*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[NB_ELEMENTS];
            for (int i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
        else if(V_CLASS_ID == mxINT16_CLASS)
	{
            short* ptr_data_tmp = static_cast<short*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[NB_ELEMENTS];
            for (int i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
        else if (V_CLASS_ID == mxUINT16_CLASS)
	{
            unsigned short* ptr_data_tmp = static_cast<unsigned short*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[NB_ELEMENTS];
            for (int i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
        else if (V_CLASS_ID == mxINT32_CLASS)
	{
            int* ptr_data_tmp = static_cast<int*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[NB_ELEMENTS];
            for (int i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
        else if (V_CLASS_ID == mxUINT32_CLASS)
	{
            unsigned int* ptr_data_tmp = static_cast<unsigned int*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[NB_ELEMENTS];
            for (int i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
        else if (V_CLASS_ID == mxINT64_CLASS)
	{
            long long* ptr_data_tmp = static_cast<long long*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[NB_ELEMENTS];
            for (int i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
        else if (V_CLASS_ID == mxUINT64_CLASS)
	{
            unsigned long long* ptr_data_tmp = static_cast<unsigned long long*> (mxGetData(prhs[2]));
            ptr_data = new faust_real[NB_ELEMENTS];
            for (int i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<faust_real> (ptr_data_tmp[i]);
	}
           
		
 
	// Si prhs[2] est un vecteur
	if(SIZE_A2 == 1) 
	{
        	faust_vec A(SIZE_A1, ptr_data);
        	faust_vec B(SIZE_B1);
		B = (*core_ptr)*A;
		
		const mwSize dims[2]={SIZE_B1,SIZE_B2};
		if(sizeof(faust_real)==sizeof(float))
			plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
		else if(sizeof(faust_real)==sizeof(double))
			plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
		else
			mexErrMsgTxt("faust_real type is neither double nor float");
        
		faust_real* ptr_out = static_cast<faust_real*> (mxGetData(plhs[0]));
		memcpy(ptr_out, B.getData(), SIZE_B1*SIZE_B2*sizeof(faust_real));
	}
	// Si prhs[2] est une matrice
	else
	{
        	faust_mat A(ptr_data, SIZE_A1, SIZE_A2);
		faust_mat B(SIZE_B1, SIZE_A2);
		B = (*core_ptr)*A;
		
		const mwSize dims[2]={SIZE_B1,SIZE_B2};
		if(sizeof(faust_real)==sizeof(float))
			plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
		else if(sizeof(faust_real)==sizeof(double))
			plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
		else
			mexErrMsgTxt("faust_real type is neither double nor float");
        
		faust_real* ptr_out = static_cast<faust_real*> (mxGetData(plhs[0]));
		memcpy(ptr_out, B.getData(), SIZE_B1*SIZE_B2*sizeof(faust_real));
	}
        	delete [] ptr_data ; ptr_data = NULL;
        
       
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
