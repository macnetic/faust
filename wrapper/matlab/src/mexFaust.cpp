#include "mex.h"
#include "class_handle.hpp"
#include "faust_Transform.h"
#include "tools_mex.h"
#include "faust_MatDense.h"
#include <stdexcept>
#include "faust_constant.h"
#include "faust_Timer.h"
// prhs[0] : name of command :
//    "delete" to delete the Faust::Transform<FFPP> object dynamically allocated previously
//    "multiply" to multiply the Faust::Transform<FFPP> object by a vector or a matrix

// prhs[1] : address of the Faust::Transform<FFPP> object dynamically allocated previously

// prhs[2] (only necessary if prhs[0] matches "multiply") : vector or matrix A to multiply by the Faust::Transform<FFPP> object



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	#ifdef FAUST_VERBOSE
		if (typeid(FFPP) == typeid(float))
		{
			std::cout<<"FFPP == float"<<std::endl;
		}

		if (typeid(FFPP) == typeid(double))
		{
			std::cout<<"FFPP == double"<<std::endl;
		}
	#endif
	try{
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
        if((nrhs<2) || (nrhs>3))
			mexErrMsgTxt("1 or 2 inputs are expected.");
        

		if(!mxIsCell(prhs[1]))
			mexErrMsgTxt("1st arg input must be a cell-array");
        
		std::vector<Faust::MatSparse<FFPP,Cpu> > vec_spmat;
		mwSize nb_element = mxGetNumberOfElements(prhs[1]);
		if (nb_element != 0)
		{
			mxArray * mxMat;
			for (int i=0;i<nb_element;i++)
			{
				mxMat=mxGetCell(prhs[1],i);
				addSpmat<FFPP>(mxMat,vec_spmat);
			}
		}
        
        FFPP lambda = 1.0;
		if (nrhs > 2)			
			lambda = (FFPP) mxGetScalar(prhs[2]);

		Faust::Transform<FFPP,Cpu>* F = new Faust::Transform<FFPP,Cpu>(vec_spmat,lambda);
		plhs[0]=convertPtr2Mat<Faust::Transform<FFPP,Cpu> >(F);

		return;
	}




	// Delete
	if (!strcmp("delete", cmd))
	{
		// Destroy the C++ object
		destroyObject<Faust::Transform<FFPP,Cpu> >(prhs[1]);
		// Warn if other commands were ignored
		if (nlhs != 0 || nrhs != 2)
			mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
		return;
	}


	// Get the class instance pointer from the second input
	Faust::Transform<FFPP,Cpu>* core_ptr = convertMat2Ptr<Faust::Transform<FFPP,Cpu> >(prhs[1]);

	
	if (!strcmp("size",cmd))
	{
		const size_t SIZE_B1 = core_ptr->getNbRow();
        	const size_t SIZE_B2 = core_ptr->getNbCol();
		const mwSize dims[2]={1,2};
		plhs[0]=mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
		double* ptr_out = (double*) mxGetData(plhs[0]);
		ptr_out[0]=(double) SIZE_B1;
		ptr_out[1]=(double) SIZE_B2;
		return;
	}

	// get_fact : return the id factor of the faust	
	if (!strcmp("get_fact", cmd))
    	{
		if (nlhs > 1 || nrhs != 3)
		{
			mexErrMsgTxt("get_fact : incorrect number of arguments.");
		}
		int id_fact = (faust_unsigned_int) (mxGetScalar(prhs[2])-1);
		Faust::MatSparse<FFPP,Cpu> const sp_factor = core_ptr->get_fact(id_fact);
		Faust::MatDense<FFPP,Cpu> const dense_factor(sp_factor);
		plhs[0]=FaustMat2mxArray(dense_factor);
		return;
		

	}


	if (!strcmp("get_nb_factor", cmd))
    	{
		if (nlhs != 1 || nrhs != 2)
		{
			mexErrMsgTxt("get_nb_factor : incorrect number of arguments.");
		}
		
		faust_unsigned_int nb_fact = core_ptr->size();
		plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
		double* ptr_out = (double*) mxGetData(plhs[0]);
		ptr_out[0]=(double) nb_fact;
		return;
		

	}


    

    	if (!strcmp("get_product",cmd))
	{	
		
		if (nlhs != 1 || nrhs != 3)
			mexErrMsgTxt("get_product: Unexpected arguments");
		if(core_ptr->size() == 0)
			mexErrMsgTxt("get_product : empty faust core");
		
		bool transpose_flag = (bool) mxGetScalar(prhs[2]);

		char op;
		if (transpose_flag)
			op='T';
		else
			op='N';
	
		faust_unsigned_int nbRowOp,nbColOp;
		(*core_ptr).setOp(op,nbRowOp,nbColOp);
		const size_t SIZE_B1 = nbRowOp;
        	const size_t SIZE_B2 = nbColOp;
		Faust::MatDense<FFPP,Cpu> prod=core_ptr->get_product(op);
		const mwSize dims[2]={SIZE_B1,SIZE_B2};
		if(sizeof(FFPP)==sizeof(float))
			plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
		else if(sizeof(FFPP)==sizeof(double))
			plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
		else
			mexErrMsgTxt("FFPP type is neither double nor float");

		FFPP* ptr_out = static_cast<FFPP*> (mxGetData(plhs[0]));
		memcpy(ptr_out, prod.getData(), SIZE_B1*SIZE_B2*sizeof(FFPP));

		return;
	}

	if (!strcmp("copy",cmd))
	{



		if (nlhs > 1)
			mexErrMsgTxt("transpose : too much output argument");
		else
		{
			//Faust::Timer t1,t2,t3;
			//t1.start();			
			Faust::Transform<FFPP,Cpu>* F = new Faust::Transform<FFPP,Cpu>((*core_ptr));			
			//t1.stop();			
			//t2.start();			
			//(*F).transpose();
			//t2.stop();
			//t3.start();
			plhs[0]=convertPtr2Mat<Faust::Transform<FFPP,Cpu> >(F);
			//t3.stop();
			//std::cout<<"t1 new : "<<t1.get_time()<<std::endl;
			//std::cout<<"t2 transpose : "<<t2.get_time()<<std::endl;
			//std::cout<<"t3 convertPtrMat : "<<t3.get_time()<<std::endl;
				
		}
		return;

	}


    if (!strcmp("multiply", cmd)) {
	
	if (nlhs > 1 ||  nrhs != 4)
            mexErrMsgTxt("Multiply: Unexpected arguments.");

	mwSize nelem = mxGetNumberOfElements(prhs[3]);
	if (nelem != 1)
		mexErrMsgTxt("invalid char argument.");
	
	// boolean flag to know if the faust si transposed
	bool transpose_flag = (bool) mxGetScalar(prhs[3]);
	char op;
	if (transpose_flag)
		op='T';
	else
		op='N';


			
        const size_t nbRowA = mxGetM(prhs[2]);
        const size_t nbColA = mxGetN(prhs[2]);
	faust_unsigned_int nbRowOp_,nbColOp_;
	(*core_ptr).setOp(op,nbRowOp_,nbColOp_);
	const size_t nbRowOp = nbRowOp_;
	const size_t nbColOp = nbColOp_;			
        const size_t nbRowB = nbRowOp;
        const size_t nbColB = nbColA;

	

        // Check parameters
        
	
        if (mxGetNumberOfDimensions(prhs[2]) != 2
                || nbRowA != nbColOp )
            mexErrMsgTxt("Multiply: Wrong number of dimensions for the input vector or matrix (third argument).");


        FFPP* ptr_data = NULL;

	const mxClassID V_CLASS_ID = mxGetClassID(prhs[2]);
	const size_t NB_ELEMENTS = mxGetNumberOfElements(prhs[2]);
        if(V_CLASS_ID == mxDOUBLE_CLASS)
        {
            double* ptr_data_tmp = static_cast<double*> (mxGetData(prhs[2]));
            ptr_data = new FFPP[NB_ELEMENTS];
            for (size_t i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<FFPP> (ptr_data_tmp[i]);
	}
        else if(V_CLASS_ID == mxSINGLE_CLASS)
	{
            float* ptr_data_tmp = static_cast<float*> (mxGetData(prhs[2]));
            ptr_data = new FFPP[NB_ELEMENTS];
            for (size_t i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<FFPP> (ptr_data_tmp[i]);
	}
        else if(V_CLASS_ID == mxINT8_CLASS)
	{
            char* ptr_data_tmp = static_cast<char*> (mxGetData(prhs[2]));
            ptr_data = new FFPP[NB_ELEMENTS];
            for (size_t i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<FFPP> (ptr_data_tmp[i]);
	}
        else if(V_CLASS_ID == mxUINT8_CLASS)
	{
            unsigned char* ptr_data_tmp = static_cast<unsigned char*> (mxGetData(prhs[2]));
            ptr_data = new FFPP[NB_ELEMENTS];
            for (size_t i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<FFPP> (ptr_data_tmp[i]);
	}
        else if(V_CLASS_ID == mxINT16_CLASS)
	{
            short* ptr_data_tmp = static_cast<short*> (mxGetData(prhs[2]));
            ptr_data = new FFPP[NB_ELEMENTS];
            for (size_t i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<FFPP> (ptr_data_tmp[i]);
	}
        else if (V_CLASS_ID == mxUINT16_CLASS)
	{
            unsigned short* ptr_data_tmp = static_cast<unsigned short*> (mxGetData(prhs[2]));
            ptr_data = new FFPP[NB_ELEMENTS];
            for (size_t i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<FFPP> (ptr_data_tmp[i]);
	}
        else if (V_CLASS_ID == mxINT32_CLASS)
	{
            int* ptr_data_tmp = static_cast<int*> (mxGetData(prhs[2]));
            ptr_data = new FFPP[NB_ELEMENTS];
            for (size_t i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<FFPP> (ptr_data_tmp[i]);
	}
        else if (V_CLASS_ID == mxUINT32_CLASS)
	{
            unsigned int* ptr_data_tmp = static_cast<unsigned int*> (mxGetData(prhs[2]));
            ptr_data = new FFPP[NB_ELEMENTS];
            for (size_t i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<FFPP> (ptr_data_tmp[i]);
	}
        else if (V_CLASS_ID == mxINT64_CLASS)
	{
            long long* ptr_data_tmp = static_cast<long long*> (mxGetData(prhs[2]));
            ptr_data = new FFPP[NB_ELEMENTS];
            for (size_t i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<FFPP> (ptr_data_tmp[i]);
	}
        else if (V_CLASS_ID == mxUINT64_CLASS)
	{
            unsigned long long* ptr_data_tmp = static_cast<unsigned long long*> (mxGetData(prhs[2]));
            ptr_data = new FFPP[NB_ELEMENTS];
            for (size_t i =0 ; i<NB_ELEMENTS ; i++)
                ptr_data[i] = static_cast<FFPP> (ptr_data_tmp[i]);
	}
	else
            mexErrMsgTxt("Unknown matlab type.");


	// Si prhs[2] est un vecteur
	if(nbColA == 1)
	{
        Faust::Vect<FFPP,Cpu> A(nbRowA, ptr_data);
        Faust::Vect<FFPP,Cpu> B(nbRowB);
	//NB        
	//B = (*core_ptr)*A;
	B = (*core_ptr).multiply(A,op);
	
		const mwSize dims[2]={nbRowB,nbColB};
		if(sizeof(FFPP)==sizeof(float))
			plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
		else if(sizeof(FFPP)==sizeof(double))
			plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
		else
			mexErrMsgTxt("FFPP type is neither double nor float");

		FFPP* ptr_out = static_cast<FFPP*> (mxGetData(plhs[0]));
		memcpy(ptr_out, B.getData(), nbRowB*nbColB*sizeof(FFPP));
	}
	// Si prhs[2] est une matrice
	else
	{
        	Faust::MatDense<FFPP,Cpu> A(ptr_data, nbRowA, nbColA);
		Faust::MatDense<FFPP,Cpu> B(nbRowB, nbColA);
		B = (*core_ptr).multiply(A,op);

		const mwSize dims[2]={nbRowB,nbColB};
		if(sizeof(FFPP)==sizeof(float))
			plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
		else if(sizeof(FFPP)==sizeof(double))
			plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
		else
			mexErrMsgTxt("FFPP type is neither double nor float");

		FFPP* ptr_out = static_cast<FFPP*> (mxGetData(plhs[0]));
		memcpy(ptr_out, B.getData(), nbRowB*nbColB*sizeof(FFPP));
	}
	if(ptr_data) {delete [] ptr_data ; ptr_data = NULL;}


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
	}catch (const std::exception& e)
	{
		 mexErrMsgTxt(e.what());
	}
}
