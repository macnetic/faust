#include "mex.h"
#include <vector>
#include <string>
#include <algorithm>

#include "faust_mat.h"
#include "faust_spmat.h"
#include "faust_core.h"
#include "faust_constant.h"

#include "tools_mex.h"

#include "class_handle.hpp"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* int nnzMax = mxGetNzmax(prhs[0]);
    int nbCol = mxGetN(prhs[0]);
    int nbRow = mxGetM(prhs[0]);
    mexPrintf("DIM (%d,%d) NNZMAX : %d\n",nbRow,nbCol,nnzMax);
    
    size_t* jc,*ir;
    double* pr;
    
    jc = (size_t *) mxCalloc(nbCol+1,sizeof(size_t));
    jc = (size_t *)mxGetJc(prhs[0]);
    ir = (size_t *) mxCalloc(nnzMax,sizeof(size_t));
    ir = (size_t *) mxGetIr(prhs[0]);
    pr = (double *) mxCalloc(nnzMax,sizeof(double));
    pr = (double *) mxGetPr(prhs[0]);
    
    for (int i=0;i<nnzMax;i++)
    {
        mexPrintf("Id_row : %d Value %f\n",ir[i],pr[i]);
    }
    
    for (int i=0;i<(nbCol+1);i++)
    {
        mexPrintf("Col_ptr: %d\n",jc[i]);
    }
    faust_spmat S(nnzMax,nbRow,nbCol,(double *)pr,(int *)ir,(int *)jc); */
	//mexPrintf("abc\n");
	if(nlhs!=1)
		mexErrMsgTxt("mexLoadFaust must have 1 output.");
	if(nrhs!=1)
		mexErrMsgTxt("mexLoadFaust must have 1 input.");


	//int nbRow,nbCol;
	if(!mxIsCell(prhs[0]))
    {
        mexErrMsgTxt("input must be a cell-array");
    }
	std::vector<faust_spmat> vec_spmat;
	mwSize nb_element = mxGetNumberOfElements(prhs[0]);

	if (nb_element == 0)
		mexWarnMsgTxt("Empty cell array.");
    else 
	{
			mxArray * mxMat;	
			for (mwSize i=0;i<nb_element;i++)
			{	
				mxMat=mxGetCell(prhs[0],i);
				addSpmat(mxMat,vec_spmat);
			}
	}
			
	
	faust_core* F = new faust_core(vec_spmat); 
	plhs[0]=convertPtr2Mat<faust_core>(F);
    
    
    
}




