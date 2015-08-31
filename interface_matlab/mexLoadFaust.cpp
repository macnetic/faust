#include <mex.h>
#include <vector>
#include <string>
#include <algorithm>

#include <faust_mat.h>
#include <faust_spmat.h>
#include <faust_core.h>
#include <faust_constant.h>

#include <mexFaustMat.h>

#include "class_handle.hpp"
void loadDenseFaust( const mxArray * Cells,std::vector<faust_spmat> &vec_spmat);
void loadSpFaust(const mxArray * Cells,std::vector<faust_spmat> &vec_spmat);


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
	mexPrintf("abc\n");
	if(nlhs!=1)
		mexErrMsgTxt("mexLoadFaust must have 1 output.");
	if(nrhs!=1)
		mexErrMsgTxt("mexLoadFaust must have 1 input.");


	int nbRow,nbCol;
	if(!mxIsCell(prhs[0]))
    {
        mexErrMsgTxt("input must be a cell-array");
    }
	std::vector<faust_spmat> vec_spmat;
	mwSize nb_element = mxGetNumberOfElements(prhs[0]);
	if (nb_element == 0)
	{
		
	}else
		if (!mxIsSparse(mxGetCell(prhs[0],0)))
		{
			 mexPrintf("Dense\n");	
			loadDenseFaust(prhs[0],vec_spmat);
		}else
		{	
			mexPrintf("Sparse\n");
			loadSpFaust(prhs[0],vec_spmat);
		}
			
	
	faust_core* F = new faust_core(vec_spmat); 
	plhs[0]=convertPtr2Mat<faust_core>(F);
    
    
    
}




void loadDenseFaust(const mxArray * Cells,std::vector<faust_spmat> &vec_spmat) 
{
	if(!mxIsCell(Cells))
    {
        mexErrMsgTxt("cons must be a cell-array");
    }
    mwSize nbRow = mxGetM(Cells);
    mwSize nbCol = mxGetN(Cells);
	
	if (nbRow !=1)
	{
		mexErrMsgTxt("input arg must have one row");
	}
	mxArray * mxMat;
	faust_mat M;
	int former_size = 0;
	faust_spmat spM;
	vec_spmat.resize(0);
	for (int i = 0 ;i<nbCol;i++)
	{	
		mxMat=mxGetCell(Cells,i);
		getFaustMat(mxMat,M);
		
		if (i != 0)
		{
			if (M.getNbRow() != former_size)
			{
				mexErrMsgTxt("invalid dimensions of the matrix");
			}	
		}
		former_size = M.getNbCol();
		spM = M;
		vec_spmat.push_back(spM);	
	}
}


void loadSpFaust(const mxArray * Cells,std::vector<faust_spmat> &vec_spmat)
{
	if(!mxIsCell(Cells))
    {
        mexErrMsgTxt("Cells must be a cell-array");
    }
    mwSize nbRow = mxGetM(Cells);
    mwSize nbCol = mxGetN(Cells);
	
	if (nbRow !=1)
	{
		mexErrMsgTxt("input arg must have one row");
	}
	mxArray * mxMat;
	int former_size = 0;
	faust_spmat spM;
	vec_spmat.resize(0);
	for (int i = 0 ;i<nbCol;i++)
	{	
		mexPrintf("A");
		mxMat=mxGetCell(Cells,i);
		getFaustspMat(mxMat,spM);
		
		if (i != 0)
		{
			if (spM.getNbRow() != former_size)
			{
				mexErrMsgTxt("invalid dimensions of the matrix");
			}	
		}
		former_size = spM.getNbCol();
		vec_spmat.push_back(spM);	
	}
	
	/*faust_mat fact;
	mxArray* mxfact;
	for (int i=0;i<nbCol;i++)
	{	
		mexPrintf("i : %d\n",i);
		fact=vec_spmat[i];
		mxfact=FaustMat2mxArray(fact);
		 mexCallMATLAB(0,NULL,1,&mxfact,"disp");
		
	}*/
}
