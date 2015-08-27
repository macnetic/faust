#include <mex.h>
#include <vector>
#include <string>
#include <algorithm>

#include <faust_mat.h>
#include <faust_spmat.h>
#include <faust_core.h>
#include <faust_constant.h>

#include <mexFaustMat.h>




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
	int nbRow,nbCol;
	
	if(!mxIsCell(prhs[0]))
    {
        mexErrMsgTxt("cons must be a cell-array");
    }
    nbRow = mxGetM(prhs[0]);
    nbCol = mxGetN(prhs[0]);
	
	if (nbRow !=1)
	{
		mexErrMsgTxt("input arg must have one row");
	}
	mxArray * mxMat;
	faust_mat M;
	int former_size = 0;
	faust_spmat spM;
	std::vector<faust_spmat> vec_spmat;
	for (int i = 0 ;i<nbCol;i++)
	{	
		mxMat=mxGetCell(prhs[0],i);
		M = getFaustMat(mxMat);
		
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
	
	faust_core F(vec_spmat); 
	faust_mat prod=F.get_product();
	mxArray* mxProd=FaustMat2mxArray(prod);
	plhs[0]=mxProd;
    
    
    
}