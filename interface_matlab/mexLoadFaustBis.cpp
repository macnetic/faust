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

if(!mxIsCell(prhs[0]))
    {
        mexErrMsgTxt("input must be a cell-array");
    }
    int nbRow = mxGetM(prhs[0]);
    int nbCol = mxGetN(prhs[0]);
	
	if (nbRow !=1)
	{
		mexErrMsgTxt("input arg must have one row");
	}
	mxArray * mxMat;
	int former_size = 0;
	faust_spmat spM;
	std::vector<faust_spmat> vec_spmat;
	for (int i = 0 ;i<nbCol;i++)
	{	
		mexPrintf("A");
		mxMat=mxGetCell(prhs[0],i);
		spM = getFaustspMat(mxMat);
		
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
	
	faust_mat fact;
	mxArray* mxfact;
	for (int i=0;i<nbCol;i++)
	{	
		mexPrintf("i : %d\n",i);
		fact=vec_spmat[i];
		mxfact=FaustMat2mxArray(fact);
		 mexCallMATLAB(0,NULL,1,&mxfact,"disp");
		
	}
	faust_core F(vec_spmat);
	faust_mat prod=F.get_product();
	mxArray* mxProd=FaustMat2mxArray(prod);
	plhs[0]=mxProd;
	
}