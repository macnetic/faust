#include <mexFaustMat.h>

faust_mat getFaustMat(mxArray* Mat_array)
{
    int  nbRow,nbCol;
    double* MatPtr;
    const mwSize *dimsMat;
    dimsMat = mxGetDimensions(Mat_array);
	nbRow = (int) dimsMat[0];
	nbCol = (int) dimsMat[1];
    MatPtr = mxGetPr(Mat_array);
    
    faust_mat Mat(nbRow,nbCol);
    memcpy(Mat.getData(),MatPtr,nbRow*nbCol*sizeof(double));
    return Mat;
}

mxArray*  FaustMat2mxArray(faust_mat M)
{		
		mxArray * mxMat;
		double * mat_ptr;
		int row,col;
		row = M.getNbRow();
        col = M.getNbCol();
        mxMat = mxCreateDoubleMatrix(row,col,mxREAL);
        mat_ptr = (double *) mxCalloc(row*col,sizeof(double));
        memcpy(mat_ptr,M.getData(),row*col*sizeof(double));
        mxSetM(mxMat, row);
        mxSetN(mxMat, col);
        mxSetPr(mxMat, mat_ptr);
		mxArray * rhs[1];
		rhs[0]=mxMat;
        mexCallMATLAB(0,NULL,1,rhs, "disp");
		return mxMat;
}