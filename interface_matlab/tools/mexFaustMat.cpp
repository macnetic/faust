#include <mexFaustMat.h>
#include <faust_constraint_real.h>
#include <faust_constraint_mat.h>
#include <faust_constraint_int.h>
faust_mat getFaustMat(mxArray* Mat_array)
{
    int  nbRow,nbCol;
    double* MatPtr;
    const mwSize *dimsMat;
    dimsMat = mxGetDimensions(Mat_array);
	nbRow = (int) dimsMat[0];
	nbCol = (int) dimsMat[1];
	if ((nbRow == 0) || (nbCol == 0))
	{
		mexErrMsgIdAndTxt("mexFaustMat.h:getFaustMat",
           "empty matrix");
	}
	if (mxIsSparse(Mat_array))
	{	
		//mexErrMsgTxt("sparse matrix entry instead of dense matrix");
		mexErrMsgIdAndTxt("a","a sparse matrix entry instead of dense matrix");
	}
    MatPtr = mxGetPr(Mat_array);
    
    faust_mat Mat(nbRow,nbCol);
    memcpy(Mat.getData(),MatPtr,nbRow*nbCol*sizeof(double));
    return Mat;
}


faust_spmat getFaustspMat(mxArray* spMat_array)
{	
	if (!mxIsSparse(spMat_array))
	{
		mexErrMsgIdAndTxt("mexFaustMat.h:getFaustspMat",
           "input array must be sparse");
	}		
	int nnzMax = mxGetNzmax(spMat_array);
    int nbCol = mxGetN(spMat_array);
    int nbRow = mxGetM(spMat_array);
    mexPrintf("DIM (%d,%d) NNZMAX : %d\n",nbRow,nbCol,nnzMax);
    
    size_t* jc,*ir;
    double* pr;
    
    jc = (size_t *) mxCalloc(nbCol+1,sizeof(size_t));
    jc = (size_t *)mxGetJc(spMat_array);
    ir = (size_t *) mxCalloc(nnzMax,sizeof(size_t));
    ir = (size_t *) mxGetIr(spMat_array);
    pr = (double *) mxCalloc(nnzMax,sizeof(double));
    pr = (double *) mxGetPr(spMat_array);
    
    for (int i=0;i<nnzMax;i++)
    {
        mexPrintf("Id_row : %d Value %f\n",ir[i],pr[i]);
    }
    
    for (int i=0;i<(nbCol+1);i++)
    {
        mexPrintf("Col_ptr: %d\n",jc[i]);
    }

    faust_spmat S(nnzMax,nbRow,nbCol,pr,ir,jc); 
	/*faust_mat A=S;
	mxArray*   mxA=FaustMat2mxArray(A);
	mexPrintf("INSIDE\n");
	 mexCallMATLAB(0,NULL,1,&mxA,"disp");*/
	return S;
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
		//rhs[0]=mxMat;
        //mexCallMATLAB(0,NULL,1,rhs, "disp");
		return mxMat;
}





void setCellFacts(mxArray **  cellFacts,std::vector<faust_mat> facts)
{   
    int rowFact,colFact;
    int nb_fact = facts.size();
    faust_mat mat;
    //mexPrintf("size CellFacts nb fact %d ",nb_fact);
    (*cellFacts) = mxCreateCellMatrix(1,nb_fact);
    mxArray * mxMat;
    double* mat_ptr;
    mxMat = mxCreateDoubleMatrix(0,0,mxREAL);
    
    for (size_t k = 0; k < nb_fact; k++)
    {
        mat = facts[k];
        rowFact = mat.getNbRow();
        colFact = mat.getNbCol();
        
        mat_ptr = (double *) mxCalloc(rowFact*colFact,sizeof(double));
        //mat_ptr = mxGetPr(mxMat);
        memcpy(mat_ptr,mat.getData(),rowFact*colFact*sizeof(double));
        mxSetM(mxMat, rowFact);
        mxSetN(mxMat, colFact);
        mxSetPr(mxMat, mat_ptr);
        
//          mexPrintf("DIM : %d %d\n",rowFact,colFact);
//         for (int i=1;i<rowFact;i++)
//         {
//             for (int j=1;j<colFact;j++)
//             {
//                 mexPrintf("%f ",mat_ptr[i+j*rowFact]);
//             }
//             mexPrintf("\n");
//         }
//          mexPrintf("\n");
//          mexPrintf("\n");
         //rhs[0]=mxMat;
         //mexCallMATLAB(0,NULL,1,rhs, "disp");
        mxSetCell((*cellFacts),k,mxDuplicateArray(mxMat));
    }
    
//     mexPrintf("Display Cell \n");
//     for (size_t l = 0; l < mxGetNumberOfElements((*cellFacts)); l++)
//     {    
//         mexPrintf("facts %d \n",l);
//          mexPrintf("\n");
//          mxArray * rhs[1]; 
//          rhs[0]=mxGetCell((*cellFacts),l);
//          mexCallMATLAB(0,NULL,1,rhs, "disp");
//     }
    
}

void setVectorFaustMat(std::vector<faust_mat> &vecMat,mxArray *Cells)
{	
	mxArray* mxMat;
	mwSize nb_fact = mxGetNumberOfElements (Cells);
	faust_mat mat;
	vecMat.resize(0);
	for (int i=0;i<nb_fact;i++)
	{
		mxMat=mxGetCell(Cells,i);
		mat=getFaustMat(mxMat);
		vecMat.push_back(mat);	
	}
}


void getConstraint(std::vector<const faust_constraint_generic*> & consS,mxArray* mxCons)
{     
    mwSize bufCharLen,nbRowCons,nbColCons;
    int status;        
    char * consName;
    double paramCons;
    mxArray * mxConsParams; 
    
    mxConsParams=mxGetCell(mxCons,0);
    bufCharLen = mxGetNumberOfElements(mxConsParams)+1;
    consName = (char *) mxCalloc(bufCharLen,sizeof(char));
    status = mxGetString(mxConsParams,consName,bufCharLen);
    paramCons = mxGetScalar(mxConsParams);
    mxConsParams=mxGetCell(mxCons,2);
    nbRowCons   = (int) mxGetScalar(mxConsParams);
    mxConsParams=mxGetCell(mxCons,3);
    nbColCons = (int) mxGetScalar(mxConsParams);
    
    bool is_const_int =((strcmp(consName,"sp") == 0) || (strcmp(consName,"sppos")==0));
	is_const_int = ((is_const_int) || ((strcmp(consName,"spcol") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"splin") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"splincol") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"lOpen") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"l1pen") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"splin") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"wav") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"blkdiag") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"splin_test") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"supp") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"normlin") == 0)));
    
    int const_type = -1;
	if (is_const_int)
	{
		const_type = 0;
	}
	/*if (is_const_real)
	{
		const_type = 1;
	}
	if (is_const_mat)
	{
		const_type = 2;
	}
     */
    int intParameter;
    switch(const_type)
	{
        case 0:
            mxConsParams=mxGetCell(mxCons,1);
             intParameter = (int) std::floor(mxGetScalar(mxConsParams)+0.5);
             mexPrintf("NAME  %s PARAMS %d DIMS : (%d,%d)\n",consName,intParameter,nbRowCons,nbColCons);
            faust_constraint_name consNameType;
            if (strcmp(consName,"sp") == 0)
            {
                consNameType = CONSTRAINT_NAME_SP;
            }
            if (strcmp(consName,"spcol") == 0)
            {
                consNameType = CONSTRAINT_NAME_SPCOL;
            }
            if (strcmp(consName,"splin") == 0)
            {
                consNameType = CONSTRAINT_NAME_SPLIN;
            }
            
            consS.push_back((new faust_constraint_int(consNameType,intParameter,nbRowCons,nbColCons)));
            break;
            
            
            default :
                mexErrMsgTxt("Unknown constraint name ");
			break;
    }

   
   
    
}










