#include <mexFaustMat.h>
#include <faust_constraint_real.h>
#include <faust_constraint_mat.h>
#include <faust_constraint_int.h>

void getFaustMat(const mxArray* Mat_array,faust_mat & Mat)
{
    int  nbRow,nbCol;
   
    const mwSize *dimsMat;
    dimsMat = mxGetDimensions(Mat_array);
    nbRow = (int) dimsMat[0];
    nbCol = (int) dimsMat[1];
    if ((nbRow == 0) || (nbCol == 0))
        mexErrMsgIdAndTxt("mexFaustMat.h:getFaustMat", "empty matrix");
    if (mxIsSparse(Mat_array))
    {	
        //mexErrMsgTxt("sparse matrix entry instead of dense matrix");
        mexErrMsgIdAndTxt("a","a sparse matrix entry instead of dense matrix");
    }
	const mxClassID V_CLASS_ID = mxGetClassID(Mat_array);
	 faust_real* MatPtr; 
	if (((V_CLASS_ID == mxDOUBLE_CLASS) && (sizeof(double) == sizeof(faust_real))) || ((V_CLASS_ID == mxSINGLE_CLASS) && (sizeof(float) == sizeof(faust_real))))
	{
		MatPtr = (faust_real*) mxGetPr(Mat_array);
	}else
	{	
		if (V_CLASS_ID == mxDOUBLE_CLASS) 
		{	
			MatPtr = (faust_real*) mxCalloc(nbRow*nbCol,sizeof(faust_real));
			double* MatPtrDouble =(double*) mxGetPr(Mat_array);
			for (int i=0;i<nbRow*nbCol;i++)
				MatPtr[i] = (faust_real) MatPtrDouble[i];
		}
		else if (V_CLASS_ID == mxSINGLE_CLASS)
		{		
			MatPtr = (faust_real*) mxCalloc(nbRow*nbCol,sizeof(faust_real));
			float* MatPtrSingle= (float*) (mxGetData(Mat_array));
			for (int i=0;i<nbRow*nbCol;i++)
				MatPtr[i] = (faust_real) MatPtrSingle[i];
		
		
		}else
		{
		 mexErrMsgTxt("getFaustMat :input matrix format must be single or double");
		}
	}
	
     Mat.resize(nbRow,nbCol);
	
    memcpy(Mat.getData(),MatPtr,nbRow*nbCol*sizeof(faust_real));
	if (((V_CLASS_ID == mxDOUBLE_CLASS) && (sizeof(double) != sizeof(faust_real))) || ((V_CLASS_ID == mxSINGLE_CLASS) && (sizeof(float) != sizeof(faust_real))))
	{
		mxFree(MatPtr);
	}	
	

    
}


void getFaustspMat(const mxArray* spMat_array,faust_spmat & S)
{	
	if (!mxIsSparse(spMat_array))
	{
		mexErrMsgIdAndTxt("mexFaustMat.h:getFaustspMat",
           "input array must be sparse");
	}		
	int nnzMax = mxGetNzmax(spMat_array);
    int nbCol = mxGetN(spMat_array);
    int nbRow = mxGetM(spMat_array);
    //mexPrintf("DIM (%d,%d) NNZMAX : %d\n",nbRow,nbCol,nnzMax);
    
    size_t* jc,*ir;
    double* pr;
    
    //jc = (size_t *) mxCalloc(nbCol+1,sizeof(size_t));
    jc = (size_t *)mxGetJc(spMat_array);
    //ir = (size_t *) mxCalloc(nnzMax,sizeof(size_t));
    ir = (size_t *) mxGetIr(spMat_array);
    //pr = (double *) mxCalloc(nnzMax,sizeof(double));
    pr = (double *) mxGetPr(spMat_array);


    S.set(nnzMax,nbRow,nbCol,pr,ir,jc); 
	//mxFree(jc);
	//mxFree(ir);
	//mxFree(pr);

    /*faust_mat A=S;
    mxArray*   mxA=FaustMat2mxArray(A);
    mexPrintf("INSIDE\n");
    mexCallMATLAB(0,NULL,1,&mxA,"disp");*/
}















mxArray*  FaustMat2mxArray(faust_mat M)
{		
		mxArray * mxMat;
		faust_real * mat_ptr;
		int row,col;
		row = M.getNbRow();
        col = M.getNbCol();
        mxMat = mxCreateDoubleMatrix(row,col,mxREAL);
        mat_ptr = (faust_real *) mxCalloc(row*col,sizeof(faust_real));
        memcpy(mat_ptr,M.getData(),row*col*sizeof(double));
        mxSetM(mxMat, row);
        mxSetN(mxMat, col);
		if (sizeof(double) == sizeof(faust_real))
		{
			mxSetPr(mxMat, (double *)mat_ptr);
		}else
		{
			double*	 mat_ptr_bis = (double *) mxCalloc(row*col,sizeof(double));
			for (int i=0;i<row*col;i++)
			{
				mat_ptr_bis[i] = (double) mat_ptr[i];
			}
			mxSetPr(mxMat, mat_ptr_bis);
		}		
        
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
        memcpy(mat_ptr,mat.getData(),rowFact*colFact*sizeof(double));
        mxSetM(mxMat, rowFact);
        mxSetN(mxMat, colFact);
        mxSetPr(mxMat, mat_ptr);
        mxSetCell((*cellFacts),k,mxDuplicateArray(mxMat));
    } 
    
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
		getFaustMat(mxMat,mat);
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
             //mexPrintf("NAME  %s PARAMS %d DIMS : (%d,%d)\n",consName,intParameter,nbRowCons,nbColCons);
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
void addSpmat(const mxArray * mxMat,std::vector<faust_spmat> &vec_spmat)
{
	
	faust_spmat spM;
	
	if (!mxIsSparse(mxMat))
	{
		faust_mat M;
		getFaustMat(mxMat,M);
		spM = M;
	}else
	{
		getFaustspMat(mxMat,spM);
	}
	
	int former_size = 0;
	if (vec_spmat.size() != 0)
	{
		if (vec_spmat[vec_spmat.size()-1].getNbCol() != spM.getNbRow())
		{
			mexErrMsgTxt("addSpmat : matrix dimensions mismatch");
		}	
	}
		vec_spmat.push_back(spM);	
}



void testCoherence(const mxArray* params,std::vector<bool> & presentFields)
{
  int nbr_field=mxGetNumberOfFields(params);
  presentFields.resize(8);
  presentFields.assign(8,false); 
  if(nbr_field < 3)
  {
      mexErrMsgTxt("The number of field of params must be at least 3 "); 
  }
  
  
  for (int i=0;i<nbr_field;i++)
  {     
        const char * fieldName;
        fieldName = mxGetFieldNameByNumber(params,i);
        mexPrintf("fieldname %d : %s\n",i,fieldName);
        
        if (strcmp(fieldName,"data") == 0)
        {
            presentFields[0] = true;
        }
        
        if (strcmp(fieldName,"nfacts") == 0)
        {
            presentFields[1] = true;
        }
        if (strcmp(fieldName,"cons") == 0)
        {
            presentFields[2] = true;
        }
        if (strcmp(fieldName,"niter1") == 0)
        {
            presentFields[3] = true;
        }
        if (strcmp(fieldName,"niter2") == 0)
        {
            presentFields[4] = true;
        }
        if (strcmp(fieldName,"verbose") == 0)
        {
            presentFields[5] = true;
        }
        if (strcmp(fieldName,"fact_side") == 0)
        {
            presentFields[6] = true;
        }
        if (strcmp(fieldName,"update_way") == 0)
        {
            presentFields[7] = true;
        }
        if (strcmp(fieldName,"init_lambda") == 0)
        {
            presentFields[8] = true;
        }
        if (strcmp(fieldName,"compute_lambda") == 0)
        {
            presentFields[9] = true;
        }
  }
  
}




void DisplayParams(const faust_params & params)
{   
    mexPrintf("/////// PARAMS //////\n");
    faust_mat data = params.data; 
   
    int nbRow = data.getNbRow();
    int nbCol = data.getNbCol();
    mexPrintf("DATA DIM : %d %d\n",nbRow,nbCol);
    if (nbRow*nbCol < 1000)
    {
        for (int i = 0;i<data.getNbRow();i++)
        {
            for (int j = 0;j<data.getNbCol();j++)
            {
                mexPrintf("%f ",data(i,j));
            }
            mexPrintf("\n");
        }
    }
    mexPrintf("\n\n");
    
    mexPrintf("NB_FACTS : %d\n",params.nb_fact);
    mexPrintf("CONSTRAINTS : nbr %d\n",params.cons[0].size());
    for (int L=0;L<params.cons[0].size();L++)
	{
		for (int jl=0;jl<params.cons.size();jl++)
		{	
			//char * type_cons =  getType((*params.cons[jl][L]).get_constraint_name());
            char * type_cons =  (*params.cons[jl][L]).getType();
			mexPrintf(" %s ",(*params.cons[jl][L]).get_constraint_name());
			mexPrintf(" DIMS (%d,%d) : ",(*params.cons[jl][L]).getRows(),(*params.cons[jl][L]).getCols());
			
			
			if (strcmp(type_cons,"INT") == 0)
			{	
				faust_constraint_int* const_int = (faust_constraint_int*)(params.cons[jl][L]);
				mexPrintf(" parameter : %d",(*const_int).getParameter());
			}
			
			if (strcmp(type_cons,"FAUST_REAL") == 0)
			{	
				faust_constraint_real* const_real = (faust_constraint_real*)(params.cons[jl][L]);
				mexPrintf(" parameter : %f",(*const_real).getParameter());
			}
			mexPrintf("\n"); 

			
		}
		
	}
}


