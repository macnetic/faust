#ifndef __FAUST_TOOLS_MEX_HPP__
#define __FAUST_TOOLS_MEX_HPP__


//#include "tools_mex.h"
#include "faust_ConstraintGeneric.h"
#include "faust_ConstraintFPP.h"
#include "faust_ConstraintMat.h"
#include "faust_ConstraintInt.h"
#include "faust_Params.h"



template<typename T>
void getFaustVec(const mxArray * vec_array,Faust::Vect<T> & vec)
{
	    int  nbRow,nbCol;

	if (mxIsEmpty(vec_array))
	{
		mexErrMsgTxt("tools_mex.h:getFaustVec :input matrix is empty.");
	}
	mwSize nb_dim=mxGetNumberOfDimensions(vec_array);
	if (nb_dim != 2)
	{
		mexErrMsgTxt("tools_mex.h:getFaustVec :input vector must be a 2D array.");
	}
    const mwSize *dimsMat;
    dimsMat = mxGetDimensions(vec_array);

	nbRow = (int) dimsMat[0];
    nbCol = (int) dimsMat[1];
    if ((nbRow == 0) || (nbCol == 0))
        mexErrMsgIdAndTxt("tools_mex.h:getFaustVec", "empty vector");
	if ((nbCol) > 1)
	{
		mexErrMsgTxt("getFaustVec : input must be a column-vector");
	}
    if (mxIsSparse(vec_array))
    {
        //mexErrMsgTxt("sparse matrix entry instead of dense matrix");
        mexErrMsgIdAndTxt("a","a sparse matrix entry instead of dense vector");
    }
	const mxClassID V_CLASS_ID = mxGetClassID(vec_array);
	 T* MatPtr;
	if (((V_CLASS_ID == mxDOUBLE_CLASS) && (sizeof(double) == sizeof(T))) || ((V_CLASS_ID == mxSINGLE_CLASS) && (sizeof(float) == sizeof(T))))
	{
		MatPtr = (T*) mxGetPr(vec_array);
	}else
	{
		if (V_CLASS_ID == mxDOUBLE_CLASS)
		{
			MatPtr = (T*) mxCalloc(nbRow,sizeof(T));
			double* MatPtrDouble =(double*) mxGetPr(vec_array);
			for (int i=0;i<nbRow*nbCol;i++)
				MatPtr[i] = (T) MatPtrDouble[i];
		}
		else if (V_CLASS_ID == mxSINGLE_CLASS)
		{
			MatPtr = (T*) mxCalloc(nbRow*nbCol,sizeof(T));
			float* MatPtrSingle= (float*) (mxGetData(vec_array));
			for (int i=0;i<nbRow*nbCol;i++)
				MatPtr[i] = (T) MatPtrSingle[i];


		}else
		{
		 mexErrMsgTxt("getFaustVec :input vector format must be single or double");
		}
	}

     vec.resize(nbRow);

    memcpy(vec.getData(),MatPtr,nbRow*sizeof(T));
	if (((V_CLASS_ID == mxDOUBLE_CLASS) && (sizeof(double) != sizeof(T))) || ((V_CLASS_ID == mxSINGLE_CLASS) && (sizeof(float) != sizeof(T))))
	{
		mxFree(MatPtr);
	}

}






template<typename T>
void getFaustMat(const mxArray* Mat_array,Faust::MatDense<T> & Mat)
{

    int  nbRow,nbCol;

	if (mxIsEmpty(Mat_array))
	{
		mexErrMsgTxt("tools_mex.h:getFaustMat :input matrix is empty.");
	}
	mwSize nb_dim=mxGetNumberOfDimensions(Mat_array);
	if (nb_dim != 2)
	{
		mexErrMsgTxt("tools_mex.h:getFaustMat :input matrix must be a 2D array.");
	}
    const mwSize *dimsMat;
    dimsMat = mxGetDimensions(Mat_array);

	nbRow = (int) dimsMat[0];
    nbCol = (int) dimsMat[1];
    if ((nbRow == 0) || (nbCol == 0))
        mexErrMsgIdAndTxt("tools_mex.h:getFaustMat", "empty matrix");
    if (mxIsSparse(Mat_array))
    {
        //mexErrMsgTxt("sparse matrix entry instead of dense matrix");
        mexErrMsgIdAndTxt("a","a sparse matrix entry instead of dense matrix");
    }
	const mxClassID V_CLASS_ID = mxGetClassID(Mat_array);
	 T* MatPtr;
	if (((V_CLASS_ID == mxDOUBLE_CLASS) && (sizeof(double) == sizeof(T))) || ((V_CLASS_ID == mxSINGLE_CLASS) && (sizeof(float) == sizeof(T))))
	{
		MatPtr = (T*) mxGetPr(Mat_array);
	}else
	{
		if (V_CLASS_ID == mxDOUBLE_CLASS)
		{
			MatPtr = (T*) mxCalloc(nbRow*nbCol,sizeof(T));
			double* MatPtrDouble =(double*) mxGetPr(Mat_array);
			for (int i=0;i<nbRow*nbCol;i++)
				MatPtr[i] = (T) MatPtrDouble[i];
		}
		else if (V_CLASS_ID == mxSINGLE_CLASS)
		{
			MatPtr = (T*) mxCalloc(nbRow*nbCol,sizeof(T));
			float* MatPtrSingle= (float*) (mxGetData(Mat_array));
			for (int i=0;i<nbRow*nbCol;i++)
				MatPtr[i] = (T) MatPtrSingle[i];


		}else
		{
		 mexErrMsgTxt("getFaustMat :input matrix format must be single or double");
		}
	}

     Mat.resize(nbRow,nbCol);

    memcpy(Mat.getData(),MatPtr,nbRow*nbCol*sizeof(T));
	if (((V_CLASS_ID == mxDOUBLE_CLASS) && (sizeof(double) != sizeof(T))) || ((V_CLASS_ID == mxSINGLE_CLASS) && (sizeof(float) != sizeof(T))))
	{
		mxFree(MatPtr);
	}



}

template<typename T>
void getFaustspMat(const mxArray* spMat_array,Faust::MatSparse<T> & S)
{
	if (!mxIsSparse(spMat_array))
	{
		mexErrMsgIdAndTxt("tools_mex.h:getFaustspMat",
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

    /*Faust::MatDense A=S;
    mxArray*   mxA=FaustMat2mxArray(A);
    mexPrintf("INSIDE\n");
    mexCallMATLAB(0,NULL,1,&mxA,"disp");*/
}

template<typename T>
mxArray*  FaustMat2mxArray(const Faust::MatDense<T>& M)
{
		mxArray * mxMat;
		T * mat_ptr;
		int row,col;
		row = M.getNbRow();
        col = M.getNbCol();
        mxMat = mxCreateDoubleMatrix(row,col,mxREAL);
        mat_ptr = (T *) mxCalloc(row*col,sizeof(T));
        memcpy(mat_ptr,M.getData(),row*col*sizeof(double));
        mxSetM(mxMat, row);
        mxSetN(mxMat, col);
		if (sizeof(double) == sizeof(T))
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

		//mxArray * rhs[1];
		//rhs[0]=mxMat;
        //mexCallMATLAB(0,NULL,1,rhs, "disp");
		return mxMat;
}



template<typename T>
void setCellFacts(mxArray **  cellFacts,std::vector<Faust::MatDense<T> >& facts)
{
    int rowFact,colFact;
    int nb_fact = facts.size();
    Faust::MatDense<T> mat;
    (*cellFacts) = mxCreateCellMatrix(1,nb_fact);
    mxArray * mxMat;
    T* mat_ptr;
	mwSize dims[2]={0,0};
	if (sizeof(T) == sizeof(double))
	{
		mxMat = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
	}else if (sizeof(T) == sizeof(float))
	{
		mxMat = mxCreateNumericArray(2,dims,mxSINGLE_CLASS,mxREAL);
	}else
	{
		mexErrMsgTxt("mexFaustMat : setCellFacts : T type must be equal to double or float");
	}

    for (size_t k = 0; k < nb_fact; k++)
    {
        mat = facts[k];
        rowFact = mat.getNbRow();
        colFact = mat.getNbCol();
        mxSetM(mxMat, rowFact);
        mxSetN(mxMat, colFact);

        mat_ptr = (T *) mxCalloc(rowFact*colFact,sizeof(T));

		memcpy(mat_ptr,mat.getData(),rowFact*colFact*sizeof(T));


        mxSetData(mxMat, mat_ptr);
        mxSetCell((*cellFacts),k,mxDuplicateArray(mxMat));
    }

}

template<typename T>
void setVectorFaustMat(std::vector<Faust::MatDense<T> > &vecMat,mxArray *Cells)
{
	mxArray* mxMat;
	mwSize nb_fact = mxGetNumberOfElements(Cells);
	Faust::MatDense<T> mat;
	vecMat.resize(0);
	mexPrintf("cells_size : %d\n",nb_fact);
	for (mwSize i=0;i<nb_fact;i++)
	{
		mexPrintf("i : %d\n",i);
		mxMat=mxGetCell(Cells,i);
		mexPrintf("mxMat set\n",i);
		if (mxMat == NULL)
		{
			mexErrMsgTxt("tools_mex.h:setVectorFaustMat :input matrix is empty.");
		}
		mexPrintf("mxMat test_empty\n",i);
		getFaustMat(mxMat,mat);
		mexPrintf("mat set\n",i);
		vecMat.push_back(mat);
	}
	mexPrintf("fin SetVectorFaustMat\n");
}

template<typename T>
void getConstraint(std::vector<const Faust::ConstraintGeneric<T>*> & consS,mxArray* mxCons)
{
    mwSize bufCharLen,nbRowCons,nbColCons,nb_params;
    int status;
    char * consName;
    double paramCons;
    mxArray * mxConsParams;

    nb_params = mxGetNumberOfElements(mxCons);
	if (!mxIsCell(mxCons))
		mexErrMsgTxt("tools_mex.h : getConstraint : constraint must be a cell-array. ");
	if (nb_params != 4)
		mexErrMsgTxt("tools_mex.h : getConstraint : size of constraint must be equal to 4. ");

    mxConsParams=mxGetCell(mxCons,0);
	if (!mxIsChar(mxConsParams))
		mexErrMsgTxt("tools_mex.h : getConstraint : constraint first cell of the constraint must be a character  ");
    bufCharLen = mxGetNumberOfElements(mxConsParams)+1;
    consName = (char *) mxCalloc(bufCharLen,sizeof(char));
    status = mxGetString(mxConsParams,consName,bufCharLen);
    if(status)
        mexErrMsgTxt("tools_mex.h : getConstraint : problem in mxGetString");
    paramCons = mxGetScalar(mxConsParams);
    mxConsParams=mxGetCell(mxCons,2);
    nbRowCons   = (int) mxGetScalar(mxConsParams);
    mxConsParams=mxGetCell(mxCons,3);
    nbColCons = (int) mxGetScalar(mxConsParams);

	int const_type = getTypeConstraint(consName);
	faust_constraint_name consNameType=getEquivalentConstraint(consName);

    switch(const_type)
	{
        case 0:
		{

            mxConsParams=mxGetCell(mxCons,1);
            int  intParameter = (int) (mxGetScalar(mxConsParams)+0.5);
             //mexPrintf("NAME  %s PARAMS %d DIMS : (%d,%d)\n",consName,intParameter,nbRowCons,nbColCons);
            consS.push_back(new Faust::ConstraintInt<T>(consNameType,intParameter,nbRowCons,nbColCons));
            break;
		}
		case 1:
		{
            mxConsParams=mxGetCell(mxCons,1);
			T realParameter = (T) mxGetScalar(mxConsParams);
			consS.push_back((new Faust::ConstraintFPP<T>(consNameType,realParameter,nbRowCons,nbColCons)));

			break;
		}
		case 2 :
		{
			mxConsParams=mxGetCell(mxCons,1);
			Faust::MatDense<T> matParameter;
			getFaustMat(mxConsParams,matParameter);
			consS.push_back((new Faust::ConstraintMat<T>(consNameType,matParameter,nbRowCons,nbColCons)));
			break;
		}
		default :
            mexErrMsgTxt("Unknown constraint name ");
		break;
    }




}

template<typename T>
void addSpmat(const mxArray * mxMat,std::vector<Faust::MatSparse<T> > &vec_spmat)
{

	Faust::MatSparse<T> spM;

	if (!mxIsSparse(mxMat))
	{
		Faust::MatDense<T> M;
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



// void testCoherence(const mxArray* params,std::vector<bool> & presentFields)
// {
  // int nbr_field=mxGetNumberOfFields(params);
  // presentFields.resize(8);
  // presentFields.assign(8,false);
  // if(nbr_field < 3)
  // {
      // mexErrMsgTxt("The number of field of params must be at least 3 ");
  // }


  // for (int i=0;i<nbr_field;i++)
  // {
        // const char * fieldName;
        // fieldName = mxGetFieldNameByNumber(params,i);
        // mexPrintf("fieldname %d : %s\n",i,fieldName);

        // if (strcmp(fieldName,"data") == 0)
        // {
            // presentFields[0] = true;
        // }

        // if (strcmp(fieldName,"nfacts") == 0)
        // {
            // presentFields[1] = true;
        // }
        // if (strcmp(fieldName,"cons") == 0)
        // {
            // presentFields[2] = true;
        // }
        // if (strcmp(fieldName,"niter1") == 0)
        // {
            // presentFields[3] = true;
        // }
        // if (strcmp(fieldName,"niter2") == 0)
        // {
            // presentFields[4] = true;
        // }
        // if (strcmp(fieldName,"verbose") == 0)
        // {
            // presentFields[5] = true;
        // }
        // if (strcmp(fieldName,"fact_side") == 0)
        // {
            // presentFields[6] = true;
        // }
        // if (strcmp(fieldName,"update_way") == 0)
        // {
            // presentFields[7] = true;
        // }
        // if (strcmp(fieldName,"init_lambda") == 0)
        // {
            // presentFields[8] = true;
        // }
        // if (strcmp(fieldName,"compute_lambda") == 0)
        // {
            // presentFields[9] = true;
        // }
  // }

// }



template<typename T>
void DisplayParams(const Faust::Params<T> & params)
{
    mexPrintf("/////// PARAMS //////\n");
    Faust::MatDense<T> data = params.data;

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
    for (unsigned int L=0;L<params.cons[0].size();L++)
	{
		for (unsigned int jl=0;jl<params.cons.size();jl++)
		{
			mexPrintf(" %s ", params.cons[jl][L]->get_constraint_name());
			mexPrintf(" DIMS (%d,%d) : ",(*params.cons[jl][L]).getRows(),(*params.cons[jl][L]).getCols());


			if (params.cons[jl][L]->isConstraintParameterInt())
			{
				Faust::ConstraintInt<T>* const_int = (Faust::ConstraintInt<T>*)(params.cons[jl][L]);
				mexPrintf(" parameter : %d",(*const_int).getParameter());
			}

			if (params.cons[jl][L]->isConstraintParameterReal())
			{
				Faust::ConstraintFPP<T>* const_real = (Faust::ConstraintFPP<T>*)(params.cons[jl][L]);
				mexPrintf(" parameter : %f",(*const_real).getParameter());
			}
			mexPrintf("\n");


		}

	}
}

#endif
