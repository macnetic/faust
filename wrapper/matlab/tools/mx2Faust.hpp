/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
/*                      Luc Le Magoarou, Remi Gribonval                     */
/*                      INRIA Rennes, FRANCE                                */
/*                      http://www.inria.fr/                                */
/*                                                                          */
/*  The FAuST Toolbox is distributed under the terms of the GNU Affero      */
/*  General Public License.                                                 */
/*  This program is free software: you can redistribute it and/or modify    */
/*  it under the terms of the GNU Affero General Public License as          */
/*  published by the Free Software Foundation.                              */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful, but     */
/*  WITHOUT ANY WARRANTY; without even the implied warranty of              */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    */
/*  See the GNU Affero General Public License for more details.             */
/*                                                                          */
/*  You should have received a copy of the GNU Affero General Public        */
/*  License along with this program.                                        */
/*  If not, see <http://www.gnu.org/licenses/>.                             */
/*                                                                          */
/*                             Contacts:                                    */
/*      Nicolas Bellot  : nicolas.bellot@inria.fr                           */
/*      Adrien Leman    : adrien.leman@inria.fr                             */
/*      Thomas Gautrais : thomas.gautrais@inria.fr                          */
/*      Luc Le Magoarou : luc.le-magoarou@inria.fr                          */
/*      Remi Gribonval  : remi.gribonval@inria.fr                           */
/*                                                                          */
/*                              References:                                 */
/*  [1] Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse       */
/*  approximations of matrices and applications", Journal of Selected       */
/*  Topics in Signal Processing, 2016.                                      */
/*  <https://hal.archives-ouvertes.fr/hal-01167948v1>                       */



#ifndef __FAUST_MX2FAUST_HPP__
#define __FAUST_MX2FAUST_HPP__



#include "faust_MatGeneric.h"
#include "faust_ConstraintGeneric.h"
#include "faust_ConstraintFPP.h"
#include "faust_ConstraintMat.h"
#include "faust_ConstraintInt.h"
#include "faust_Params.h"
#include "faust_MatDense.h"
#include "faust_MatSparse.h"
#include "faust_Vect.h"
#include "faust_Transform.h"
#include "faust_LinearOperator.h"



template<typename FPP>
bool isScalarCompatible(Faust::LinearOperator<FPP,Cpu> & L, const mxArray * Matlab_Mat)
{
	bool isMatlabComplex = mxIsComplex(Matlab_Mat);
	bool isTransformComplex = !(L.isReal());
	return (isMatlabComplex == isTransformComplex);

}







/*
template<typename FPP>
void mxArray2FaustMat(const mxArray* Mat_array,Faust::MatDense<FPP,Cpu> & Mat)
{

    int  nbRow,nbCol;
    



	
	if (mxIsEmpty(Mat_array))
	{
		mexErrMsgTxt("tools_mex.h:mxArray2FaustMat :input matrix is empty.");
	}
	mwSize nb_dim=mxGetNumberOfDimensions(Mat_array);
	if (nb_dim != 2)
	{
		mexErrMsgTxt("tools_mex.h:mxArray2FaustMat :input matrix must be a 2D array.");
	}
    const mwSize *dimsMat;
    dimsMat = mxGetDimensions(Mat_array);

	nbRow = (int) dimsMat[0];
    nbCol = (int) dimsMat[1];
    if ((nbRow == 0) || (nbCol == 0))
        mexErrMsgIdAndTxt("tools_mex.h:mxArray2FaustMat", "empty matrix");
    if (mxIsSparse(Mat_array))
    {
        //mexErrMsgTxt("sparse matrix entry instead of dense matrix");
        mexErrMsgIdAndTxt("a","a sparse matrix entry instead of dense matrix");
    }

    //check scalar compayibility 
    if (!isScalarCompatible(Mat,Mat_array))		
	mexErrMsgTxt("mxArray2FaustMat scalar type (complex/real) are not compatible");
	
	const mxClassID V_CLASS_ID = mxGetClassID(Mat_array);
	 FPP* MatPtr;
	if (((V_CLASS_ID == mxDOUBLE_CLASS) && (sizeof(double) == sizeof(FPP))) || ((V_CLASS_ID == mxSINGLE_CLASS) && (sizeof(float) == sizeof(FPP))))
	{
		MatPtr = (FPP*) mxGetPr(Mat_array);
	}else
	{
		if (V_CLASS_ID == mxDOUBLE_CLASS)
		{
			MatPtr = (FPP*) mxCalloc(nbRow*nbCol,sizeof(FPP));
			double* MatPtrDouble =(double*) mxGetPr(Mat_array);
			for (int i=0;i<nbRow*nbCol;i++)
				MatPtr[i] = (FPP) MatPtrDouble[i];
		}
		else if (V_CLASS_ID == mxSINGLE_CLASS)
		{
			MatPtr = (FPP*) mxCalloc(nbRow*nbCol,sizeof(FPP));
			float* MatPtrSingle= (float*) (mxGetData(Mat_array));
			for (int i=0;i<nbRow*nbCol;i++)
				MatPtr[i] = (FPP) MatPtrSingle[i];


		}else
		{
		 mexErrMsgTxt("mxArray2FaustMat :input matrix format must be single or double");
		}
	}

     Mat.resize(nbRow,nbCol);

    memcpy(Mat.getData(),MatPtr,nbRow*nbCol*sizeof(FPP));
	if (((V_CLASS_ID == mxDOUBLE_CLASS) && (sizeof(double) != sizeof(FPP))) || ((V_CLASS_ID == mxSINGLE_CLASS) && (sizeof(float) != sizeof(FPP))))
	{
		mxFree(MatPtr);
	}



}*/




template<typename FPP>
void mxArray2FaustMat(const mxArray* Mat_array,Faust::MatDense<FPP,Cpu> & Mat)
{
	int  nbRow,nbCol;
    	if (mxIsEmpty(Mat_array))
	{
		mexErrMsgTxt("tools_mex.h:mxArray2FaustMat :input matrix is empty.");
	}
	mwSize nb_dim=mxGetNumberOfDimensions(Mat_array);
	if (nb_dim != 2)
	{
		mexErrMsgTxt("tools_mex.h:mxArray2FaustMat :input matrix must be a 2D array.");
	}
    const mwSize *dimsMat;
    dimsMat = mxGetDimensions(Mat_array);

	nbRow = (int) dimsMat[0];
    nbCol = (int) dimsMat[1];
    if ((nbRow == 0) || (nbCol == 0))
        mexErrMsgIdAndTxt("tools_mex.h:mxArray2FaustMat", "empty matrix");
    if (mxIsSparse(Mat_array))
    {
        //mexErrMsgTxt("sparse matrix entry instead of dense matrix");
        mexErrMsgIdAndTxt("a","a sparse matrix entry instead of dense matrix");
    }

    Mat.resize(nbRow,nbCol);

    FPP* MatPtr;
    mxArray2Ptr(Mat_array,MatPtr);			
    memcpy(Mat.getData(),MatPtr,nbRow*nbCol*sizeof(FPP));	
    		
    if(MatPtr) {delete [] MatPtr ; MatPtr = NULL;}





}




template<typename FPP>
void mxArray2FaustspMat(const mxArray* spMat_array,Faust::MatSparse<FPP,Cpu> & S)
{
	
	if (!mxIsSparse(spMat_array))
	{
		mexErrMsgIdAndTxt("tools_mex.h:mxArray2FaustspMat",
           "input array must be sparse");
	}

    	
    int nnzMax = mxGetNzmax(spMat_array);
    int nbCol = mxGetN(spMat_array);
    int nbRow = mxGetM(spMat_array);
    //mexPrintf("DIM (%d,%d) NNZMAX : %d\n",nbRow,nbCol,nnzMax);

    size_t* jc,*ir;
    FPP* ptr_data;

    //jc = (size_t *) mxCalloc(nbCol+1,sizeof(size_t));
    jc = (size_t *)mxGetJc(spMat_array);
    //ir = (size_t *) mxCalloc(nnzMax,sizeof(size_t));
    ir = (size_t *) mxGetIr(spMat_array);
    //pr = (double *) mxCalloc(nnzMax,sizeof(double));
    
   		
    mxArray2Ptr(spMat_array,ptr_data);  
  	

    S.set(nnzMax,nbRow,nbCol,ptr_data,ir,jc);
    
    if(ptr_data) {delete [] ptr_data ; ptr_data = NULL;}
   	

}



template<typename FPP,class FUNCTOR>
void mxArray2PtrBase(const mxArray* mxMat, FPP* & ptr_data,FUNCTOR & mxGetDataFunc)
{
		if ((mxGetDataFunc != mxGetData) && (mxGetDataFunc != mxGetImagData))
			mexErrMsgTxt("mxArrayPtr(mxArray*,FPP* &, FUNCTOR) : invalid functor passed to mxArrayPtr(mxArray*,FPP* &, FUNCTOR &) must be mxGetData or mxGetImagData");
				
		

		if ((mxGetDataFunc == mxGetImagData) && (!mxIsComplex(mxMat)) )
			mexErrMsgTxt("mxArrayPtr(mxArray*,FPP* &, FUNCTOR) : can't get imaginary part of a real matrix ");

		const mxClassID V_CLASS_ID = mxGetClassID(mxMat);

		size_t nb_element_tmp;
		
		if (mxIsSparse(mxMat))
			nb_element_tmp = mxGetNzmax(mxMat);
		else
			nb_element_tmp = mxGetNumberOfElements(mxMat);		
		
		const size_t NB_ELEMENTS = nb_element_tmp;
		

		if(V_CLASS_ID == mxDOUBLE_CLASS)

		{

		    double* ptr_data_tmp = static_cast<double*> (mxGetDataFunc(mxMat));
		    ptr_data = new FPP[NB_ELEMENTS];
		    for (size_t i =0 ; i<NB_ELEMENTS ; i++)
		        ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);

		}

		else if(V_CLASS_ID == mxSINGLE_CLASS)

		{

		    float* ptr_data_tmp = static_cast<float*> (mxGetDataFunc(mxMat));
		    ptr_data = new FPP[NB_ELEMENTS];
		    for (size_t i =0 ; i<NB_ELEMENTS ; i++)
		        ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);

		}

		else if(V_CLASS_ID == mxINT8_CLASS)

		{

		    char* ptr_data_tmp = static_cast<char*> (mxGetDataFunc(mxMat));
		    ptr_data = new FPP[NB_ELEMENTS];
		    for (size_t i =0 ; i<NB_ELEMENTS ; i++)
		        ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);

		}

		else if(V_CLASS_ID == mxUINT8_CLASS)

		{

		    unsigned char* ptr_data_tmp = static_cast<unsigned char*> (mxGetDataFunc(mxMat));
		    ptr_data = new FPP[NB_ELEMENTS];
		    for (size_t i =0 ; i<NB_ELEMENTS ; i++)
		        ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);

		}

		else if(V_CLASS_ID == mxINT16_CLASS)

		{

		    short* ptr_data_tmp = static_cast<short*> (mxGetDataFunc(mxMat));
		    ptr_data = new FPP[NB_ELEMENTS];
		    for (size_t i =0 ; i<NB_ELEMENTS ; i++)
		        ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);

		}

		else if (V_CLASS_ID == mxUINT16_CLASS)

		{

		    unsigned short* ptr_data_tmp = static_cast<unsigned short*> (mxGetDataFunc(mxMat));
		    ptr_data = new FPP[NB_ELEMENTS];
		    for (size_t i =0 ; i<NB_ELEMENTS ; i++)
		        ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);
		}
		else if (V_CLASS_ID == mxINT32_CLASS)
		{
		    int* ptr_data_tmp = static_cast<int*> (mxGetDataFunc(mxMat));
		    ptr_data = new FPP[NB_ELEMENTS];
		    for (size_t i =0 ; i<NB_ELEMENTS ; i++)
		        ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);
		}
		else if (V_CLASS_ID == mxUINT32_CLASS)
		{
		    unsigned int* ptr_data_tmp = static_cast<unsigned int*> (mxGetDataFunc(mxMat));
		    ptr_data = new FPP[NB_ELEMENTS];
		    for (size_t i =0 ; i<NB_ELEMENTS ; i++)
		        ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);
		}
		else if (V_CLASS_ID == mxINT64_CLASS)
		{
		    long long* ptr_data_tmp = static_cast<long long*> (mxGetDataFunc(mxMat));
		    ptr_data = new FPP[NB_ELEMENTS];
		    for (size_t i =0 ; i<NB_ELEMENTS ; i++)
		        ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);
		}
		else if (V_CLASS_ID == mxUINT64_CLASS)
		{
		    unsigned long long* ptr_data_tmp = static_cast<unsigned long long*> (mxGetDataFunc(mxMat));
		    ptr_data = new FPP[NB_ELEMENTS];
		    for (size_t i =0 ; i<NB_ELEMENTS ; i++)
		        ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);
		}
		else
		    mexErrMsgTxt("Unknown matlab type.");

}





template<typename FPP>
void mxArray2Ptr(const mxArray* mxMat, FPP* & ptr_data)
{

	mxArray2PtrBase(mxMat,ptr_data,mxGetData);

}





template<typename FPP>
void mxArray2Ptr(const mxArray* mxMat, std::complex<FPP>* & ptr_data)
{
				
		const mxClassID V_CLASS_ID = mxGetClassID(mxMat);

		size_t nb_element_tmp;
		
		if (mxIsSparse(mxMat))
			nb_element_tmp = mxGetNzmax(mxMat);
		else
			nb_element_tmp = mxGetNumberOfElements(mxMat);		
		
		const size_t NB_ELEMENTS = nb_element_tmp;
		
		// get the real part of the Matlab Matrix
		FPP* ptr_real_part_data;
		mxArray2Ptr(mxMat,ptr_real_part_data);
		
		ptr_data = new std::complex<FPP>[NB_ELEMENTS];
		
		if (mxIsComplex(mxMat))
		{
			
			// Complex Matlab Matrix			
			// get the imaginary part of the Matlab Matrix
			FPP* ptr_imag_part_data;
			mxArray2PtrBase(mxMat,ptr_imag_part_data,mxGetImagData);
			
			// copy the values in the output vector	
			for (int i=0;i < NB_ELEMENTS;i++)
				ptr_data[i]=std::complex<FPP>(ptr_real_part_data[i],ptr_imag_part_data[i]);			
		

			if(ptr_imag_part_data) {delete [] ptr_imag_part_data ; ptr_imag_part_data = NULL;}				

			
		}else
		{
			// Real Matlab Matrix
			// copy only the real part of the matrix (the imaginary part is set to zero
			for (int i=0;i < NB_ELEMENTS;i++)
				ptr_data[i]=std::complex<FPP>(ptr_real_part_data[i],(FPP) 0.0);
	
			
		}

		if(ptr_real_part_data) {delete [] ptr_real_part_data ; ptr_real_part_data = NULL;}
}








template<typename FPP>
void concatMatGeneric(const mxArray * mxMat,std::vector<Faust::MatGeneric<FPP,Cpu> *> &list_mat)
{
	
	if (mxMat == NULL)
	   mexErrMsgTxt("concatMatGeneric : empty matlab matrix"); 

	Faust::MatGeneric<FPP,Cpu> *  M;
	
	



	if (!mxIsSparse(mxMat))
	{	
		
		Faust::MatDense<FPP,Cpu> denseM;
		mxArray2FaustMat(mxMat,denseM);
		M=denseM.Clone();
		
	}else
	{
				
		Faust::MatSparse<FPP,Cpu> spM;		
		mxArray2FaustspMat(mxMat,spM);
		M=spM.Clone();
	}

		list_mat.push_back(M);
		
	
}



template<typename FPP>
void setVectorFaustMat(std::vector<Faust::MatDense<FPP,Cpu> > &vecMat,mxArray *Cells)
{
	mxArray* mxMat;
	mwSize nbFact = mxGetNumberOfElements(Cells);
	Faust::MatDense<FPP,Cpu> mat;
	vecMat.resize(0);
	mexPrintf("cells_size : %d\n",nbFact);
	for (mwSize i=0;i<nbFact;i++)
	{
		mexPrintf("i : %d\n",i);
		mxMat=mxGetCell(Cells,i);
		mexPrintf("mxMat set\n",i);
		if (mxMat == NULL)
		{
			mexErrMsgTxt("tools_mex.h:setVectorFaustMat :input matrix is empty.");
		}
		mexPrintf("mxMat test_empty\n",i);
		mxArray2FaustMat(mxMat,mat);
		mexPrintf("mat set\n",i);
		vecMat.push_back(mat);
	}
	mexPrintf("fin SetVectorFaustMat\n");
}






template<typename FPP>
void getConstraint(std::vector<const Faust::ConstraintGeneric<FPP,Cpu>*> & consS,mxArray* mxCons)
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

	int const_type = get_type_constraint(consName);
	faust_constraint_name consNameType=get_equivalent_constraint(consName);

    switch(const_type)
	{
        case 0:
		{

            mxConsParams=mxGetCell(mxCons,1);
            int  intParameter = (int) (mxGetScalar(mxConsParams)+0.5);
             //mexPrintf("NAME  %s PARAMS %d DIMS : (%d,%d)\n",consName,intParameter,nbRowCons,nbColCons);
            consS.push_back(new Faust::ConstraintInt<FPP,Cpu>(consNameType,intParameter,nbRowCons,nbColCons));
            break;
		}
		case 1:
		{
            mxConsParams=mxGetCell(mxCons,1);
			FPP realParameter = (FPP) mxGetScalar(mxConsParams);
			consS.push_back((new Faust::ConstraintFPP<FPP,Cpu>(consNameType,realParameter,nbRowCons,nbColCons)));

			break;
		}
		case 2 :
		{
			mxConsParams=mxGetCell(mxCons,1);
			Faust::MatDense<FPP,Cpu> matParameter;
			mxArray2FaustMat(mxConsParams,matParameter);
			consS.push_back((new Faust::ConstraintMat<FPP,Cpu>(consNameType,matParameter,nbRowCons,nbColCons)));
			break;
		}
		default :
            mexErrMsgTxt("Unknown constraint name ");
		break;
    }




}





















































#endif
