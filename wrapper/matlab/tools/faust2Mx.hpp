/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.gforge.inria.fr>                         */
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



#ifndef __FAUST_FAUST2MX_HPP__
#define __FAUST_FAUST2MX_HPP__



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
mxArray*  FaustMat2mxArray(const Faust::MatDense<FPP,Cpu>& M)
{
	if (!M.isReal())
		mexErrMsgTxt("FaustMat2mxArray : Faust::MatDense must be real");
	
	mxArray * mxMat;
	int row,col;
	row = M.getNbRow();
        col = M.getNbCol();
		
        /*mxMat = mxCreateDoubleMatrix(row,col,mxREAL);
        double* mxMatdata = mxGetPr(mxMat);
	
		if (typeid(double) == typeid(FPP))
			memcpy(mxMatdata,M.getData(),sizeof(double)*row*col);			
			
		else
		{
			double*	 mat_ptr = (double *) mxCalloc(row*col,sizeof(double));
			for (int i=0;i<row*col;i++)
			{
				mat_ptr[i] = (double) M.getData()[i];
			}
			memcpy(mxMatdata,mat_ptr,sizeof(double)*row*col);
			mxFree(mat_ptr);
			
		}*/
	const mwSize dims[2]={row,col};
	if(typeid(FPP)==typeid(float))
	{
		mxMat = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
	}else if(sizeof(FPP)==sizeof(double))
	{
		mxMat = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);		
	}else
	{
		mexErrMsgTxt("FaustMat2mxArray : unsupported type of float");
	}
	
	FPP*    ptr_out = static_cast<FPP*> (mxGetData(mxMat));
		memcpy(ptr_out, M.getData(),row*col*sizeof(FPP));	
		
	
		return mxMat;
		
}



template<typename FPP>
mxArray*  FaustMat2mxArray(const Faust::MatDense<std::complex<FPP>,Cpu>& M)
{
	mxArray * mxMat;
	int row,col;
	row = M.getNbRow();
        col = M.getNbCol();
		
  
	const mwSize dims[2]={row,col};
	if(typeid(FPP)==typeid(float))
	{
		mxMat = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxCOMPLEX);
	}else if(sizeof(FPP)==sizeof(double))
	{
		mxMat = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);		
	}else
	{
		mexErrMsgTxt("FaustMat2mxArray (complex) : unsupported type of float");
	}
	
	FPP*    ptr_real_data = static_cast<FPP*> (mxGetData(mxMat));
	FPP*    ptr_imag_data = static_cast<FPP*> (mxGetImagData(mxMat));
			
	splitComplexPtr(M.getData(),row*col,ptr_real_data,ptr_imag_data);	
	
	return mxMat;


}


template<typename FPP>
mxArray*  FaustVec2mxArray(const Faust::Vect<std::complex<FPP>,Cpu>& v)
{
	mxArray * mxMat;
	int row,col;
	row = v.size();
        col = 1;
		
  
	const mwSize dims[2]={row,col};
	if(typeid(FPP)==typeid(float))
	{
		mxMat = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxCOMPLEX);
	}else if(sizeof(FPP)==sizeof(double))
	{
		mxMat = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);		
	}else
	{
		mexErrMsgTxt("FaustMat2mxArray (complex) : unsupported type of float");
	}
	
	FPP*    ptr_real_data = static_cast<FPP*> (mxGetData(mxMat));
	FPP*    ptr_imag_data = static_cast<FPP*> (mxGetImagData(mxMat));
			
	splitComplexPtr(v.getData(),row*col,ptr_real_data,ptr_imag_data);	
	
	return mxMat;


}





template<typename FPP>
void splitComplexPtr(const std::complex<FPP>*  cpx_ptr, int nb_element, FPP* & real_ptr, FPP* & imag_ptr)
{
	std::complex<FPP> cpxValue;	
	for (int i=0;i<nb_element;i++)
	{
		cpxValue=cpx_ptr[i];
		real_ptr[i] = cpxValue.real();
		imag_ptr[i] = cpxValue.imag();
	}
	


}



template<typename FPP>
mxArray*  FaustVec2mxArray(const Faust::Vect<FPP,Cpu>& M)
{
	if (!M.isReal())
		mexErrMsgTxt("FaustMat2mxArray : Faust::MatDense must be real");
	
	mxArray * mxMat;
	int row,col;
	row = M.size();
        col = 1;
	
		
	

	const mwSize dims[2]={row,col};
	if(typeid(FPP)==typeid(float))
	{
		mxMat = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
	}else if(sizeof(FPP)==sizeof(double))
	{
		mxMat = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);		
	}else
	{
		mexErrMsgTxt("FaustVec2mxArray : unsupported type of float");
	}
	
	FPP*    ptr_out = static_cast<FPP*> (mxGetData(mxMat));
		memcpy(ptr_out, M.getData(),row*col*sizeof(FPP));	
		return mxMat;
		
}




template<typename FPP>
void setCellFacts(mxArray **  cellFacts,std::vector<Faust::MatDense<FPP,Cpu> >& facts)
{
    int rowFact,colFact;
    int nbFact = facts.size();
    Faust::MatDense<FPP,Cpu> mat;
    (*cellFacts) = mxCreateCellMatrix(1,nbFact);
    mxArray * mxMat;
    FPP* mat_ptr;
	mwSize dims[2]={0,0};
	if (sizeof(FPP) == sizeof(double))
	{
		mxMat = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
	}else if (sizeof(FPP) == sizeof(float))
	{
		mxMat = mxCreateNumericArray(2,dims,mxSINGLE_CLASS,mxREAL);
	}else
	{
		mexErrMsgTxt("mexFaustMat : setCellFacts : FPP type must be equal to double or float");
	}

    for (size_t k = 0; k < nbFact; k++)
    {
        mat = facts[k];
        rowFact = mat.getNbRow();
        colFact = mat.getNbCol();
        mxSetM(mxMat, rowFact);
        mxSetN(mxMat, colFact);

        mat_ptr = (FPP *) mxCalloc(rowFact*colFact,sizeof(FPP));

		memcpy(mat_ptr,mat.getData(),rowFact*colFact*sizeof(FPP));


        mxSetData(mxMat, mat_ptr);
        mxSetCell((*cellFacts),k,mxDuplicateArray(mxMat));
    }

}








#endif
