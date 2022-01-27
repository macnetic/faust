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
#include "faust_ParamsFGFT.h"
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

template<typename FPP, FDevice DEV>
void mxArray2FaustMat(const mxArray* Mat_array,Faust::MatDense<FPP,DEV> & Mat)
{
	int  nbRow,nbCol;
	if (mxIsEmpty(Mat_array))
	{
		mexErrMsgTxt("tools_mex.h:mxArray2FaustMat :input matrix is empty.");
	}
	mwSize nb_dim = mxGetNumberOfDimensions(Mat_array);
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
//	memcpy(Mat.getData(),MatPtr,nbRow*nbCol*sizeof(FPP));

	Mat.setData(MatPtr, nbRow, nbCol);
	if(MatPtr) {delete [] MatPtr ; MatPtr = NULL;}
}

template<typename FPP, FDevice DEV>
Faust::MatDense<FPP,DEV>* mxArray2FaustMat(const mxArray* Mat_array)
{
	Faust::MatDense<FPP,DEV> *M = new Faust::MatDense<FPP,DEV>();
	mxArray2FaustMat(Mat_array, *M);
	return M;
}

template<typename FPP, FDevice DEV>
void mxArray2FaustspMat(const mxArray* spMat_array,Faust::MatSparse<FPP,DEV> & S)
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


	mxArray2Ptr(spMat_array,ptr_data); // NOTE: it's ok that FPP == float while spMat_array is a double (the conversion is handled mxArray2PtrBase)


	S.set(nnzMax,nbRow,nbCol,ptr_data,ir,jc);

	if(ptr_data) {delete [] ptr_data ; ptr_data = NULL;}


}

#ifdef MX_HAS_INTERLEAVED_COMPLEX
template<typename FPP>
void mxArray2PtrBase(const mxArray* mxMat, FPP* & ptr_data)
{


	const mxClassID V_CLASS_ID = mxGetClassID(mxMat);

	size_t nb_element_tmp;

	if (mxIsSparse(mxMat))
		nb_element_tmp = mxGetNzmax(mxMat);
	else
		nb_element_tmp = mxGetNumberOfElements(mxMat);

	const size_t nb_elts = nb_element_tmp;


	if(V_CLASS_ID == mxDOUBLE_CLASS)

	{
		if(mxIsComplex(mxMat))
		{
			if(! is_same<FPP, complex<double>>::value)
				mexErrMsgTxt("mxMat is complex double, the output buffer must be complex<double>");
			ptr_data = new FPP[nb_elts];
			mxComplexDouble* mx_cplx_doubles = mxGetComplexDoubles(mxMat);
			memcpy(ptr_data, mx_cplx_doubles, sizeof(FPP)*nb_elts);
		}
		else
		{
			if(! is_same<FPP, double>::value)
			{
				mexErrMsgTxt("mxMat is double, the output buffer must be double");
			}
			ptr_data = new FPP[nb_elts];
			mxDouble* mx_doubles = mxGetDoubles(mxMat);
			memcpy(ptr_data, mx_doubles, sizeof(FPP)*nb_elts);
		}

	}

	else if(V_CLASS_ID == mxSINGLE_CLASS)

	{

		if(mxIsComplex(mxMat))
		{
			if(! is_same<FPP, complex<float>>::value)
				mexErrMsgTxt("mxMat is complex float, the output buffer must be complex<float>");
			ptr_data = new FPP[nb_elts];
			mxComplexSingle* mx_cplx_floats = mxGetComplexSingles(mxMat);
			memcpy(ptr_data, mx_cplx_floats, sizeof(FPP)*nb_elts);
		}
		else
		{
			if(! is_same<FPP, float>::value)
				mexErrMsgTxt("mxMat is float, the output buffer must be float");
			ptr_data = new FPP[nb_elts];
			mxSingle* mx_floats = mxGetSingles(mxMat);
			memcpy(ptr_data, mx_floats, sizeof(FPP)*nb_elts);
		}
	}

	else if(mxIsComplex(mxMat))
	{
		mexErrMsgTxt("Don't handle complex mxArray if the class is not double or single.");
	}

	else if(V_CLASS_ID == mxINT8_CLASS)

	{


		if(! is_same<FPP, char>::value)
			mexErrMsgTxt("mxMat is int8s, the output buffer must be char");
		ptr_data = new FPP[nb_elts];
		mxInt8* mx_chars = mxGetInt8s(mxMat);
		memcpy(ptr_data, mx_chars, sizeof(FPP)*nb_elts);

	}

	else if(V_CLASS_ID == mxUINT8_CLASS)

	{
		if(! is_same<FPP, unsigned char>::value)
			mexErrMsgTxt("mxMat is uint8, the output buffer must be unsigned char");
		ptr_data = new FPP[nb_elts];
		mxUint8* mx_chars = mxGetUint8s(mxMat);
		memcpy(ptr_data, mx_chars, sizeof(FPP)*nb_elts);
	}

	else if(V_CLASS_ID == mxINT16_CLASS)

	{
		if(! is_same<FPP, short>::value)
			mexErrMsgTxt("mxMat is int16s, the output buffer must be short");
		ptr_data = new FPP[nb_elts];
		mxInt16* mx_chars = mxGetInt16s(mxMat);
		memcpy(ptr_data, mx_chars, sizeof(FPP)*nb_elts);
	}

	else if (V_CLASS_ID == mxUINT16_CLASS)

	{
		if(! is_same<FPP, unsigned short>::value)
			mexErrMsgTxt("mxMat is uint16, the output buffer must be unsigned short");
		ptr_data = new FPP[nb_elts];
		mxUint16* mx_ints = mxGetUint16s(mxMat);
		memcpy(ptr_data, mx_ints, sizeof(FPP)*nb_elts);
	}
	else if (V_CLASS_ID == mxINT32_CLASS)
	{
		if(! is_same<FPP, int>::value)
			mexErrMsgTxt("mxMat is int32s, the output buffer must be int");
		ptr_data = new FPP[nb_elts];
		mxInt32* mx_ints = mxGetInt32s(mxMat);
		memcpy(ptr_data, mx_ints, sizeof(FPP)*nb_elts);
	}
	else if (V_CLASS_ID == mxUINT32_CLASS)
	{
		if(! is_same<FPP, unsigned int>::value)
			mexErrMsgTxt("mxMat is uint32, the output buffer must be unsigned int");
		ptr_data = new FPP[nb_elts];
		mxUint32* mx_ints = mxGetUint32s(mxMat);
		memcpy(ptr_data, mx_ints, sizeof(FPP)*nb_elts);
	}
	else if (V_CLASS_ID == mxINT64_CLASS)
	{
		if(! is_same<FPP, long int>::value)
			mexErrMsgTxt("mxMat is int64s, the output buffer must be long int");
		ptr_data = new FPP[nb_elts];
		mxInt64* mx_ints = mxGetInt64s(mxMat);
		memcpy(ptr_data, mx_ints, sizeof(FPP)*nb_elts);
	}
	else if (V_CLASS_ID == mxUINT64_CLASS)
	{
		if(! is_same<FPP, unsigned long int>::value)
			mexErrMsgTxt("mxMat is uint64, the output buffer must be unsigned long int");
		ptr_data = new FPP[nb_elts];
		mxUint64* mx_ints = mxGetUint64s(mxMat);
		memcpy(ptr_data, mx_ints, sizeof(FPP)*nb_elts);
	}
	else
		mexErrMsgTxt("Unknown matlab type.");

}

#else // MX_HAS_INTERLEAVED_COMPLEX not def
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

	const size_t nb_elts = nb_element_tmp;


	if(V_CLASS_ID == mxDOUBLE_CLASS)

	{
		double* ptr_data_tmp = static_cast<double*> (mxGetDataFunc(mxMat));
		ptr_data = new FPP[nb_elts];
		for (size_t i =0 ; i<nb_elts ; i++)
			ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);

	}

	else if(V_CLASS_ID == mxSINGLE_CLASS)

	{


		float* ptr_data_tmp = static_cast<float*> (mxGetDataFunc(mxMat));
		ptr_data = new FPP[nb_elts];
		for (size_t i =0 ; i<nb_elts ; i++)
			ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);

	}


	else if(V_CLASS_ID == mxINT8_CLASS)

	{


		char* ptr_data_tmp = static_cast<char*> (mxGetDataFunc(mxMat));
		ptr_data = new FPP[nb_elts];
		for (size_t i =0 ; i<nb_elts ; i++)
			ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);

	}

	else if(V_CLASS_ID == mxUINT8_CLASS)

	{
		unsigned char* ptr_data_tmp = static_cast<unsigned char*> (mxGetDataFunc(mxMat));
		ptr_data = new FPP[nb_elts];
		for (size_t i =0 ; i<nb_elts ; i++)
			ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);

	}

	else if(V_CLASS_ID == mxINT16_CLASS)

	{


		short* ptr_data_tmp = static_cast<short*> (mxGetDataFunc(mxMat));
		ptr_data = new FPP[nb_elts];
		for (size_t i =0 ; i<nb_elts ; i++)
			ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);
	}

	else if (V_CLASS_ID == mxUINT16_CLASS)

	{
		unsigned short* ptr_data_tmp = static_cast<unsigned short*> (mxGetDataFunc(mxMat));
		ptr_data = new FPP[nb_elts];
		for (size_t i =0 ; i<nb_elts ; i++)
			ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);
	}
	else if (V_CLASS_ID == mxINT32_CLASS)
	{

		int* ptr_data_tmp = static_cast<int*> (mxGetDataFunc(mxMat));
		ptr_data = new FPP[nb_elts];
		for (size_t i =0 ; i<nb_elts ; i++)
			ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);
	}
	else if (V_CLASS_ID == mxUINT32_CLASS)
	{

		unsigned int* ptr_data_tmp = static_cast<unsigned int*> (mxGetDataFunc(mxMat));
		ptr_data = new FPP[nb_elts];
		for (size_t i =0 ; i<nb_elts ; i++)
			ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);
	}
	else if (V_CLASS_ID == mxINT64_CLASS)
	{

		long long* ptr_data_tmp = static_cast<long long*> (mxGetDataFunc(mxMat));
		ptr_data = new FPP[nb_elts];
		for (size_t i =0 ; i<nb_elts ; i++)
			ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);
	}
	else if (V_CLASS_ID == mxUINT64_CLASS)
	{

		unsigned long long* ptr_data_tmp = static_cast<unsigned long long*> (mxGetDataFunc(mxMat));
		ptr_data = new FPP[nb_elts];
		for (size_t i =0 ; i<nb_elts ; i++)
			ptr_data[i] = static_cast<FPP> (ptr_data_tmp[i]);
	}
	else
		mexErrMsgTxt("Unknown matlab type.");

}
#endif

template<typename FPP>
void mxArray2Ptr(const mxArray* mxMat, FPP* & ptr_data)
{

#ifdef MX_HAS_INTERLEAVED_COMPLEX
	mxArray2PtrBase(mxMat,ptr_data);
#else
	mxArray2PtrBase(mxMat,ptr_data,mxGetData);
#endif

}




#ifdef MX_HAS_INTERLEAVED_COMPLEX
template<typename FPP>
void mxArray2Ptr(const mxArray* mxMat, std::complex<FPP>* & ptr_data)
{
	mxArray2PtrBase(mxMat, ptr_data);
}
#else
template<typename FPP>
void mxArray2Ptr(const mxArray* mxMat, std::complex<FPP>* & ptr_data)
{

	const mxClassID V_CLASS_ID = mxGetClassID(mxMat);

	size_t nb_element_tmp;

	if (mxIsSparse(mxMat))
		nb_element_tmp = mxGetNzmax(mxMat);
	else
		nb_element_tmp = mxGetNumberOfElements(mxMat);

	const size_t nb_elts = nb_element_tmp;

	// get the real part of the Matlab Matrix
	FPP* ptr_real_part_data;
	mxArray2Ptr(mxMat,ptr_real_part_data);

	ptr_data = new std::complex<FPP>[nb_elts];

	if (mxIsComplex(mxMat))
	{

		// Complex Matlab Matrix
		// get the imaginary part of the Matlab Matrix
		FPP* ptr_imag_part_data;
		mxArray2PtrBase(mxMat,ptr_imag_part_data,mxGetImagData);

		// copy the values in the output vector
		for (int i=0;i < nb_elts;i++)
			ptr_data[i]=std::complex<FPP>(ptr_real_part_data[i],ptr_imag_part_data[i]);


		if(ptr_imag_part_data) {delete [] ptr_imag_part_data ; ptr_imag_part_data = NULL;}


	}else
	{
		// Real Matlab Matrix
		// copy only the real part of the matrix (the imaginary part is set to zero
		for (int i=0;i < nb_elts;i++)
			ptr_data[i]=std::complex<FPP>(ptr_real_part_data[i],(FPP) 0.0);


	}

	if(ptr_real_part_data) {delete [] ptr_real_part_data ; ptr_real_part_data = NULL;}
}
#endif








template<typename FPP, FDevice DEV>
void concatMatGeneric(const mxArray * mxMat,std::vector<Faust::MatGeneric<FPP,DEV> *> &list_mat)
{

	if (mxMat == NULL)
		mexErrMsgTxt("concatMatGeneric : empty matlab matrix");

	Faust::MatGeneric<FPP,DEV> *  M = nullptr;
	Faust::MatDense<FPP,DEV> *  dM = nullptr;
	Faust::MatSparse<FPP,DEV> *  sM = nullptr;

	if (!mxIsSparse(mxMat))
		M = dM = mxArray2FaustMat<FPP, DEV>(mxMat);
	else
	{
		M = sM = new Faust::MatSparse<FPP,DEV>();
		mxArray2FaustspMat(mxMat, *sM);
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
//	mexPrintf("cells_size : %d\n",nbFact);
	for (mwSize i=0;i<nbFact;i++)
	{
//		mexPrintf("i : %d\n",i);
		mxMat=mxGetCell(Cells,i);
//		mexPrintf("mxMat set\n",i);
		if (mxMat == NULL)
		{
			mexErrMsgTxt("tools_mex.h:setVectorFaustMat :input matrix is empty.");
		}
//		mexPrintf("mxMat test_empty\n",i);
		mxArray2FaustMat(mxMat,mat);
//		mexPrintf("mat set\n",i);
		vecMat.push_back(mat);
	}
//	mexPrintf("fin SetVectorFaustMat\n");
}






template<typename FPP, typename FPP2>
void getConstraint(std::vector<const Faust::ConstraintGeneric*> & consS,mxArray* mxCons)
{
	mwSize bufCharLen,nbRowCons,nbColCons,nb_params, param_sz;
	int status;
	char * consName;
	double paramCons;
	mxArray * mxConsParams;
	if (!mxIsCell(mxCons))
		mexErrMsgTxt("tools_mex.h : getConstraint : constraint must be a cell-array. ");
	nb_params = mxGetNumberOfElements(mxCons);
	if (nb_params < 4 || nb_params > 5)
		mexErrMsgTxt("mx2Faust.hpp: getConstraint : size of constraint must be 4 or 5. ");

//	mexPrintf("getConstraint() nb_params=%d\n", nb_params);
	mxConsParams=mxGetCell(mxCons,0);
	if (!mxIsChar(mxConsParams))
		mexErrMsgTxt("tools_mex.h : getConstraint : constraint first cell of the constraint must be a character  ");
	//    bufCharLen = mxGetNumberOfElements(mxConsParams)+1;
	bufCharLen = mxGetN(mxConsParams)+1; //equivalent to mxGetNumberOfElements() (chars are organized in col)
	// but mxGetString() doc. indicates to use mxGetN/M()
	//    mexPrintf("getConstraint() bufCharLen= %d\n", bufCharLen);
	consName = (char *) mxCalloc(bufCharLen,sizeof(char));
	status = mxGetString(mxConsParams,consName,bufCharLen);
	if(status)
		mexErrMsgTxt("tools_mex.h : getConstraint : problem in mxGetString");
	//	mexPrintf("getConstraint() consName = %s\n", consName);
	paramCons = mxGetScalar(mxConsParams);
	mxConsParams=mxGetCell(mxCons,2);
	nbRowCons   = (int) mxGetScalar(mxConsParams);
	mxConsParams=mxGetCell(mxCons,3);
	nbColCons = (int) mxGetScalar(mxConsParams);

	int const_type = get_type_constraint(consName);
	faust_constraint_name consNameType=get_equivalent_constraint(consName);
	if(const_type != 2 && nb_params != 4)
		mexErrMsgTxt("mx2Faust.hpp: getConstraint for this constraint type (non-matrix) must be 4.");
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
				FPP2 realParameter = (FPP2) mxGetScalar(mxConsParams);
				consS.push_back((new Faust::ConstraintFPP<FPP,Cpu, FPP2>(consNameType,realParameter,nbRowCons,nbColCons)));

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

	mxFree(consName); // deallocate even if the buf is in matlab managed memory
	// the documentation advise to free managed buf in an opt. goal


}


template<typename SCALAR, typename FPP2>
const Params<SCALAR, Cpu, FPP2>* mxArray2FaustParams(const mxArray* matlab_params)
{

	//TODO: this function should be refactored into several small functions (one per field to parse from matlab)

	std::vector<bool> presentFields;
	// test if all the needed parameters are present
	testCoherence(matlab_params,presentFields);


	mxArray    *mxCurrentField,*mxCurrentCons;



	// nrow of the matrix that will be factorized
	int nb_row = -1;
	if (presentFields[NROW])
	{
		mxCurrentField = mxGetField(matlab_params,0,"nrow");
		nb_row =(int)  mxGetScalar(mxCurrentField);
	}else
	{
		mexErrMsgTxt("params.nrow must be specified");
	}

	//nb column of the matrix that will be factorized
	int nb_col = -1;
	if (presentFields[NCOL])
	{
		mxCurrentField = mxGetField(matlab_params,0,"ncol");
		nb_col =(int)  mxGetScalar(mxCurrentField);
	}else
	{
		mexErrMsgTxt("params.nrow must be specified");
	}


	//nbFact initialization
	int nbFact = -1;
	if (presentFields[NFACTS])
	{

		mxCurrentField = mxGetField(matlab_params,0,"nfacts");
		nbFact =(int)  mxGetScalar(mxCurrentField);

	}else
	{
		mexErrMsgTxt("params.nfacts must be specified");
	}


	//constraints
	std::vector<std::vector<const Faust::ConstraintGeneric*> > consSS;
	if (presentFields[CONS])
	{
		mwSize nbRowCons,nbColCons;
		mxCurrentField = mxGetField(matlab_params,0,"cons");

		if(!mxIsCell(mxCurrentField))
		{
			mexErrMsgTxt("cons must be a cell-array");
		}
		nbRowCons = mxGetM(mxCurrentField);
		nbColCons = mxGetN(mxCurrentField);

		/*if(nbRowCons !=2)
		  {
		  mexPrintf("\n cons has %d rows \n",nbRowCons);
		  mexErrMsgTxt("cons must have 2 rows");
		  }*/
		/*if(nbColCons != (nbFact-1))
		  {
		  mexPrintf("\n cons has %d cols and nbFact = %d\n",nbColCons,nbFact);
		  mexErrMsgTxt("incoherence between the number of columns of cons and nfacts ");
		  }*/
		//mexPrintf("\n cons has %d rows and %d cols \n",nbRowCons,nbColCons);
		//Faust::ConstraintGeneric * consToAdd;
		std::vector<const Faust::ConstraintGeneric*> consS;

		for (mwSize i=0;i<nbRowCons;i++)
		{
			//            mexPrintf("%d / %d\n",i,nbRowCons);
			for (mwSize j=0;j<nbColCons;j++)
			{
				//mexPrintf("cons(%d , %d)\n",i,j);
				mxCurrentCons=mxGetCell(mxCurrentField,i+(j*nbRowCons));
				getConstraint<SCALAR,FPP2>(consS,mxCurrentCons);
				//consS.push_back(consToAdd);
			}
			consSS.push_back(consS);
			consS.resize(0);
		}

	}else
	{
		mexErrMsgTxt("params.cons must be specified");
	}

	//TODO: replace by default values as constants from StoppingCriterion class
	bool is_criterion_error = false;
	int num_its = 500;
	FPP2 error_treshold = 0.3;
	int max_num_its = 10000;
	if(presentFields[NITER1])
	{
		mxCurrentField = mxGetField(matlab_params, 0,  mat_field_type2str(NITER1).c_str());
		num_its = (int) mxGetScalar(mxCurrentField);
	}
	if(presentFields[SC_IS_CRITERION_ERROR]){
		mxCurrentField = mxGetField(matlab_params, 0, mat_field_type2str(SC_IS_CRITERION_ERROR).c_str());
		is_criterion_error =  (bool) mxGetScalar(mxCurrentField);
	}
	if(presentFields[SC_ERROR_TRESHOLD])
	{
		mxCurrentField = mxGetField(matlab_params, 0, mat_field_type2str(SC_ERROR_TRESHOLD).c_str());
		error_treshold = (FPP2) mxGetScalar(mxCurrentField);
	}
	if(presentFields[SC_MAX_NUM_ITS])
	{
		mxCurrentField = mxGetField(matlab_params, 0, mat_field_type2str(SC_MAX_NUM_ITS).c_str());
		max_num_its = (int) mxGetScalar(mxCurrentField);
	}
	Faust::StoppingCriterion<FPP2> crit1(num_its, is_criterion_error, error_treshold, max_num_its);
	//TODO: replace by default values as constants from StoppingCriterion class
	is_criterion_error = false;
	num_its = 500;
	error_treshold = 0.3;
	max_num_its = 10000;
	if(presentFields[NITER2])
	{
		mxCurrentField = mxGetField(matlab_params, 0, mat_field_type2str(NITER2).c_str());
		num_its = (int) mxGetScalar(mxCurrentField);
	}
	if(presentFields[SC_IS_CRITERION_ERROR2]){
		mxCurrentField = mxGetField(matlab_params, 0, mat_field_type2str(SC_IS_CRITERION_ERROR2).c_str());
		is_criterion_error =  (bool) mxGetScalar(mxCurrentField);
	}
	if(presentFields[SC_ERROR_TRESHOLD2])
	{
		mxCurrentField = mxGetField(matlab_params, 0, mat_field_type2str(SC_ERROR_TRESHOLD2).c_str());
		error_treshold = (FPP2) mxGetScalar(mxCurrentField);
	}
	if(presentFields[SC_MAX_NUM_ITS2])
	{
		mxCurrentField = mxGetField(matlab_params, 0, mat_field_type2str(SC_MAX_NUM_ITS2).c_str());
		max_num_its = (int) mxGetScalar(mxCurrentField);
	}
	Faust::StoppingCriterion<FPP2> crit2(num_its, is_criterion_error, error_treshold, max_num_its);
	//init_facts
	std::vector<Faust::MatDense<SCALAR,Cpu> > init_facts;
	if (presentFields[INIT_FACTS])
	{
		mxCurrentField = mxGetField(matlab_params,0,mat_field_type2str(INIT_FACTS).c_str());
		setVectorFaustMat(init_facts,mxCurrentField);

	}
	//verbosity
	bool isVerbose = false;
	if (presentFields[VERBOSE])
	{
		mxCurrentField = mxGetField(matlab_params,0,mat_field_type2str(VERBOSE).c_str());
		isVerbose =(bool)  mxGetScalar(mxCurrentField);
	}

	//fact_side
	bool factside = false;
	if (presentFields[FACT_SIDE])
	{
		mxCurrentField = mxGetField(matlab_params,0,mat_field_type2str(FACT_SIDE).c_str());
		factside =(bool)  mxGetScalar(mxCurrentField);
	}

	//update_way
	bool updateway = false;
	if (presentFields[UPDATE_WAY])
	{
		mxCurrentField = mxGetField(matlab_params,0,mat_field_type2str(UPDATE_WAY).c_str());
		updateway =(bool)  mxGetScalar(mxCurrentField);
	}


	//init_lambda
	FPP2 init_lambda = 1.0;
	if (presentFields[INIT_LAMBDA])
	{
		mxCurrentField = mxGetField(matlab_params,0,mat_field_type2str(INIT_LAMBDA).c_str());
//		SCALAR* tmp_ptr = &init_lambda;
		// it works whatever mxCurrentField class is (complex or not)
//		mxArray2Ptr<SCALAR>(const_cast<const mxArray*>(mxCurrentField), tmp_ptr);
		//       init_lambda = (SCALAR) mxGetScalar(mxCurrentField);
		init_lambda = (FPP2) mxGetScalar(mxCurrentField);
	}

	Faust::Params<SCALAR,Cpu,FPP2>* params;
	if(presentFields[INIT_D])
	{
		//get the diagonal vector to define the init_D matrix (cf. FactHierarchicalF(G)FT
		SCALAR* init_D = new SCALAR[nb_row]; //nb_col == nb_row when using FactHierarchicalF(G)FT
		mxCurrentField = mxGetField(matlab_params,0,mat_field_type2str(INIT_D).c_str());
		mxArray2Ptr<SCALAR>(const_cast<const mxArray*>(mxCurrentField), init_D);
		params = new ParamsFGFT<SCALAR,Cpu,FPP2>(nb_row,nb_col,nbFact,consSS, init_facts, init_D, crit1,crit2,isVerbose,updateway,factside,init_lambda);
		delete [] init_D;
	}
	else
	{
		params = new Params<SCALAR,Cpu,FPP2>(nb_row,nb_col,nbFact,consSS,/*std::vector<Faust::MatDense<SCALAR,Cpu> >()*/ init_facts,crit1,crit2,isVerbose,updateway,factside,init_lambda);
	}


	FactorsFormat factors_format = Params<SCALAR, Cpu, FPP2>::defaultFactorsFormat;
	if(presentFields[FACTOR_FORMAT])
	{
		mxCurrentField = mxGetField(matlab_params, 0, mat_field_type2str(FACTOR_FORMAT).c_str());
		factors_format = static_cast<FactorsFormat>((int)mxGetScalar(mxCurrentField));
	}

	bool packing_RL = Params<SCALAR, Cpu, FPP2>::defaultPackingRL;
	if(presentFields[PACKING_RL])
	{
		mxCurrentField = mxGetField(matlab_params, 0, mat_field_type2str(PACKING_RL).c_str());
		packing_RL = (bool) mxGetScalar(mxCurrentField);
	}
	FPP2 norm2_threshold =	FAUST_PRECISION;
	int norm2_max_iter = FAUST_NORM2_MAX_ITER;
	if(presentFields[NORM2_THRESHOLD])
	{
		mxCurrentField = mxGetField(matlab_params, 0, mat_field_type2str(NORM2_THRESHOLD).c_str());
		norm2_threshold = (FPP2) mxGetScalar(mxCurrentField);
	}
	if(presentFields[NORM2_MAX_ITER])
	{
		mxCurrentField = mxGetField(matlab_params, 0, mat_field_type2str(NORM2_MAX_ITER).c_str());
		norm2_max_iter = (int) mxGetScalar(mxCurrentField);
	}
	if(norm2_max_iter)
		params->norm2_max_iter = norm2_max_iter;
	if(norm2_threshold != FPP2(0))
		params->norm2_threshold = norm2_threshold;
	params->packing_RL = packing_RL;
	params->factors_format = factors_format;
	bool constant_step_size = false;
	if(presentFields[CONSTANT_STEP_SIZE])
	{
		mxCurrentField = mxGetField(matlab_params, 0, mat_field_type2str(CONSTANT_STEP_SIZE).c_str());
		constant_step_size = (bool) mxGetScalar(mxCurrentField);
	}
	FPP2 step_size = FAUST_PRECISION;
	if(presentFields[STEP_SIZE])
	{
		mxCurrentField = mxGetField(matlab_params, 0, mat_field_type2str(STEP_SIZE).c_str());
		step_size = (FPP2) mxGetScalar(mxCurrentField);

	}
	params->step_size = step_size;
	params->isConstantStepSize = constant_step_size;
	bool no_normalization = Params<SCALAR, Cpu, FPP2>::defaultNoNormalization;
	if(presentFields[NO_NORMALIZATION])
	{
		mxCurrentField = mxGetField(matlab_params, 0, mat_field_type2str(NO_NORMALIZATION).c_str());
		no_normalization = (bool) mxGetScalar(mxCurrentField);
	}
	params->no_normalization = no_normalization;
	bool no_lambda = Params<SCALAR, Cpu, FPP2>::defaultNoLambda;
	if(presentFields[NO_LAMBDA])
	{
		mxCurrentField = mxGetField(matlab_params, 0, mat_field_type2str(NO_LAMBDA).c_str());
		no_lambda = (bool) mxGetScalar(mxCurrentField);
	}
	params->no_lambda = no_lambda;
	return params;
}

template<typename SCALAR, typename FPP2>
const ParamsPalm<SCALAR,Cpu,FPP2>* mxArray2FaustParamsPALM4MSA(const mxArray* matlab_params, std::vector<bool> & presentFields)
{
	Faust::ParamsPalm<SCALAR, Cpu, FPP2>* params = nullptr;
	mxArray    *mxCurrentField,*mxCurrentCons;


	// data initialisation
	Faust::MatDense<SCALAR,Cpu> data;
	if (presentFields[0])
	{
		mxCurrentField = mxGetField(matlab_params,0,"data");

		mxArray2FaustMat(  mxCurrentField,data ) ;
		/*mexPrintf("DATA");
		  for (int i = 0;i<data.getNbRow();i++)
		  {
		  for (int j = 0;j<data.getNbCol();j++)
		  {
		//bidon = std::snprintf(coeff,10,"%d",A(i,j));
		mexPrintf("%f ",data(i,j));
		}
		mexPrintf("\n")
		};*/
	}else
	{
		mexErrMsgTxt("params.data must be specified");
	}

	//nbFact initialisation
	int nbFact=0;
	if (presentFields[1])
	{

		mxCurrentField = mxGetField(matlab_params,0,"nfacts");
		nbFact =(int)  mxGetScalar(mxCurrentField);
		//        mexPrintf("NB FACT : %d\n",nbFact);
	}else
	{
		mexErrMsgTxt("params.nfacts must be specified");
	}



	//constraints
	std::vector<const Faust::ConstraintGeneric*> consS;
	if (presentFields[2])
	{
		mwSize nbRowCons,nbColCons;
		mxCurrentField = mxGetField(matlab_params,0,"cons");
		if(!mxIsCell(mxCurrentField))
		{
			mexErrMsgTxt("cons must be a cell-array");
		}
		nbRowCons = mxGetM(mxCurrentField);
		nbColCons = mxGetN(mxCurrentField);

		if(nbRowCons !=1)
		{

			mexErrMsgTxt("cons must have 1 rows");
		}
		//		mexPrintf("cons nbColCons=%d\n", nbColCons);
		if(nbColCons != (nbFact))
		{
			//mexPrintf("\n cons has %d cols and nbFact = %d\n",nbColCons,nbFact);
			//mexErrMsgTxt("incoherence between the number of columns of cons and nfacts ");
		}
		//mexPrintf("\n cons has %d rows and %d cols \n",nbRowCons,nbColCons);
		//Faust::ConstraintGeneric * consToAdd;



		for (mwSize j=0;j<nbColCons;j++)
		{
			//                mexPrintf("cons(%d)\n",j);
			mxCurrentCons=mxGetCell(mxCurrentField,j);
			getConstraint<SCALAR,FPP2>(consS,mxCurrentCons);
			//consS.push_back(consToAdd);
		}

	}else
	{
		mexErrMsgTxt("params.cons must be specified");
	}
	//niter1
	//    Faust::StoppingCriterion<SCALAR> crit1;
	//    if (presentFields[3])
	//    {
	//         mxCurrentField = mxGetField(matlab_params,0,"niter");
	//        int nb_iter1 =(int)  mxGetScalar(mxCurrentField);
	//        Faust::StoppingCriterion<SCALAR> newCrit1(nb_iter1);
	//        crit1 = newCrit1;
	//    }
	//    mexPrintf("\n crit1 nb_it = %d\n",crit1.get_crit());

	//TODO: replace by default values as constants from StoppingCriterion class
	bool is_criterion_error = false;
	int num_its = 500;
	FPP2 error_treshold = 0.3;
	int max_num_its = 10000;
	if(presentFields[3])
	{
		mxCurrentField = mxGetField(matlab_params, 0, "niter");
		num_its = (int) mxGetScalar(mxCurrentField);
	}
	if(presentFields[8]){
		mxCurrentField = mxGetField(matlab_params, 0, "sc_is_criterion_error");
		is_criterion_error =  (bool) mxGetScalar(mxCurrentField);
	}
	if(presentFields[9])
	{
		mxCurrentField = mxGetField(matlab_params, 0, "sc_error_treshold");
		error_treshold = (FPP2) mxGetScalar(mxCurrentField);
	}
	if(presentFields[10])
	{
		mxCurrentField = mxGetField(matlab_params, 0, "sc_max_num_its");
		max_num_its = (int) mxGetScalar(mxCurrentField);
	}
	Faust::StoppingCriterion<FPP2> crit1(num_its, is_criterion_error, error_treshold, max_num_its);
	//	crit1.Display();
	//init_facts
	std::vector<Faust::MatDense<SCALAR,Cpu> > init_facts;
	if (presentFields[4])
	{
		mxCurrentField = mxGetField(matlab_params,0,"init_facts");
		setVectorFaustMat(init_facts,mxCurrentField);

	}else
	{
		mexErrMsgTxt("init_facts must be must be specified");
	}
	//verbosity
	bool isVerbose = false;
	if (presentFields[5])
	{
		mxCurrentField = mxGetField(matlab_params,0,"verbose");
		isVerbose =(bool)  mxGetScalar(mxCurrentField);
	}
	//update_way
	bool updateway = false;
	if (presentFields[7])
	{
		mxCurrentField = mxGetField(matlab_params,0,"update_way");
		updateway =(bool)  mxGetScalar(mxCurrentField);
	}


	//init_lambda
	FPP2 init_lambda = 1.0;
	if (presentFields[6])
	{
		mxCurrentField = mxGetField(matlab_params,0,"init_lambda");
		FPP2* tmp_ptr = &init_lambda;
		// it works whatever mxCurrentField class is (complex or not)
		mxArray2Ptr<FPP2>(const_cast<const mxArray*>(mxCurrentField),tmp_ptr);
		//       init_lambda = (SCALAR) mxGetScalar(mxCurrentField);
	}

	GradientCalcOptMode grad_calc_opt_mode = Params<SCALAR,Cpu, FPP2>::defaultGradCalcOptMode;
	if(presentFields[11])
	{
		mxCurrentField = mxGetField(matlab_params,0,"grad_calc_opt_mode");
		grad_calc_opt_mode = static_cast<GradientCalcOptMode>((int)mxGetScalar(mxCurrentField));
	}

	bool constant_step_size = false;
	if (presentFields[12])
	{
		mxCurrentField = mxGetField(matlab_params,0,"constant_step_size");
		constant_step_size =(bool)  mxGetScalar(mxCurrentField);
	}
	FPP2 step_size = Params<SCALAR,Cpu, FPP2>::defaultStepSize;
	if(presentFields[13])
	{
		mxCurrentField = mxGetField(matlab_params, 0, "step_size");
		step_size = (FPP2) mxGetScalar(mxCurrentField);
	}
	FPP2 norm2_threshold =	FAUST_PRECISION;
	int norm2_max_iter = FAUST_NORM2_MAX_ITER;
	if(presentFields[14])
	{
		mxCurrentField = mxGetField(matlab_params, 0, "norm2_max_iter");
		norm2_max_iter = (int) mxGetScalar(mxCurrentField);
	}
	if(presentFields[15])
	{
		mxCurrentField = mxGetField(matlab_params, 0, "norm2_threshold");
		norm2_threshold = (FPP2) mxGetScalar(mxCurrentField);
	}
	FactorsFormat factors_format = Params<SCALAR, Cpu, FPP2>::defaultFactorsFormat;
	bool packing_RL = Params<SCALAR, Cpu, FPP2>::defaultPackingRL;
	if(presentFields[16])
	{
		mxCurrentField = mxGetField(matlab_params, 0, "factor_format");
		factors_format = static_cast<FactorsFormat>((int)mxGetScalar(mxCurrentField));
	}
	if(presentFields[17])
	{
		mxCurrentField = mxGetField(matlab_params, 0, "packing_RL");
		packing_RL = (bool) mxGetScalar(mxCurrentField);
	}
	bool no_normalization = Params<SCALAR, Cpu, FPP2>::defaultNoNormalization;
	if(presentFields[18])
	{
		mxCurrentField = mxGetField(matlab_params, 0, "no_normalization");
		no_normalization = (bool) mxGetScalar(mxCurrentField);
	}
	bool no_lambda = Params<SCALAR, Cpu, FPP2>::defaultNoLambda;
	if(presentFields[19])
	{
		mxCurrentField = mxGetField(matlab_params, 0, "no_lambda");
		no_lambda = (bool) mxGetScalar(mxCurrentField);
	}
	//compute_lambda
	// bool compute_lambda = true;
	// if (presentFields[8])
	// {
	// mxCurrentField = mxGetField(matlab_params,0,"compute_lambda");
	// compute_lambda = (bool) mxGetScalar(mxCurrentField);
	// }

	params = new Faust::ParamsPalm<SCALAR,Cpu, FPP2>(data,
														nbFact,
														consS,
														init_facts,
														crit1,
														isVerbose,
														updateway,
														init_lambda,
														constant_step_size,
														step_size,
														grad_calc_opt_mode);

	if(norm2_max_iter) params->norm2_max_iter = norm2_max_iter;
	if(norm2_threshold != FPP2(0)) params->norm2_threshold = norm2_threshold;
	params->factors_format = factors_format;
	params->packing_RL = packing_RL;
	params->no_normalization = no_normalization;
	params->no_lambda = no_lambda;

	return params;
}

template<typename SCALAR>
void mxArray2FaustMHTPParams(const mxArray* matlab_params, Faust::MHTPParams<SCALAR>& params)
{
	// all fields are optional
	mxArray *mx_field = mxGetField(matlab_params, 0, "mhtp_num_its");
	if(params.used = (mx_field != nullptr))
		params.sc = Faust::StoppingCriterion<SCALAR>((int) mxGetScalar(mx_field));

	mx_field = mxGetField(matlab_params, 0, "mhtp_constant_step_size");
	if(params.used = (params.used || (mx_field != nullptr)))
		params.constant_step_size = (bool) mxGetScalar(mx_field);

	mx_field = mxGetField(matlab_params, 0, "mhtp_step_size");
	if(params.used = (params.used || (mx_field != nullptr)))
		params.step_size = (Real<SCALAR>) mxGetScalar(mx_field);

	mx_field = mxGetField(matlab_params, 0, "mhtp_palm4msa_period");
	if(params.used = (params.used || (mx_field != nullptr)))
		params.palm4msa_period = (int) mxGetScalar(mx_field);


	mx_field = mxGetField(matlab_params, 0, "mhtp_updating_lambda");
	if(params.used = (params.used || (mx_field != nullptr)))
		params.updating_lambda = (bool) mxGetScalar(mx_field);

}

#endif
