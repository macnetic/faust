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
#include "matrix.h"






template<typename FPP, FDevice DEV>
mxArray*  FaustMat2mxArray(const Faust::MatDense<FPP,DEV>& M)
{
#ifndef MX_HAS_INTERLEAVED_COMPLEX
	if (!M.isReal())
		mexErrMsgTxt("FaustMat2mxArray : Faust::MatDense must be real");
#endif

	mxArray * mxMat;
	int row,col;
	row = M.getNbRow();
	col = M.getNbCol();

	const mwSize dims[3]={(mwSize)row,(mwSize)col};
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

	FPP*    ptr_out = nullptr;
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	newMxGetData(ptr_out, mxMat);
#else
	ptr_out = static_cast<FPP*> (mxGetData(mxMat));
#endif
//	memcpy(ptr_out, M.getData(), row*col*sizeof(FPP));

	M.copyBuf(ptr_out);

	return mxMat;

}

#ifdef MX_HAS_INTERLEAVED_COMPLEX
template<>
void newMxGetData<float>(float*& ptr_out, const mxArray* mxMat)
{
	if(mxGetClassID(mxMat) != mxSINGLE_CLASS || mxIsComplex(mxMat))
		mexErrMsgTxt("newMxGetData: the mex matrix must be double as the ptr is.");
	ptr_out	= static_cast<float*> (mxGetSingles(mxMat));
}

template<>
void newMxGetData<double>(double*& ptr_out, const mxArray* mxMat)
{
	if(mxGetClassID(mxMat) != mxDOUBLE_CLASS && mxIsComplex(mxMat))
		mexErrMsgTxt("newMxGetData: the mex matrix must be double as the ptr.");
	ptr_out	= static_cast<double*> (mxGetDoubles(mxMat));
}

template<>
void newMxGetData<std::complex<double>>(complex<double>*& ptr_out, const mxArray* mxMat)
{
	if(mxGetClassID(mxMat) != mxDOUBLE_CLASS || ! mxIsComplex(mxMat))
		mexErrMsgTxt("newMxGetData: the mex matrix must be complex double as the ptr.");
	ptr_out	= reinterpret_cast<std::complex<double>*> (mxGetComplexDoubles(mxMat));
}

template<>
void newMxGetData<std::complex<float>>(complex<float>*& ptr_out, const mxArray* mxMat)
{
	if(mxGetClassID(mxMat) != mxSINGLE_CLASS || ! mxIsComplex(mxMat))
		mexErrMsgTxt("newMxGetData: the mex matrix must be complex float as the ptr.");
	ptr_out	= reinterpret_cast<std::complex<float>*> (mxGetComplexDoubles(mxMat));
}
#endif

template<typename FPP>
mxArray* helperCreateNumMxArray(const mwSize dims[2])
{
	mxArray *mxMat;
	if(std::is_same<FPP, complex<float>>::value)
		mxMat = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxCOMPLEX);
	else if(std::is_same<FPP, complex<double>>::value)
		mxMat = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);
	else if(std::is_same<FPP, float>::value)
		mxMat = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
	else if(std::is_same<FPP, double>::value)
		mxMat = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
	else
		mexErrMsgTxt("helperCreateNumMxArray: unsupported scalar type of matrix");
	return mxMat;
}

template<typename FPP>
mxArray*  FaustMat2mxArray(const Faust::MatDense<std::complex<FPP>,Cpu>& M)
{
	mxArray * mxMat;
	int row,col;
	row = M.getNbRow();
        col = M.getNbCol();


	const mwSize dims[2]={(mwSize)row,(mwSize)col};
	mxMat = helperCreateNumMxArray<complex<FPP>>(dims);

#ifdef MX_HAS_INTERLEAVED_COMPLEX
// TODO: refactor this code (it is used many times in this file)
	complex<FPP>* ptr_data;
	if(is_same<FPP, float>::value)
		ptr_data = reinterpret_cast<complex<FPP>*>(mxGetComplexSingles(mxMat));
	else
		ptr_data = reinterpret_cast<complex<FPP>*>(mxGetComplexDoubles(mxMat));
	memcpy(ptr_data, M.getData(), sizeof(complex<FPP>)*row*col);
#else
	FPP*    ptr_real_data = static_cast<FPP*> (mxGetData(mxMat));
	FPP*    ptr_imag_data = static_cast<FPP*> (mxGetImagData(mxMat));

	splitComplexPtr(M.getData(),row*col,ptr_real_data,ptr_imag_data);
#endif
	return mxMat;
}


template<typename FPP>
mxArray*  FaustVec2mxArray(const Faust::Vect<std::complex<FPP>,Cpu>& v)
{
	mxArray * mxMat;
	int row,col;
	row = v.size();
        col = 1;
		
  
	const mwSize dims[2]={(mwSize)row,(mwSize)col};
	mxMat = helperCreateNumMxArray<FPP>(dims);
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	copyComplexDataToMxArray(v.getData(), sizeof(complex<FPP>)*row*col, mxMat);
#else

	FPP*    ptr_real_data = static_cast<FPP*> (mxGetData(mxMat));
	FPP*    ptr_imag_data = static_cast<FPP*> (mxGetImagData(mxMat));
	splitComplexPtr(v.getData(),row*col,ptr_real_data,ptr_imag_data);
#endif

	return mxMat;


}

#ifdef MX_HAS_INTERLEAVED_COMPLEX
template<typename FPP>
void copyComplexDataToMxArray(const complex<FPP> *data, size_t data_sz, mxArray* mxMat, bool conjugate/*=false*/)
{
	complex<FPP>* ptr_data;
	if(is_same<FPP, float>::value)
		ptr_data = reinterpret_cast<complex<FPP>*>(mxGetComplexSingles(mxMat));
	else
		ptr_data = reinterpret_cast<complex<FPP>*>(mxGetComplexDoubles(mxMat));
	if(conjugate)
	{
		for(int i=0;i<data_sz;i++)

			ptr_data[i]  = std::complex<FPP>(std::real(data[i]), - std::imag(data[i]));
	}
	else
		memcpy(ptr_data, data, sizeof(complex<FPP>)*data_sz);
}
#endif



#ifndef MX_HAS_INTERLEAVED_COMPLEX
template<typename FPP>
void splitComplexPtr(const std::complex<FPP>*  cpx_ptr, int nb_element, FPP* & real_ptr, FPP* & imag_ptr, const bool conjugate)
{
	std::complex<FPP> cpxValue;
	for (int i=0;i<nb_element;i++)
	{
		cpxValue=cpx_ptr[i];
		real_ptr[i] = cpxValue.real();
		imag_ptr[i] = conjugate?-cpxValue.imag():cpxValue.imag();
	}
}
#endif



template<typename FPP>
mxArray* FaustVec2mxArray(const Faust::Vect<FPP,Cpu>& M)
{
	if (!M.isReal())
		mexErrMsgTxt("FaustMat2mxArray : Faust::MatDense must be real");

	mxArray * mxMat;
	int row,col;
	row = M.size();
	col = 1;

	const mwSize dims[2]={(mwSize)row,(mwSize)col};
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

	FPP* ptr_out;
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	newMxGetData(ptr_out, mxMat);
#else
	ptr_out = static_cast<FPP*> (mxGetData(mxMat));
#endif
	memcpy(ptr_out, M.getData(), row*col*sizeof(FPP));
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

template<typename FPP>
mxArray*  FaustSpMat2mxArray(const Faust::MatSparse<FPP,Cpu>& M)
{
	faust_unsigned_int nnz = M.getNonZeros();
	mxArray* sparseMat = mxCreateSparse(M.getNbRow(),
			M.getNbCol(),
			nnz,
			mxREAL);
	mwIndex* ir = mxGetIr(sparseMat);
	mwIndex* jc = mxGetJc(sparseMat);
	FPP* pr;
	//TODO: and complex case ?
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	pr = static_cast<FPP*>(mxGetDoubles(sparseMat));
#else
	pr = static_cast<FPP*>(mxGetPr(sparseMat));
#endif
	// sadly we can't do a direct copy into ir and jc because MatSparse uses int type for indices
	// (not mwIndex which is in fact a size_t)
	// so we need to copy in intermediate buffers and then affect their elements
	// into jc, ir
	auto tM = M;
	// transpose M because matlab is in CSC format while Faust is in CSR
	tM.transpose();
//	tM.copyRowPtr(jc);
//	tM.copyColInd(ir);
//	tM.copyValuePtr(pr);
	tM.copyBufs(jc, ir, pr);
	return sparseMat;
}

	template<FDevice DEV>
mxArray* transformFact2SparseMxArray(faust_unsigned_int id, Faust::TransformHelper<float, DEV>* core_ptr)
{
	mxArray* sparseMat = mxCreateSparse(core_ptr->get_fact_nb_rows(id),
			core_ptr->get_fact_nb_cols(id),
			core_ptr->get_fact_nnz(id),
			mxREAL);
	mwIndex* ir = mxGetIr(sparseMat);
	mwIndex* jc = mxGetJc(sparseMat);
	//	std::cout << "transformFact2SparseMxArray()" << std::endl;
	//	FPP* pr = static_cast<FPP*>(mxGetDoubles(sparseMat)); //given up because fails to compile with template FPP
	double* pr;
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	pr = static_cast<double*>(mxGetDoubles(sparseMat));
#else
	pr = static_cast<double*>(mxGetPr(sparseMat));
#endif
	faust_unsigned_int nnz, num_rows, num_cols;
	// sadly we can't do a direct copy into ir and jc because MatSparse uses int type for indices
	// (not mwIndex which is in fact a size_t)
	// so we need to copy in intermediate buffers and then affect their elements
	// into jc, ir
	int * i_ir = (int*) malloc(sizeof(int)*mxGetNzmax(sparseMat)); //TODO: use new[]/delete[]
	int * i_jc = (int*) malloc(sizeof(int)*(mxGetN(sparseMat)+1));
	// last arg. transpose == true, because Faust works with csr/row-major order matrices
	// while matlab works with csc/column-major order matrices, so we ask transposition to obtain the proper repr.
	// the case when the Faust is transposed is handled internally (by not transposing the buffers in csr to get transpose in csc)
	// the reordering operation costs additional calculation time
	//TODO: remove nnz, num_rows, num_cols when get_fact() will authorize NULL
	float* flt_pr = new float[mxGetNzmax(sparseMat)];
	core_ptr->get_fact(id, i_jc, i_ir, flt_pr, &nnz, &num_rows, &num_cols, true);
	for(int i=0;i<nnz;i++)
	{
		ir[i] = (mwIndex)i_ir[i];
		//		cout << ir[i] << " ";
		pr[i] = static_cast<double>(flt_pr[i]);
	}
	delete [] flt_pr;
	//	cout << endl;
	//	for(int i=0;i<nnz;i++)
	//	{
	//		cout << pr[i] << " ";
	//	}
	//	cout << endl;
	for(int i=0;i<mxGetN(sparseMat)+1;i++)
	{
		jc[i] = (mwIndex)i_jc[i];
		//		cout << jc[i] << " ";
	}
	//	cout << endl;
	//	std::cout << "transformFact2SparseMxArray()" << std::endl;
	free(i_ir);
	free(i_jc);
	return sparseMat;
}
	template<typename FPP, FDevice DEV>
mxArray* transformFact2SparseMxArray(faust_unsigned_int id, Faust::TransformHelper<FPP, DEV>* core_ptr)
{
	mxArray* sparseMat = mxCreateSparse(core_ptr->get_fact_nb_rows(id),
			core_ptr->get_fact_nb_cols(id),
			core_ptr->get_fact_nnz(id),
			mxREAL);
	mwIndex* ir = mxGetIr(sparseMat);
	mwIndex* jc = mxGetJc(sparseMat);
	//	std::cout << "transformFact2SparseMxArray()" << std::endl;
	//	FPP* pr = static_cast<FPP*>(mxGetDoubles(sparseMat)); //given up because fails to compile with template FPP
	FPP* pr = nullptr;
#if MX_HAS_INTERLEAVED_COMPLEX
	pr = static_cast<FPP*>(mxGetDoubles(sparseMat));
#else
	pr = static_cast<FPP*>(mxGetPr(sparseMat));
#endif
	faust_unsigned_int nnz, num_rows, num_cols;
	// sadly we can't do a direct copy into ir and jc because MatSparse uses int type for indices
	// (not mwIndex which is in fact a size_t)
	// so we need to copy in intermediate buffers and then affect their elements
	// into jc, ir
	int * i_ir = (int*) malloc(sizeof(int)*mxGetNzmax(sparseMat)); //TODO: use new[]/delete[]
	int * i_jc = (int*) malloc(sizeof(int)*(mxGetN(sparseMat)+1));
	// last arg. transpose == true, because Faust works with csr/row-major order matrices
	// while matlab works with csc/column-major order matrices, so we ask transposition to obtain the proper repr.
	// the case when the Faust is transposed is handled internally (by not transposing the buffers in csr to get transpose in csc)
	// the reordering operation costs additional calculation time
	//TODO: remove nnz, num_rows, num_cols when get_fact() will authorize NULL
	core_ptr->get_fact(id, i_jc, i_ir, pr, &nnz, &num_rows, &num_cols, true);
	for(int i=0;i<nnz;i++)
	{
		ir[i] = (mwIndex)i_ir[i];
		//		cout << ir[i] << " ";
	}
	//	cout << endl;
	//	for(int i=0;i<nnz;i++)
	//	{
	//		cout << pr[i] << " ";
	//	}
	//	cout << endl;
	for(int i=0;i<mxGetN(sparseMat)+1;i++)
	{
		jc[i] = (mwIndex)i_jc[i];
		//		cout << jc[i] << " ";
	}
	//	cout << endl;
	//	std::cout << "transformFact2SparseMxArray()" << std::endl;
	free(i_ir);
	free(i_jc);
	return sparseMat;
}

template<typename FPP, FDevice DEV>
mxArray* transformFact2SparseMxArray(faust_unsigned_int id, Faust::TransformHelper<complex<FPP>,DEV>* core_ptr)
{
	faust_unsigned_int nnz, num_rows, num_cols;
	mxArray* sparseMat = mxCreateSparse(core_ptr->get_fact_nb_rows(id),
			core_ptr->get_fact_nb_cols(id),
			core_ptr->get_fact_nnz(id),
			mxCOMPLEX);
	mwIndex* ir = mxGetIr(sparseMat);
	mwIndex* jc = mxGetJc(sparseMat);
	int * i_ir = (int*) malloc(sizeof(int)*mxGetNzmax(sparseMat)); //TODO: use new[]/delete[]
	int * i_jc = (int*) malloc(sizeof(int)*mxGetN(sparseMat)+1);
	complex<FPP>* cplxData  = (complex<FPP>*) malloc(sizeof(complex<FPP>)*mxGetNzmax(sparseMat));
	core_ptr->get_fact(id, i_jc, i_ir, cplxData, &nnz, &num_rows, &num_cols, true);
	assert(nnz == mxGetNzmax(sparseMat));
#if MX_HAS_INTERLEAVED_COMPLEX
	copyComplexDataToMxArray(cplxData, sizeof(complex<FPP>)*nnz, sparseMat);
#else
	FPP*    ptr_real_data = static_cast<FPP*> (mxGetPr(sparseMat));
	FPP*    ptr_imag_data = static_cast<FPP*> (mxGetPi(sparseMat));
	splitComplexPtr(cplxData, nnz, ptr_real_data, ptr_imag_data);
#endif
	for(int i=0;i<nnz;i++)
		ir[i] = i_ir[i];
	for(int i=0;i<mxGetN(sparseMat)+1;i++)
		jc[i] = i_jc[i];
	free(i_ir);
	free(i_jc);
	free(cplxData);
	return sparseMat;
}

template<typename FPP, FDevice DEV>
mxArray* transformFact2FullMxArray(faust_unsigned_int id, Faust::TransformHelper<FPP,DEV>* core_ptr)
{
	const mwSize dim_sizes[2] = {core_ptr->get_fact_nb_rows(id),
			core_ptr->get_fact_nb_cols(id)};
//	cout << "transformFact2FullMxArray() start" << endl;
	mxClassID classId = typeid(FPP)==typeid(float)?mxSINGLE_CLASS:mxDOUBLE_CLASS;
	mxArray* fullMat = mxCreateNumericArray(2, dim_sizes, classId, mxREAL);
	FPP* data_ptr;
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	newMxGetData(data_ptr, fullMat);
#else
	data_ptr = static_cast<FPP*>(mxGetData(fullMat));
#endif
	faust_unsigned_int num_rows, num_cols;
	//transposition is handled internally
	core_ptr->get_fact(id, data_ptr, &num_rows, &num_cols);
	//data_ptr now contains a copy of id-th factor elements
	//then the matlab matrix is set in just one copy
//	cout << "transformFact2FullMxArray() end" << endl;
	return fullMat;
}


template<typename FPP, FDevice DEV>
mxArray* transformFact2FullMxArray(faust_unsigned_int id, Faust::TransformHelper<complex<FPP>,DEV>* core_ptr)
{
	const mwSize dim_sizes[2] = {core_ptr->get_fact_nb_rows(id),
		core_ptr->get_fact_nb_cols(id)};
	mxClassID classId = typeid(FPP)==typeid(float)?mxSINGLE_CLASS:mxDOUBLE_CLASS;
	mxArray* fullMat = mxCreateNumericArray(2, dim_sizes, classId, mxCOMPLEX);
#ifndef MX_HAS_INTERLEAVED_COMPLEX
	FPP* data_ptr = static_cast<FPP*>(mxGetData(fullMat));
	FPP* img_data_ptr = static_cast<FPP*>(mxGetImagData(fullMat));
#endif
	faust_unsigned_int num_rows, num_cols;
	if(core_ptr->isTransposed())
	{
		// the factor is transposed, we can't directly access the data
		// cause it needs to be reordered
		complex<FPP>* src_data_ptr = (complex<FPP>*) malloc(sizeof(complex<FPP>)*dim_sizes[0]*dim_sizes[1]); // TODO: alloc-free in c++ style (new/delete)
		core_ptr->get_fact(id, src_data_ptr, &num_rows, &num_cols);
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		copyComplexDataToMxArray(src_data_ptr, num_cols*num_rows, fullMat);
#else
		splitComplexPtr(src_data_ptr, num_cols*num_rows, data_ptr, img_data_ptr/*, core_ptr->isConjugate()*/); //let conjugate to false because it's handled internally by get_fact()
		// it's different to real transformFact2FullMxArray() case here
		// because with complex we need to split im/real (until we stop to use separated complex API)
		// when the split will be removed, we'll do as transformFact2SparseMxArray() for FPP==double/float
		// (handling complex and real the same way)
#endif
		free(src_data_ptr);
	}
	else
	{
		complex<FPP>* src_data_ptr;
		//not calling the same prototype of get_fact() called in transpose case (and real version of this function)
		core_ptr->get_fact(id, src_data_ptr, &num_rows, &num_cols);
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		copyComplexDataToMxArray(src_data_ptr, num_rows*num_cols, fullMat, core_ptr->isConjugate());
#else
		splitComplexPtr(src_data_ptr, num_rows*num_cols, data_ptr, img_data_ptr, core_ptr->isConjugate());//let conjugate to false because it's handled internally by get_fact()
#endif

	}
	//TODO: matlab 2018 setComplexDoubles()/setComplexSingles() to avoid separating
	// in real and imaginary part to copy, directly copy into mxArray with get_fact() in the pointer returned by mxComplexSingles()/mxComplexDoubles()
	return fullMat;
}

template<class FPP>
void  mxArray2Scalar(const mxArray* scalar, typename std::enable_if<std::is_floating_point<FPP>::value, FPP>::type* out)
{
//	cout << "mxArray2Scalar() FPP real" << endl;
	*out = mxGetScalar(scalar);
}

template<typename FPP>
void mxArray2Scalar(const mxArray* scalar, complex<FPP>* out)
{
	double real_part = mxGetScalar(scalar);
	complex<FPP> cplx;
	if(mxIsComplex(scalar))
	{
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		if(is_same<FPP, float>::value)
		{
			mxComplexSingle *data = mxGetComplexSingles(scalar);
			cplx = complex<FPP>(real_part, data->imag);
		}
		else
		{
			mxComplexDouble *data = mxGetComplexDoubles(scalar);
			cplx = complex<FPP>(real_part, data->imag);
		}
#else
		double* im_ptr;
		im_ptr = (FPP*) mxGetPi(scalar);
//		cout << "im_ptr[0]: " << im_ptr[0] << endl;
		cplx = complex<FPP>(real_part, im_ptr[0]);
#endif
	}
	else
		cplx = complex<FPP>(real_part);
	*out = cplx;
}


#endif
