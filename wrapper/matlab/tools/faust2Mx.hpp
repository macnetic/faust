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
		
  
	const mwSize dims[2]={(mwSize)row,(mwSize)col};
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
		
  
	const mwSize dims[2]={(mwSize)row,(mwSize)col};
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



template<typename FPP>
mxArray*  FaustVec2mxArray(const Faust::Vect<FPP,Cpu>& M)
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


	template<typename FPP>
mxArray* transformFact2SparseMxArray(faust_unsigned_int id, Faust::TransformHelper<FPP,Cpu>* core_ptr)
{
	mxArray* sparseMat = mxCreateSparse(core_ptr->get_fact_nb_rows(id),
			core_ptr->get_fact_nb_cols(id),
			core_ptr->get_fact_nnz(id),
			mxREAL);
	mwIndex* ir = mxGetIr(sparseMat);
	mwIndex* jc = mxGetJc(sparseMat);
	//	std::cout << "transformFact2SparseMxArray()" << std::endl;
	//	FPP* pr = static_cast<FPP*>(mxGetDoubles(sparseMat)); //given up because fails to compile with template FPP
	FPP* pr = static_cast<FPP*>(mxGetPr(sparseMat));
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

template<typename FPP>
mxArray* transformFact2SparseMxArray(faust_unsigned_int id, Faust::TransformHelper<complex<FPP>,Cpu>* core_ptr)
{
	mxArray* sparseMat = mxCreateSparse(core_ptr->get_fact_nb_rows(id),
			core_ptr->get_fact_nb_cols(id),
			core_ptr->get_fact_nnz(id),
			mxCOMPLEX);
	mwIndex* ir = mxGetIr(sparseMat);
	mwIndex* jc = mxGetJc(sparseMat);
	FPP*    ptr_real_data = static_cast<FPP*> (mxGetPr(sparseMat));
	FPP*    ptr_imag_data = static_cast<FPP*> (mxGetPi(sparseMat));
	faust_unsigned_int nnz, num_rows, num_cols;
	complex<FPP>* cplxData  = (complex<FPP>*) malloc(sizeof(complex<FPP>)*mxGetNzmax(sparseMat));
	int * i_ir = (int*) malloc(sizeof(int)*mxGetNzmax(sparseMat)); //TODO: use new[]/delete[]
	int * i_jc = (int*) malloc(sizeof(int)*mxGetN(sparseMat)+1);
	core_ptr->get_fact(id, i_jc, i_ir, cplxData, &nnz, &num_rows, &num_cols, true);
	assert(nnz == mxGetNzmax(sparseMat));
	splitComplexPtr(cplxData, nnz, ptr_real_data, ptr_imag_data);
	for(int i=0;i<nnz;i++)
		ir[i] = i_ir[i];
	for(int i=0;i<mxGetN(sparseMat)+1;i++)
		jc[i] = i_jc[i];
	free(i_ir);
	free(i_jc);
	free(cplxData);
	return sparseMat;
}

template<typename FPP>
mxArray* transformFact2FullMxArray(faust_unsigned_int id, Faust::TransformHelper<FPP,Cpu>* core_ptr)
{
	const mwSize dim_sizes[2] = {core_ptr->get_fact_nb_rows(id),
			core_ptr->get_fact_nb_cols(id)};
//	cout << "transformFact2FullMxArray() start" << endl;
	mxClassID classId = typeid(FPP)==typeid(float)?mxSINGLE_CLASS:mxDOUBLE_CLASS;
	mxArray* fullMat = mxCreateNumericArray(2, dim_sizes, classId, mxREAL);
	FPP* data_ptr = static_cast<FPP*>(mxGetData(fullMat));
	faust_unsigned_int num_rows, num_cols;
	//transposition is handled internally
	core_ptr->get_fact(id, data_ptr, &num_rows, &num_cols);
	//data_ptr now contains a copy of id-th factor elements
	//then the matlab matrix is set in just one copy
//	cout << "transformFact2FullMxArray() end" << endl;
	return fullMat;
}


template<typename FPP>
mxArray* transformFact2FullMxArray(faust_unsigned_int id, Faust::TransformHelper<complex<FPP>,Cpu>* core_ptr)
{
	const mwSize dim_sizes[2] = {core_ptr->get_fact_nb_rows(id),
		core_ptr->get_fact_nb_cols(id)};
	mxClassID classId = typeid(FPP)==typeid(float)?mxSINGLE_CLASS:mxDOUBLE_CLASS;
	mxArray* fullMat = mxCreateNumericArray(2, dim_sizes, classId, mxCOMPLEX);
	FPP* data_ptr = static_cast<FPP*>(mxGetData(fullMat));
	FPP* img_data_ptr = static_cast<FPP*>(mxGetImagData(fullMat));
	faust_unsigned_int num_rows, num_cols;
	if(core_ptr->isTransposed())
	{
		// the factor is transposed, we can't directly access the data
		// cause it needs to be reordered
		complex<FPP>* src_data_ptr = (complex<FPP>*) malloc(sizeof(complex<FPP>)*dim_sizes[0]*dim_sizes[1]); // TODO: alloc-free in c++ style (new/delete)
		core_ptr->get_fact(id, src_data_ptr, &num_rows, &num_cols);
		splitComplexPtr(src_data_ptr, num_cols*num_rows, data_ptr, img_data_ptr, core_ptr->isConjugate());
		free(src_data_ptr);
		// it's different to real transformFact2FullMxArray() case here
		// because with complex we need to split im/real (until we stop to use separated complex API)
		// when the split will be removed, we'll do as transformFact2SparseMxArray() for FPP==double/float
		// (handling complex and real the same way)
	}
	else
	{
		const complex<FPP>* src_data_ptr;
		//not calling the same prototype of get_fact() called in transpose case (and real version of this function)
		core_ptr->get_fact(id, &src_data_ptr, &num_rows, &num_cols);
		splitComplexPtr(src_data_ptr, num_rows*num_cols, data_ptr, img_data_ptr, core_ptr->isConjugate());
	}
	//TODO: matlab 2018 setComplexDoubles()/setComplexSingles() to avoid separating
	// in real and imaginary part to copy, directly copy into mxArray with get_fact() in the pointer returned by mxComplexSingles()/mxComplexDoubles()
	return fullMat;
}





#endif
