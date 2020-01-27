/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2019):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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
/*      Hakim H. hakim.hadj-djilani@inria.fr                                */
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
/****************************************************************************/
#ifndef __FAUSTCORE_HPP__
#define __FAUSTCORE_HPP__


#include "faust_Vect.h"
#include <complex>
//#include "faust_HierarchicalFact.h"
//#include "faust_Params.h"
#include "faust_linear_algebra.h"
#include "faust_transform_algebra.h"
#include <iostream>
#include "faust_exception.h"
#include <fstream>
#include "faust_BlasHandle.h"
#include "faust_SpBlasHandle.h"
#include "matio.h"




template<typename FPP>
void Faust::Transform<FPP,Cpu>::faust_gemm(const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char  typeA, char  typeB) const
{
	handleError(m_className,"faust_gemm : not yet implemented");
	/*
	   faust_unsigned_int nbRowOpA,nbRowOpB,nbColOpA,nbColOpB;

	   if (size() == 0)
	   {
	   handleError(m_className,"faust_gemm : empty Faust::Transform");
	   }
	   if (typeA == 'T')
	   {
	   nbRowOpA = this->getNbCol();
	   nbColOpA = this->getNbRow();
	   }else
	   {
	   nbRowOpA = this->getNbRow();
	   nbColOpA = this->getNbCol();
	   }


	   if (typeB == 'T')
	   {
	   nbRowOpB = B.getNbCol();
	   nbColOpB = B.getNbRow();
	   }else
	   {
	   nbRowOpB = B.getNbRow();
	   nbColOpB = B.getNbCol();
	   }


	   if (nbColOpA != nbRowOpB)
	   {
	   handleError(m_className, "faust_gemm : dimension conflict  between matrix op((*this)) and matrix op(B)");

	   }

	   int* ind_ptr = new int[size()];
	   for (int j=0 ; j<size(); j++)
	   {
	   if (typeA == 'T')
	   ind_ptr[j] = j;
	   else
	   ind_ptr[j] =size()-1-j;
	   }

	   if (beta == 0)
	   {
	   data[ind_ptr[0]].faust_gemm(B,C,alpha,beta,typeA,typeB);
	   if (size() > 1)
	   {
	   Faust::MatDense<FPP,Cpu> tmp1(C);
	   for (int i=1;i<size();i++)
	   {
	   data[ind_ptr[i]].faust_gemm(tmp1,C,1.0,0.0,typeA,'N');
	   tmp1=C;
	   }
	   }

	   }else
	   {
	   if (( (C.getNbRow() != nbRowOpA)	|| (C.getNbCol() != nbColOpB) ) )
	   {
	   handleError(m_className, "faust_gemm : invalid dimension for output matrix C");
	   }


	   if (size() ==1)
	   {
	   data[0].faust_gemm(B,C,alpha,beta,typeA,typeB);
	   }else
	   {
	   Faust::MatDense<FPP,Cpu> tmp2,tmp1;
	data[ind_ptr[0]].faust_gemm(B,tmp1,alpha,0.0,typeA,typeB);

	for (int i=1;i<(size()-1);i++)
	{
		data[ind_ptr[i]].faust_gemm(tmp1,tmp2,1.0,0.0,typeA,'N');
		tmp1=tmp2;
	}

	data[ind_ptr[size()-1]].faust_gemm(tmp1,C,1.0,beta,typeA,'N');
}

}
delete[] ind_ptr;
ind_ptr = NULL;
*/
}


template<typename FPP>
const char * Faust::Transform<FPP,Cpu>::m_className="Faust::Transform<FPP,Cpu>";


template<typename FPP>
Faust::Transform<FPP,Cpu>::Transform() :
	data(std::vector<Faust::MatGeneric<FPP,Cpu>*>()),
	totalNonZeros(0), dtor_delete_data(false)
{

#ifdef __COMPILE_TIMERS__
	t_multiply_vector.resize(0);
#endif

}


template<typename FPP>
Faust::Transform<FPP,Cpu>::Transform(const Faust::Transform<FPP,Cpu> & A) :
	data(std::vector<Faust::MatGeneric<FPP,Cpu>*>()),
	totalNonZeros(A.totalNonZeros), dtor_delete_data(false)
{
	//std::cout<<"debut inside copy constructor"<<std::endl;
	data.resize(A.size());
	for (int i=0 ; i<data.size() ; i++)
	{
		//std::cout<<"iter "<<i<<"/"<<data.size()<<std::endl;
		data[i]=A.data[i]->Clone();
		if(!dtor_delete_data) ref_man.acquire(data[i]);
	}
	//std::cout<<"fin inside copy constructor"<<std::endl;
#ifdef __COMPILE_TIMERS__
	this->t_multiply_vector.resize(data.size());
#endif

}

template<typename FPP>
Faust::Transform<FPP,Cpu>::Transform(const std::vector<Faust::MatGeneric<FPP,Cpu> *> & facts, const FPP lambda_ /*default value = 1.0 */,const bool optimizedCopy /*default value = false*/, const bool cloning_fact /* default to true */) :
	data(std::vector<Faust::MatGeneric<FPP,Cpu>*>()),
	totalNonZeros(0), dtor_delete_data(false)
{
	data.resize(facts.size());
	if(data.size() > 0) {
		if(cloning_fact || lambda_ != FPP(1.0))
		{
			data[0] = facts[0]->Clone(optimizedCopy);
		}
		else
			data[0] = facts[0];
		totalNonZeros += data[0]->getNonZeros();
		for (int i=1 ; i<data.size() ; i++)
		{
			if(cloning_fact)
				data[i]=facts[i]->Clone(optimizedCopy);
			else
				data[i] = facts[i];
			totalNonZeros += data[i]->getNonZeros();
			if(!dtor_delete_data) ref_man.acquire(data[i]);
		}

		if(lambda_ != FPP(1.0))
			(*data[0]) *= lambda_;

		if(!dtor_delete_data) ref_man.acquire(data[0]);

	}
#ifdef __COMPILE_TIMERS__
	this->t_multiply_vector.resize(data.size());
#endif
	this->check_factors_validity();
}


template<typename FPP>
Faust::Transform<FPP,Cpu>::Transform(Faust::Transform<FPP,Cpu>&& T) : totalNonZeros(T.totalNonZeros), dtor_delete_data(T.dtor_delete_data)
{
	data = std::move(T.data);
	totalNonZeros = T.totalNonZeros;
	T.data.resize(0);
}




template<typename FPP>
Faust::Transform<FPP,Cpu>::Transform(const std::vector<Faust::MatDense<FPP,Cpu> >&facts, const bool optimizedCopy /*default value = false*/ ):	data(std::vector<Faust::MatGeneric<FPP,Cpu>*>()),
	totalNonZeros(0), dtor_delete_data(false)
{
	data.resize(facts.size());
	for (int i=0 ; i<data.size() ; i++)
	{
		data[i]=facts[i].Clone(optimizedCopy);
		if(!dtor_delete_data) ref_man.acquire(data[i]);
		totalNonZeros += data[i]->getNonZeros();
	}

}

template<typename FPP>
Faust::Transform<FPP,Cpu>::Transform(const std::vector<Faust::MatSparse<FPP,Cpu> >& facts, const bool optimizedCopy /*default value = false*/ ):	data(std::vector<Faust::MatGeneric<FPP,Cpu>*>()),
	totalNonZeros(0), dtor_delete_data(false)
{
	data.resize(facts.size());
	for (int i=0 ; i<data.size() ; i++)
	{
		data[i]=facts[i].Clone(optimizedCopy);
		totalNonZeros += data[i]->getNonZeros();
		if(!dtor_delete_data) ref_man.acquire(data[i]);
	}
}

template<typename FPP>
Faust::Transform<FPP,Cpu>::Transform(const Transform<FPP, Cpu>* A, const bool transpose_A, const bool conj_A, const Transform<FPP, Cpu>* B, const bool transpose_B, const bool conj_B):
	data(std::vector<Faust::MatGeneric<FPP,Cpu>*>()), totalNonZeros(0), dtor_delete_data(false)
{
	data.resize(A->size()+B->size());
	int i = transpose_A?A->size()-1:0;
	int j = 0;
	// verification that mul. is defined between A and B is done afterward by check_factors_validity()
	bool copying = A->size()>0;
	//TODO: factorize the two loops
	while(copying)
	{
		data[j] = A->data[i]->Clone(false);
		if(transpose_A)
		{
			data[j]->transpose();
			i--;
			copying = i >= 0;
		}
		else
		{
			i++;
			copying = i < A->size();
		}
		if(conj_A) data[j]->conjugate();
		totalNonZeros += data[j]->getNonZeros();
		if(!dtor_delete_data) ref_man.acquire(data[j]);
		j++;
	}
	i = transpose_B?B->size()-1:0;
	copying = B->size()>0;
	while(copying)
	{
		data[j] = B->data[i]->Clone(false);
		if(transpose_B)
		{
			data[j]->transpose();
			i--;
			copying = i >= 0;
		}
		else
		{
			i++;
			copying = i < B->size();
		}
		if(conj_B) data[j]->conjugate();
		totalNonZeros += data[j]->getNonZeros();
		if(!dtor_delete_data) ref_man.acquire(data[j]);
		j++;
	}
	this->check_factors_validity();
}

template<typename FPP>
void Faust::Transform<FPP,Cpu>::check_factors_validity() const
{

	if (size() > 1)
	{

		for (int i=0;i<=size()-2;i++)
		{
			if (data[i]->getNbCol() != data[i+1]->getNbRow())
				handleError(m_className,"check_factors_validity : dimensions of the factors mismatch");

		}


	}

}


template<typename FPP>
void Faust::Transform<FPP,Cpu>::update_total_nnz()
{
	this->totalNonZeros = 0;
	for(auto fac : data)
	{
		totalNonZeros += fac->getNonZeros();
	}
}


template<typename FPP>
Faust::RefManager Faust::Transform<FPP,Cpu>::ref_man(Faust::Transform<FPP,Cpu>::delete_fact);

template<typename FPP>
void Faust::Transform<FPP,Cpu>::print_file(const char* filename) const
{
	if (size() > 0)
	{
		ofstream fichier;
		fichier.open(filename);
		fichier<<size()<<endl<<endl;
		fichier.close();

		for (int i=0;i<size();i++)
		{
			data[i].print_file(filename,std::fstream::app);
		}
	}
}

template<typename FPP>
void Faust::Transform<FPP, Cpu>::save_mat_file(const char* filename, bool transpose, bool conjugate) const
{
	// save the FAuST as a matlab cell array (respecting what is done in matlab wrapper)
	matvar_t *faust_matvar;
	size_t dims[2];
	int i, i2, ret;
	mat_t *matfp;
	if(!size()) handleError("Faust::Transform", "save_mat_file(): can't save an empty Faust.");
	matvar_t **faust_factor_matvars = new matvar_t*[size()];
	for(i=0; i < size(); i++){
		// revert factors order if we are in transpose case
		i2 = transpose?size()-i-1:i;
		faust_factor_matvars[i] = data[i2]->toMatIOVar(transpose, conjugate);
		if(faust_factor_matvars[i] == NULL)
			handleError("Faust::Transform", "save_mat_file(): failed to create i-th factor MatIO variable");
	}
	// write the faust cell array
	dims[0] = 1;
	dims[1] = size();
	faust_matvar = Mat_VarCreate("faust_factors", MAT_C_CELL, MAT_T_CELL, 2, dims,
			faust_factor_matvars, MAT_F_DONT_COPY_DATA);
	if(faust_matvar == NULL)
		handleError("Faust:Transform::save_mat_file()", "Failed to create FAuST MatIO variable");
	matfp = Mat_CreateVer(filename, NULL, MAT_FT_MAT5);
	if(matfp == NULL)
		handleError("Faust::Transform::save_mat_file()", "Failed creating file");
	//	Mat_VarPrint(faust_factor_matvars[1], 1);
	ret = Mat_VarWrite(matfp, faust_matvar, MAT_COMPRESSION_NONE); //TODO: enable compression ?
	//	Mat_VarPrint(faust_matvar,1);
	if(ret)
		handleError("Faust::Transform::save_mat_file()", "Failed writing the FAuST to Matlab file.");
	for(i=0; i < size(); i++)
		Mat_VarFree(faust_factor_matvars[i]);
	// if we didn't use MAT_F_DONT_COPY_DATA flag above we would get a double-free corruption error
	// because it'd have freed also data variable (here the faust factor underlying matrices)
	Mat_VarFree(faust_matvar);
	Mat_Close(matfp);
	delete[] faust_factor_matvars;
}

	template<typename FPP>
void Faust::Transform<FPP,Cpu>::scalarMultiply(const FPP scalar)
{
	if (size() > 0)
		*(data[0])*=scalar;
	else
		handleError(m_className,"scalarMultiply : empty faust can't be multiplied ");
}

	template<typename FPP>
void Faust::Transform<FPP,Cpu>::init_from_file(const char* filename)
{
	clear();
	FILE* fp=fopen(filename,"r");
	int size_tmp;
	fscanf(fp,"%d\n", &size_tmp);
	fscanf(fp,"\n");
	Faust::MatSparse<FPP,Cpu> spmat;
	std::cout<<"size_tmp"<<size_tmp<<std::endl;
	for (int i=0;i<size_tmp;i++)
	{
		spmat.init_from_file(fp);
		data.push_back(spmat);
	}
	fscanf(fp,"\n");
	updateNonZeros();

#ifdef __COMPILE_TIMERS__
	this.t_multiply_vector.resize(data.size());
#endif
}

	template<typename FPP>
void Faust::Transform<FPP,Cpu>::updateNonZeros()
{
	int totalNonZeros_tmp=0;
	for (int i=0;i<size();i++)
		totalNonZeros_tmp+=data[i].getNonZeros();
	totalNonZeros = totalNonZeros_tmp;
}


template<typename FPP>
Faust::MatDense<FPP,Cpu> Faust::Transform<FPP,Cpu>::get_product(const char opThis, const bool isConj)const
{


	if (size() == 0)
	{
		handleError(m_className,"get_product : empty Faust::Transform");
	}

	// modif NB v1102 : new method compatible with factor as MatGeneric
	faust_unsigned_int dim;
	if (opThis == 'T' || opThis == 'H')
		dim = data[0]->getNbRow();
	else
		dim = data[size()-1]->getNbCol();

	Faust::MatDense<FPP,Cpu> prod(dim);

	prod.setEyes();

	Faust::MatDense<FPP,Cpu> p = this->multiply(prod, opThis);

	if(isConj && opThis != 'H') p.conjugate();

	return p;


	/* modif NB v1102 : factor are no longer MatSparse, they are MatGeneric now
	   Faust::MatDense<FPP,Cpu> prod(data[0].getNbRow());

	   if(getNbRow()<getNbCol())
	   {
	   prod.resize(getNbRow());
	   prod.setEyes();
	   for(int i=0 ; i<data.size() ; i++)
	   prod *= data[i];
	   }else
	   {
	   prod.resize(getNbCol());
	   prod.setEyes();
	   for(int i=data.size()-1 ; i>=0 ; i--)
	   prod.multiplyLeft(data[i]);
	   }
	   if (opThis == 'T')
	   prod.transpose();
	   */
	//complexity of evaluating a Faust::Transform
	// from left to right is (dim1*total_nnz)
	// from right to left is (dim2*total_nnz)

	/*Faust::MatDense<FPP,Cpu> prod;
	  if ( (data[0].getNonZeros()+getNbRow()*(totalNonZeros-data[0].getNonZeros())) < (data[size()-1].getNonZeros()+getNbCol()*(totalNonZeros-data[size()-1].getNonZeros())) )
	  {
	  prod = data[0];
	  for(int i=1 ; i<data.size() ; i++)
	  prod *= data[i];
	  }else
	  {
	  prod = data[size()-1];
	  for(int i=data.size()-2 ; i>=0 ; i--)
	  prod.multiplyLeft(data[i]);
	  }	*/



	//return prod;
}

template<typename FPP>
Faust::MatDense<FPP,Cpu> Faust::Transform<FPP,Cpu>::get_product(Faust::BlasHandle<Cpu> blas_handle, Faust::SpBlasHandle<Cpu> spblas_handle)const
{
	return (*this).get_product();
}


template<typename FPP>
faust_unsigned_int Faust::Transform<FPP,Cpu>::getNbRow() const
{
	if (size() != 0)
	{
		return data[0]->getNbRow();
	}else
	{
		return 0;
	}

}

template<typename FPP>
faust_unsigned_int Faust::Transform<FPP,Cpu>::getNbCol() const
{
	if (size() != 0)
	{
		return data[size()-1]->getNbCol();
	}else
	{
		return 0;
	}
}

template<typename FPP>
double Faust::Transform<FPP,Cpu>::spectralNorm(const int nbr_iter_max, double threshold, int &flag) const
{
	if (size() == 0)
	{
		return 1;
	}else
	{
		//std::cout<<"Faust debut spectralNorm"<<std::endl;
		//std::cout<<"copy constructor"<<std::endl;
		Faust::Transform<FPP,Cpu> AtA((*this)); // modif AL <FPP,Cpu>
		//std::cout<<"transposition"<<std::endl;
		AtA.transpose();
		if (getNbCol() < getNbRow())
		{
			AtA.multiply((*this));
		}else
		{
			AtA.multiplyLeft((*this));
		}
		//std::cout<<"Faust fin spectralNorm"<<std::endl;
		FPP maxAbsValue = std::sqrt(power_iteration(AtA,nbr_iter_max,threshold,flag));
		return absValue(maxAbsValue);

	}
}


template<typename FPP>
double Faust::Transform<FPP,Cpu>::normL1(const bool transpose /* = false */) const
{
	double norm;
	MatDense<FPP, Cpu> full = get_product(transpose?'T':'N');
	norm = std::abs(full.normL1(/*transpose*/)); //transpose not necessary because full is already transposed if needed
	return norm;
}

template<typename FPP>
double Faust::Transform<FPP,Cpu>::normL1(const char opThis, const bool isConj) const
{
	double norm;
	MatDense<FPP, Cpu> fact = get_product(opThis, isConj);
	norm = std::abs(fact.normL1(opThis=='T'?true:false));
	return norm;
}

template<typename FPP>
double Faust::Transform<FPP,Cpu>::normFro(const char opThis, const bool isConj) const
{
	double norm;
	MatDense<FPP, Cpu> fact = get_product(opThis, isConj);
	norm = std::abs(fact.norm());
	return norm;
}

template<typename FPP>
double Faust::Transform<FPP,Cpu>::normFro() const
{
	double norm;
	MatDense<FPP, Cpu> fact = get_product();
	norm = std::abs(fact.norm());
	return norm;
}

	template<typename FPP>
void Faust::Transform<FPP,Cpu>::operator=(const Transform<FPP,Cpu>&  f)
{
	for (int i=0;i<size();i++)
	{
		if(dtor_delete_data)
			delete data[i];
		else
		{
			ref_man.release(data[i]);
		}
	}

	data.resize(0);
	data.resize(f.size());

	for (int i=0;i<f.size();i++)
	{
		data[i]=f.data[i]->Clone();
		if(!dtor_delete_data) ref_man.acquire(data[i]);
	}

	this->totalNonZeros=f.totalNonZeros;

}




template<typename FPP>
Faust::Transform<FPP,Cpu>& Faust::Transform<FPP,Cpu>::operator=(Faust::Transform<FPP,Cpu>&& T)
{
	data = std::move(T.data);
	totalNonZeros = T.totalNonZeros;
	dtor_delete_data = T.dtor_delete_data;
	T.data.resize(0);
}

	template<typename FPP>
void Faust::Transform<FPP,Cpu>::multiply(const Faust::Transform<FPP,Cpu> & A) // modif AL <FPP,Cpu>
{
	if (A.size() == 0)
	{

	}else
	{
		if (size() == 0)
		{
			(*this)=A;
		}
		else
		{
			if (getNbCol() != A.getNbRow())
			{
				handleError(m_className,"multiply : dimensions of the 2 faust_transform are in conflict");
			}
			//data.insert(data.end(),A.data.begin(),A.data.end());totalNonZeros+=A.totalNonZeros;
			for (int i=0;i<A.size();i++)
			{
				this->push_back(A.data[i]);
			}
		}
	}
#ifdef __COMPILE_TIMERS__
	this->t_multiply_vector.resize(data.size());
#endif
}


	template<typename FPP>
void Faust::Transform<FPP,Cpu>::multiplyLeft(const Faust::Transform<FPP,Cpu> & A) // modif AL <FPP,Cpu>
{
	if (A.size() == 0)
	{

	}else
	{
		if (size() == 0)
		{
			(*this)=A;
		}
		else
		{
			if (getNbRow() != A.getNbCol())
			{
				handleError(m_className,"multiplyLeft : dimensions of the 2 faustcore are in conflict");
			}
			//data.insert(data.begin(),A.data.begin(),A.data.end());totalNonZeros+=A.totalNonZeros;
			for (int i=A.size()-1;i>=0;i--)
			{
				this->push_first(A.data[i]);
			}

		}
	}
}



template<typename FPP>
Faust::MatGeneric<FPP,Cpu>* Faust::Transform<FPP,Cpu>::get_fact(faust_unsigned_int id, const bool cloning_fact /* default to true */) const
{
	if(id>=size())
	{
		handleError(m_className,"get_fact : id exceeds Faust::Transform size or id < 0");
	}

	if(cloning_fact)
		return data[id]->Clone();
	return data[id];

	//handleError(m_className,"get_fact : not anymore implemented");
	// v1105 : method not anymore compatible since Faust::Transform<FPP,Cpu> is now a std::vector<Faust::MatGeneric<FPP,Cpu>*>, not std::vector<Faust::MatSparse<FPP,Cpu> >
	/*if(id>=size())
	  {
	  handleError(m_className,"get_fact : id exceed Faust::Transform size or id < 0");
	  }


	  return data[id];*/

}

template<typename FPP>
void Faust::Transform<FPP,Cpu>::get_fact(const faust_unsigned_int id,
		const int** rowptr,
		const int** col_ids,
		const FPP** elts,
		faust_unsigned_int* nnz,
		faust_unsigned_int* num_rows,
		faust_unsigned_int* num_cols) const
{
	MatGeneric<FPP,Cpu>* mat_ptr = get_fact(id, false);
	if(mat_ptr->getType() != MatType::Sparse)
		handleError(m_className, "get_fact(uint,uint**,uint**,FPP**,uint*,uint*,uint*): this prototype must be called only on sparse factors.");
	MatSparse<FPP,Cpu>* smat_ptr = dynamic_cast<MatSparse<FPP,Cpu>*>(mat_ptr);

	*rowptr = smat_ptr->getRowPtr();
	*col_ids = smat_ptr->getColInd();
	*elts = smat_ptr->getValuePtr();
	*nnz = smat_ptr->getNonZeros();
	*num_rows = smat_ptr->getNbRow();
	*num_cols = smat_ptr->getNbCol();
}

template<typename FPP>
void Faust::Transform<FPP,Cpu>::get_fact(const faust_unsigned_int id,
		int* d_outer_count_ptr, int* d_inner_ptr, FPP* d_elts,
		faust_unsigned_int* nnz,
		faust_unsigned_int* num_rows, faust_unsigned_int* num_cols,
		bool transpose /* default to false */) const
{
	const int* s_outer_count_ptr, *s_inner_ptr;
	const FPP* s_elts;
	this->get_fact(id, &s_outer_count_ptr, &s_inner_ptr, &s_elts,
			nnz, num_rows, num_cols);
	if(transpose)
	{
		MatSparse<FPP,Cpu> tmat(*nnz, *num_rows, *num_cols, s_elts, s_outer_count_ptr,
				s_inner_ptr, transpose);
		s_outer_count_ptr = tmat.getRowPtr();
		s_inner_ptr = tmat.getColInd();
		s_elts = tmat.getValuePtr();
		//do the copy here, otherwise we'll lose tmat and its buffers when out of scope
		memcpy(d_outer_count_ptr, s_outer_count_ptr, sizeof(int)*(*num_cols+1));
		memcpy(d_inner_ptr, s_inner_ptr, sizeof(int)**nnz);
		memcpy(d_elts, s_elts, sizeof(FPP)**nnz);
		// swap num_cols and num_rows
		// (with only these 2 variables -- F2 arithmetic trick)
		*num_cols = *num_cols^*num_rows;
		*num_rows = *num_cols^*num_rows;
		*num_cols = *num_cols^*num_rows;
	}
	else
	{
		memcpy(d_outer_count_ptr, s_outer_count_ptr, sizeof(int)*(*num_rows+1));
		memcpy(d_inner_ptr, s_inner_ptr, sizeof(int)**nnz);
		memcpy(d_elts, s_elts, sizeof(FPP)**nnz);
	}
}

template<typename FPP>
void Faust::Transform<FPP,Cpu>::get_fact(const faust_unsigned_int id,
		 const FPP ** elts,
		 faust_unsigned_int* num_rows,
		 faust_unsigned_int* num_cols) const
{
	MatGeneric<FPP,Cpu>* mat_ptr = get_fact(id, false);
	if(mat_ptr->getType() != MatType::Dense)
		handleError(m_className, "get_fact(uint,FPP**,uint*,uint*,uint*): this prototype must be called only on dense factors.");
	MatDense<FPP,Cpu>* dmat_ptr = dynamic_cast<MatDense<FPP,Cpu>*>(mat_ptr);
	*elts = dmat_ptr->getData();
	if(num_rows)
		*num_rows = dmat_ptr->getNbRow();
	if(num_cols)
		*num_cols = dmat_ptr->getNbCol();
}

template<typename FPP>
void Faust::Transform<FPP,Cpu>::get_fact(const faust_unsigned_int id,
		FPP* elts,
		faust_unsigned_int* num_rows,
		faust_unsigned_int* num_cols,
		const bool transpose /*=false*/) const
{
	const FPP* s_elts;
	get_fact(id, &s_elts, num_rows, num_cols);
	if(transpose)
	{
		// change the order of s_elts to transpose
		// the matrices are in Column-major order
		for(int j=0;j<*num_cols;j++)
			for(int i=0; i<*num_rows;i++)
				elts[i**num_cols+j] = s_elts[j**num_rows+i];
		// swap num_cols and num_rows
		// (with only these 2 variables -- F2 arithmetic trick)
		*num_cols = *num_cols^*num_rows;
		*num_rows = *num_cols^*num_rows;
		*num_cols = *num_cols^*num_rows;
	}
	else
	{
		memcpy(elts, s_elts, sizeof(FPP)*(*num_rows)*(*num_cols));
	}
}


template<typename FPP>
bool Faust::Transform<FPP,Cpu>::is_fact_sparse(const faust_unsigned_int id) const
{
	return get_fact(id, false)->getType() == MatType::Sparse;
}


template<typename FPP>
bool Faust::Transform<FPP,Cpu>::is_fact_dense(const faust_unsigned_int id) const
{
	return get_fact(id, false)->getType() == MatType::Dense;
}

template<typename FPP>
faust_unsigned_int Faust::Transform<FPP,Cpu>::get_fact_nnz(const faust_unsigned_int id) const
{
	return this->get_fact(id, false)->getNonZeros();
}

	template<typename FPP>
void Faust::Transform<FPP,Cpu>::push_back(const Faust::MatGeneric<FPP,Cpu>* M,const bool optimizedCopy /*default value = false */, const bool conjugate)
{
	if (size()>0)
	{
		if(data[size()-1]->getNbCol()!= M->getNbRow() || M->getNbRow()<1)
		{
			handleError(m_className,"push_back : incorrect dimensions");
		}
	}
	Faust::MatGeneric<FPP,Cpu>* M_copy = M->Clone(optimizedCopy);
	if(conjugate) M_copy->conjugate();
	data.push_back(M_copy);
	if(!dtor_delete_data) ref_man.acquire(M_copy);
	totalNonZeros += M_copy->getNonZeros();

#ifdef __COMPILE_TIMERS__
	this->t_multiply_vector.push_back(Faust::Timer());
#endif

}


	template<typename FPP>
void Faust::Transform<FPP,Cpu>::push_first(const Faust::MatGeneric<FPP,Cpu>* M, const bool optimizedCopy /*default value = false */)
{
	if (size()>0)
		if(this->getNbRow()!=M->getNbCol() || M->getNbRow()<1)
		{
			handleError(m_className,"push_first : incorrect dimensions");
		}

	Faust::MatGeneric<FPP,Cpu>* M_copy = M->Clone(optimizedCopy);
	data.insert(data.begin(),M_copy);
	if(!dtor_delete_data) ref_man.acquire(M_copy);
	totalNonZeros += M_copy->getNonZeros();

#ifdef __COMPILE_TIMERS__
	this->t_multiply_vector.insert(this->t_multiply_vector().begin(),Faust::Timer());
#endif
}

/*
   template<typename FPP>
   void Faust::Transform<FPP,Cpu>::pop_back(Faust::MatGeneric<FPP,Cpu>* M)
   {
   if (size()>0)
   {
   M = data[size()-1];
   data.pop_back();
   totalNonZeros -= M->getNonZeros();
#ifdef __COMPILE_TIMERS__
this->t_multiply_vector.pop_back();
#endif
}
handleWarning("Faust::Transform<FPP,Cpu>::pop_back : empty Faust::Transform");
}

template<typename FPP>
void Faust::Transform<FPP,Cpu>::pop_first(Faust::MatGeneric<FPP,Cpu>* M)
{

/// WARNING : not tested
if (size()>0)
{
M = data[0];
data.erase(data.begin());
totalNonZeros -= M->getNonZeros();
#ifdef __COMPILE_TIMERS__
this->t_multiply_vector.erase(this->t_multiply_vector.begin());
#endif
}
handleWarning("Faust::Transform<FPP,Cpu>::pop_back : empty Faust::Transform");
}

template<typename FPP>
void Faust::Transform<FPP,Cpu>::pop_first(Faust::MatGeneric<FPP,Cpu>* M) const
{
if (size()>0)
{
M = data[0]->Clone();

}

}*/

	template<typename FPP>
void Faust::Transform<FPP,Cpu>::transpose()
{
	int nbFact=size();
	reverse(data.begin(),data.end());
	for (int i=0;i<nbFact;i++)
	{
		data[i]->transpose();
	}
}

template<typename FPP>
void Faust::Transform<FPP,Cpu>::conjugate()
{
	typename vector<Faust::MatGeneric<FPP,Cpu>*>::iterator it;
	for(it = data.begin(); it != data.end(); it++){
		it->conjugate();
	}
}

template<typename FPP>
Faust::Vect<FPP,Cpu> Faust::Transform<FPP,Cpu>::multiply(const Faust::Vect<FPP,Cpu> x,const char opThis) const{
#ifdef __COMPILE_TIMERS__
	if (this->t_multiply_vector.size() != this->data.size())
		handleError(m_className,"multiply : invalid t_multiply_vector size to measure time");
#endif

	if (size() == 0)
		handleWarning("Faust::Transform<FPP,Cpu> : multiply : empty Faust::Transform<FPP,Cpu>");

	Faust::Vect<FPP,Cpu> vec(x);

	if (opThis == 'N')
	{
		for (int i=this->size()-1 ; i >= 0 ; i--)
		{
#ifdef __COMPILE_TIMERS__
			this->t_multiply_vector[i].start();
#endif

			data[i]->multiply(vec);

#ifdef __COMPILE_TIMERS__
			this->t_multiply_vector[i].stop();
#endif
		}
	}else
	{
		for (int i=0 ; i < this->size() ; i++)
		{
#ifdef __COMPILE_TIMERS__
			this->t_multiply_vector[i].start();
#endif

			data[i]->multiply(vec,opThis);

#ifdef __COMPILE_TIMERS__
			this->t_multiply_vector[i].stop();
#endif
		}

	}
	return vec;



}





template<typename FPP>
void Faust::Transform<FPP,Cpu>::get_nonortho_interior_prod_ids(int &start_id, int &end_id)
{
	start_id = -1;
	end_id = -1;
	int i = 0;
	while (i < this->size() && (data[i]->is_orthogonal() /*|| data[i]-this->is_identity*/))
		i++;
	if(i<data.size())
	{
		start_id = i;
		end_id = start_id;
		i = this->size()-1;
		while (i > start_id && (data[i]->is_orthogonal() /*|| data[i]-this->is_identity*/))
			i--;
		if(i > start_id)
			end_id = i;
	}
}

template<typename FPP>
Faust::MatSparse<FPP,Cpu> Faust::Transform<FPP,Cpu>::multiply(const Faust::MatSparse<FPP,Cpu> A,const char opThis) const
{
	if (size() == 0)
		handleWarning("Faust::Transform<FPP,Cpu> : multiply : empty Faust::Transform<FPP,Cpu>");

	Faust::MatSparse<FPP,Cpu> mat(A);

	if (opThis == 'N')
		for (int i=this->size()-1; i >= 0; i--)
			data[i]->multiply(mat,opThis);
	else
		for (int i=0; i < this->size(); i++)
			data[i]->multiply(mat,opThis);

	return mat;
}

template<typename FPP>
Faust::MatDense<FPP,Cpu> Faust::Transform<FPP,Cpu>::multiply(const Faust::MatDense<FPP,Cpu> A,const char opThis) const
{


	if (size() == 0)
		handleWarning("Faust::Transform<FPP,Cpu> : multiply : empty Faust::Transform<FPP,Cpu>");

	Faust::MatDense<FPP,Cpu> mat(A);


	if (opThis == 'N')
	{
		for (int i=this->size()-1; i >= 0; i--)
		{
			//#ifdef __COMPILE_TIMERS__
			//	this->t_multiply_mat[i].start();
			//#endif
			data[i]->multiply(mat,opThis);
				//#ifdef __COMPILE_TIMERS__
				//	this->t_multiply_mat[i].stop();
				//#endif
		}
	}else
	{
		for (int i=0; i < this->size(); i++)
		{
			//#ifdef __COMPILE_TIMERS__
			//	this->t_multiply_mat[i].start();
			//#endif
			data[i]->multiply(mat,opThis);
			//#ifdef __COMPILE_TIMERS__
			//	this->t_multiply_mat[i].start();
			//#endif
		}

	}

	return mat;

}

template<typename FPP>
Faust::MatDense<FPP,Cpu> Faust::Transform<FPP,Cpu>::multiply_par(const Faust::MatDense<FPP,Cpu> A, const char opThis) const
{
	int nth = std::thread::hardware_concurrency(); // https://en.cppreference.com/w/cpp/thread/thread/hardware_concurrency
	int barrier_count = nth;
	Faust::MatDense<FPP, Cpu> *M;
	std::vector<std::thread*> threads;
	std::vector<Faust::MatDense<FPP,Cpu>*> mats(nth);
	std::vector<Faust::MatGeneric<FPP, Cpu>*> data(this->data.size()+1);
	std::mutex barrier_mutex;
	std::condition_variable barrier_cv;
	int data_size = this->data.size()+1; // counting A in addition to Transform factors
	int num_per_th = data_size;
	num_per_th /= nth;
	num_per_th = num_per_th<2?2:num_per_th;
	// recompute effective needed nth
	nth = data_size / num_per_th;
	nth = nth%2?nth-1:nth; // nth must be even
	nth = nth == 0?1:nth;
//	cout << "nth=" << nth << endl;
	barrier_count = nth;
	int i = 0;
	if(opThis != 'N')
	{
		data[0] = const_cast<Faust::MatDense<FPP,Cpu>*>(&A);
		i = data_size-1;
	}
	else
		data[data_size-1] = const_cast<Faust::MatDense<FPP,Cpu>*>(&A);
	for(auto ptr: this->data)
	{
//		cout << "fac i=" << i << ": " << ptr << endl;
		if(opThis == 'N')
			data[i++] = ptr;
		else
			data[i--] = ptr;
	}
//	cout << "A ptr: " << &A << endl;
//	cout << A.norm() << endl;
	std::thread *th;
	for(i=0; i < nth; i++)
	{
//		th = new std::thread(&Faust::Transform<FPP,Cpu>::multiply_par_run, nth, i, num_per_th, data_size, opThis, data, mats, threads);
		th = new std::thread(Faust::multiply_par_run<FPP>, nth, i, num_per_th, data_size, opThis, std::ref(data), std::ref(mats), std::ref(threads), std::ref(barrier_mutex), std::ref(barrier_cv), std::ref(barrier_count));
		threads.push_back(th);
	}
	for(auto t: threads)
		if(t->get_id() != std::this_thread::get_id() )
			t->join();

	M = dynamic_cast<Faust::MatDense<FPP,Cpu>*>(mats[0]);
//	cout << "M:" << M << endl;
	MatDense<FPP,Cpu> M_;
	if(! M)
	{
		M_ = Faust::MatDense<FPP,Cpu>(*dynamic_cast<Faust::MatSparse<FPP,Cpu>*>(mats[0]));
		M = &M_;
	}
	//TODO: delete other thread mats
	//TODO: return a ptr instead of a copy
	//TODO: delete threads
	return *M;

}

template<typename FPP>
void Faust::multiply_par_run(int nth, int thid, int num_per_th, int data_size, char opThis, std::vector<Faust::MatGeneric<FPP, Cpu>*>& data, std::vector<Faust::MatDense<FPP,Cpu>*>& mats, std::vector<std::thread*>& threads, std::mutex & barrier_mutex, std::condition_variable & barrier_cv, int & barrier_count)
{
	Faust::MatSparse<FPP, Cpu> * sM;
	Faust::MatDense<FPP, Cpu> *M;
	Faust::MatDense<FPP,Cpu>* tmp = nullptr, *tmp2 = nullptr; // (data[end_id-1]);
//	std::cout << "num_per_th: " << num_per_th << std::endl;
	int start_nth = nth;
	while(num_per_th > 1)
	{
		if(thid < nth)
		{
			int first_id, end_id, id;
			first_id = num_per_th*thid;
			end_id = num_per_th*(thid+1);
			if(thid == start_nth-1 && end_id < data_size) 
				end_id = data_size;
			else
				end_id = end_id>data_size?data_size:end_id;
			id = 'N' == opThis? end_id-1: first_id;
//			cout << "thid=" << thid << " end orig adr.:" << data[id] << endl;
			if(tmp != nullptr) tmp2 = tmp; //keep track of previous tmp to delete it afterward
			if(sM = dynamic_cast<Faust::MatSparse<FPP,Cpu>*>(data[id]))
				tmp = new Faust::MatDense<FPP,Cpu>(*sM);
			else
				tmp = new Faust::MatDense<FPP,Cpu>(*(Faust::MatDense<FPP,Cpu>*)(data[id]));
			mats[thid] = tmp;
//			cout << "thid=" << thid << "mats[thid]=" << mats[thid] << "tmp=" << tmp << endl;
//			cout << "tmp.norm()" << tmp->norm() << endl;
			if(opThis == 'N')
				for(int i=end_id-2;i>=first_id; i--)
				{
//					cout << "mul:" << data[i] << " thid: " << thid << endl;
					if(sM = dynamic_cast<Faust::MatSparse<FPP,Cpu>*>(data[i]))
						sM->multiply(*mats[thid], opThis);
					else
						dynamic_cast<Faust::MatDense<FPP,Cpu>*>(data[i])->multiply(*mats[thid], opThis);
				}
			else
				for(int i=first_id+1;i < end_id; i++)
					if(sM = dynamic_cast<Faust::MatSparse<FPP,Cpu>*>(data[i]))
						sM->multiply(*mats[thid],opThis);
					else
						dynamic_cast<Faust::MatDense<FPP,Cpu>*>(data[i])->multiply(*mats[thid], opThis);
			//				mats[thid]->Display();
			//				cout << mats[thid]->norm() << endl;
//			cout << "thid=" << thid << "mats[thid]=" << mats[thid] << "data[thid]=" << data[thid] << endl;
			if(tmp2 != nullptr) delete tmp2;
		}
#define barrier() \
		{ \
			unique_lock<mutex> l1(barrier_mutex);  /* mutex block */ \
			barrier_count--;\
/*			string ite_str = i>=0?string(" (ite = ") + to_string(i) + ") ": "";*/ \
			if(barrier_count > 0)\
			{\
				/*	 the threads wait here for the last one to finish its iteration i*/ \
/*				std::cout << "thread " << thid << ite_str  << " is waiting on barrier (" << barrier_count << ")." << std::endl;*/ \
				barrier_cv.wait(l1/*, [&barrier_count, &num_threads]{return barrier_count == num_threads;}*/);\
/*				std::cout << "thread " << thid << ite_str <<" continues." << std::endl;*/\
			}\
			else {\
				/* the last thread to finish the iteration i wakes the other to continue the execution */ \
				/*				std::cout << "thread " << thid << ite_str << " reinits barrier." << std::endl;*/\
				barrier_count = start_nth;\
				barrier_cv.notify_all();\
			}\
		}
		barrier();
		if(thid < nth)
		{
			data[thid] = mats[thid];
		}
		barrier();
		if(nth > 1)
		{
			num_per_th = 2;
		}
		else
			num_per_th = 1;

		nth >>= 1;

		//		for(auto t: threads)
		//			if(t->get_id() != std::this_thread::get_id() )
		//				t->join();
//		cout <<"num_per_th: " <<  num_per_th << " thid:" << thid << endl;


	}
}


template<typename FPP>
Faust::MatDense<FPP,Cpu> Faust::Transform<FPP,Cpu>::multiply_omp(const Faust::MatDense<FPP,Cpu> A, const char opThis) const
{

	Faust::MatDense<FPP, Cpu> *M = nullptr;
#ifdef _OPENMP_
	int nth, thid, num_per_th, data_size;
	Faust::MatDense<FPP,Cpu>* mats[8];
	Faust::MatGeneric<FPP, Cpu> * data[this->data.size()+1];
	Faust::MatSparse<FPP, Cpu> * sM;
	Faust::MatDense<FPP,Cpu>* tmp; // (data[end_id-1]);
	int i = 0;
	for(auto ptr: this->data)
	{
//		cout << "i=" << i << " data[i]=" << ptr << endl;
		data[i++] = ptr;
	}
	data[i] = const_cast<Faust::MatDense<FPP,Cpu>*>(&A);
#pragma omp parallel private(nth, thid, tmp, num_per_th, data_size)
	{
		data_size = this->data.size()+1; // data + 1
		nth	= omp_get_num_threads();
		num_per_th = data_size;
		num_per_th = num_per_th%2?num_per_th-1:num_per_th;
		num_per_th /= nth;
		num_per_th = num_per_th<2?2:num_per_th;
//		std::cout << "omp th: " << nth << " num_per_th: " << num_per_th << std::endl;
		thid = omp_get_thread_num();
//		std::cout << "num_per_th: " << num_per_th << std::endl;
		while(num_per_th > 1)
		{
			int first_id, end_id, id;
			first_id = num_per_th*thid;
			end_id = num_per_th*(thid+1);
			if(first_id < data_size)
			{
				id = 'N' == opThis? end_id-1: first_id;
				if(sM = dynamic_cast<Faust::MatSparse<FPP,Cpu>*>(data[id]))
					tmp = new Faust::MatDense<FPP,Cpu>(*sM);
				else
					tmp = new Faust::MatDense<FPP,Cpu>(*(Faust::MatDense<FPP,Cpu>*)(data[id]));
				mats[thid] = tmp;
//				cout << "thid=" << thid << "mats[thid]=" << mats[thid] << "tmp=" << tmp << endl;
				if(opThis == 'N')
					for(int i=end_id-2;i>=first_id; i--)
					{
//						cout << "mul:" << data[i] << endl;
						if(sM = dynamic_cast<Faust::MatSparse<FPP,Cpu>*>(data[i]))
							sM->multiply(*mats[thid], opThis);
						else
							dynamic_cast<Faust::MatDense<FPP,Cpu>*>(data[i])->multiply(*mats[thid], opThis);
					}
				else
					for(int i=first_id+1;i < end_id; i++)
						if(sM = dynamic_cast<Faust::MatSparse<FPP,Cpu>*>(data[i]))
							sM->multiply(*mats[thid],opThis);
						else
							dynamic_cast<Faust::MatDense<FPP,Cpu>*>(data[i])->multiply(*mats[thid], opThis);
				data[thid] = mats[thid];
//				mats[thid]->Display();
//				cout << mats[thid]->norm() << endl;
//				cout << "thid=" << thid << "mats[thid]=" << mats[thid] << "data[thid]=" << data[thid] << endl;
				data_size = nth;
			}
			if(nth > 1)
				num_per_th = nth/(nth>>1);
			else
				num_per_th = 1;
			nth >>= 1;

#pragma omp barrier
		}
		//		if(thid == 0) M = dynamic_cast<Faust::MatDense<FPP,Cpu>*>(data[0]);
	}
	M = dynamic_cast<Faust::MatDense<FPP,Cpu>*>(mats[0]);
//	cout << M << endl;
	MatDense<FPP,Cpu> M_;
	if(! M)
	{
		M_ = Faust::MatDense<FPP,Cpu>(*dynamic_cast<Faust::MatSparse<FPP,Cpu>*>(mats[0]));
		M = &M_;
	}
	//TODO: delete other thread mats
	//TODO: return a ptr instead of a copy
#else
	throw runtime_error("It's not possible to call Faust::Transform::multiply_omp because the library hasn't been compiled with OpenMP enabled.");
#endif
	return *M;
}




template<typename FPP>
void Faust::Transform<FPP,Cpu>::setOp(const char op, faust_unsigned_int& nbRowOp, faust_unsigned_int& nbColOp)const
{
	if (size() > 0)
	{
		if(op == 'N')
		{
			nbRowOp=data[0]->getNbRow();
			nbColOp=data[size()-1]->getNbCol();
		}
		else if(op == 'T')
		{
			nbRowOp=data[size()-1]->getNbCol();
			nbColOp=data[0]->getNbRow();
		}
		else
			handleError(m_className,"setOp : invalid character");
	}else
	{
		nbRowOp = 0;
		nbColOp = 0;
		handleWarning("Faust::Transform<FPP,Cpu>::setOp : empty Faust::Transform");
	}
}

template<typename FPP>
void Faust::Transform<FPP,Cpu>::Display(const bool transpose /* default to false */,const bool displaying_small_mat_elts /*false by default*/)const
{
	std::cout << to_string(transpose,displaying_small_mat_elts);
}


template<typename FPP>
std::string Faust::Transform<FPP,Cpu>::to_string(const bool transpose /* default to false */, const bool displaying_small_mat_elts/* false by default */) const
{
	std::ostringstream str;

	if (size() == 0)
		str<<"empty Faust"<<std::endl;
	else
	{
		str<<"Faust size ";
		if(transpose)
			str << this->getNbCol() << "x" << this->getNbRow();
		else
			str << this->getNbRow()<<"x"<<this->getNbCol();
		str <<", density "<<1.0/getRCG()<< ", nnz_sum "<<this->get_total_nnz() << ", " << size() << " factor(s): "<< std::endl;
		int j;
		for (int i=0 ; i<size() ; i++)
		{
			if(transpose)
				j = size()-1-i;
			else
				j = i;
			str << "- FACTOR " << i;
			str << data[j]->to_string(transpose, displaying_small_mat_elts);
		}
	}
	return str.str();
}




#ifdef __COMPILE_TIMERS__
template<typename FPP>
void Faust::Transform<FPP,Cpu>::print_timers() const
{
	std::cout << "timers in Faust::Transform :" << endl;
	std::cout << "t_multiply_vector :" << endl;
	float sum_tps=0;
	for (int i=0;i<t_multiply_vector.size();i++)
		sum_tps+=t_multiply_vector[i].get_time();

	float current_time;
	float nb_call;
	faust_unsigned_int nnz;
	float density;
	for (int i=0;i<t_multiply_vector.size();i++)
	{
		current_time = t_multiply_vector[i].get_time();
		nb_call = t_multiply_vector[i].get_nb_call();
		nnz = data[i]->getNonZeros();
		density = 100.0 * data[i]->density(); //percentage 0 < value < 100

		std::cout<<"- FACTOR "<<i<<": TYPE ";
		if (data[i]->getType() == MatType::Dense)
			cout<<"dense ";
		else if (data[i]->getType() == Sparse)
			cout<<"sparse ";
		else
			cout<<"unknown ";
		cout<<", size "<<data[i]->getNbRow()<<"x"<<data[i]->getNbCol()<<" ,nnz "<<nnz<<", density % "<<density<<std::endl;
		std::cout<<" TPS "<<current_time<<" ,% "<<100*current_time/sum_tps<< " ,nb_calls  "<<nb_call<<std::endl;
		std::cout<<" 	, mean "<<current_time/((float) nb_call)<<std::endl;
	}
	std::cout<<"	total tps (for all factors)"<<sum_tps<<std::endl;



}
#endif








#endif

