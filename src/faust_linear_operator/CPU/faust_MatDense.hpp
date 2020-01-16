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
#ifndef __FAUST_MAT_HPP__
#define __FAUST_MAT_HPP__



#include "faust_linear_algebra.h"

//#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/SparseCore>
#include "faust_exception.h"
#include "faust_constant.h"
#include <cassert>
#ifdef __GEMM_WITH_OPENBLAS__
#include "faust_cblas_algebra.h"
#endif




template<typename FPP>
const char * Faust::MatDense<FPP,Cpu>::m_className = "Faust::MatDense<FPP,Cpu>::";


	template<typename FPP>
Faust::MatDense<FPP,Cpu>::MatDense(const FPP  *data_,const faust_unsigned_int nbRow, const faust_unsigned_int nbCol ) : Faust::MatGeneric<FPP,Cpu>(nbRow,nbCol), mat(nbRow,nbCol), isZeros(false)
{

#ifdef __COMPILE_TIMERS__
	t_constr.start();
#endif

	memcpy(getData(), data_, nbRow*nbCol*sizeof(FPP));


#ifdef __COMPILE_TIMERS__
	t_constr.stop();
#endif

}

template<typename FPP>
Faust::MatGeneric<FPP,Cpu>* Faust::MatDense<FPP,Cpu>::Clone(const bool isOptimize /*default value = false*/) const
{

	if (isOptimize)
	{
		Faust::MatSparse<FPP,Cpu> S((*this));
		return optimize((*this),S);
	}
	else
	{
		return new MatDense<FPP,Cpu>((*this));
	}
}


	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::resize(const faust_unsigned_int nbRow,const faust_unsigned_int nbCol)
{
#ifdef __COMPILE_TIMERS__
	t_resize.start();
#endif


	if ((this->dim1 != nbRow) || (this->dim2 != nbCol))
	{
		Faust::MatGeneric<FPP,Cpu>::resize(nbRow,nbCol);
		mat.resize(nbRow,nbCol);
	}

	isZeros = false;
this->is_identity = false;

#ifdef __COMPILE_TIMERS__
	t_resize.stop();
#endif
}


template<typename FPP>
faust_unsigned_int Faust::MatDense<FPP,Cpu>::getNonZeros()const
{	
	faust_unsigned_int nnz = 0;
	for (int i=0; i < this->getNbRow()*this->getNbCol(); i++)
	{
		if ( (*this)[i] != FPP(0.0) )
			nnz++; 
	}

	return nnz;	
}



	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::check_dim_validity()
{

#ifdef __COMPILE_TIMERS__
	t_check_dim.start();
#endif

	bool verifSize = (this->getNbCol() == mat.cols()) &&  (this->getNbRow() == mat.rows());

	if (!verifSize)
		handleError(m_className, "check_dim_validity : Size incompatibility in the Faust::MatDense");
#ifdef __COMPILE_TIMERS__
	t_check_dim.stop();
#endif

}


template<typename FPP>
Faust::MatDense<FPP,Cpu> Faust::MatDense<FPP,Cpu>::lower_tri(const bool diag) const
{
	Faust::MatDense<FPP,Cpu> tri = Faust::MatDense<FPP,Cpu>(this->dim1, this->dim2);
	if(diag)
		tri.mat = mat.template triangularView<Eigen::Lower>();
	else
		tri.mat = mat.template triangularView<Eigen::StrictlyLower>();
#ifdef DEBUG_TRI
	std::cout << "MatDense::lower_tri(" << diag << ")" <<std::endl;
	std::cout << "orig. mat.:" <<std::endl;
	std::cout << mat <<std::endl;
	std::cout << "tri. mat.:" <<std::endl;
	std::cout << tri.mat <<std::endl;
#endif
	return tri;
}

template<typename FPP>
Faust::MatDense<FPP,Cpu> Faust::MatDense<FPP,Cpu>::upper_tri(const bool diag) const
{
	Faust::MatDense<FPP,Cpu> tri = Faust::MatDense<FPP,Cpu>(this->dim1, this->dim2);
	if(diag)
		tri.mat = mat.template triangularView<Eigen::Upper>();
	else
		tri.mat = mat.template triangularView<Eigen::StrictlyUpper>();
#ifdef DEBUG_TRI
	std::cout << "MatDense::upper_tri(" << diag << ")" <<std::endl;
	std::cout << "orig. mat.:" <<std::endl;
	std::cout << mat <<std::endl;
	std::cout << "tri. mat.:" <<std::endl;
	std::cout << tri.mat <<std::endl;
#endif
	return tri;
}


template<typename FPP>
std::list<std::pair<int,int>> Faust::MatDense<FPP,Cpu>::nonzeros_indices() const
{
	std::list<std::pair<int,int>> nz_inds;
	if(this->is_identity)
#ifdef _MSC_VER
		for(int i=0;i<min(this->dim1, this->dim2);i++) // VS14 strange issue with std::min //C2589 or C2059
#else
		for(int i=0;i<std::min(this->dim1, this->dim2);i++)
#endif
			nz_inds.push_back(std::make_pair(i,i));
	else if(! isZeros)
	{
		int i,j;
		for(int k=0;k<this->dim1*this->dim2;k++)
			if(mat(k) != FPP(0))
			{
				j = k/this->dim1;
				i = k-j*this->dim1;
				nz_inds.push_back(std::make_pair(i,j));
			}
	}
	return nz_inds;
}


	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::setOnes()
{
	FPP* ptr_data = getData();
	for (int i=0 ; i<this->dim1*this->dim2; i++)
		ptr_data[i] = FPP(1.0);
	isZeros = false;
this->is_identity = false;
}

	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::setZeros()
{
	memset(getData(), 0, sizeof(FPP) * this->dim1*this->dim2);
	isZeros = true;
this->is_identity = false;
}

	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::setEyes()
{
	setZeros();
	FPP* ptr_data = getData();
#ifdef _MSC_VER
	for (int i=0;i<min(this->dim1,this->dim2); i++)
#else
	for (int i=0;i<std::min(this->dim1,this->dim2); i++)
#endif
		ptr_data[i*this->dim1+i] = FPP(1.0);
	if (this->dim1 == this->dim2)
	this->is_identity = true;
	isZeros = false;
}

template<typename FPP>
Faust::MatDense<FPP,Cpu> Faust::MatDense<FPP,Cpu>::eye(faust_unsigned_int nrows, faust_unsigned_int ncols)
{
	Faust::MatDense<FPP,Cpu> mat(nrows, ncols);
	mat.setEyes();
	return mat;
}

// EGALITE //
template<typename FPP>
bool Faust::MatDense<FPP,Cpu>::isEqual(const Faust::MatDense<FPP,Cpu> & B) const
{
	if ((this->getNbCol() != B.getNbCol()) || (this->getNbRow() != B.getNbRow()))
		handleError(m_className, "isEqual : dimension of the 2 matrix are not the same\n");

	if (isZeros)
		return B.isZeros;
	else if (B.isZeros)
		return isZeros;
	else
		return mat.isApprox(B.mat,FAUST_PRECISION);
}


template<typename FPP>
bool Faust::MatDense<FPP,Cpu>::isEqual(const Faust::MatDense<FPP,Cpu> & B, FPP threshold) const
{
	if ((this->getNbCol() != B.getNbCol()) || (this->getNbRow() != B.getNbRow()))
	{
		handleError(m_className, "isEqual : dimension of the 2 matrix are not the same\n");
	}
	bool egalite =true;
	for (int i=0;i<this->getNbRow();i++)
	{
		for (int j=0;j<this->getNbCol();j++)
		{
			if (std::abs(mat(i,j)==0))
			{
				if (	(std::abs(mat(i,j)-B.mat(i,j))) > threshold )
				{
					egalite = false;
					std::cout<<" i : "<<i<<" j : "<<j<<std::endl;
				}
			}else
			{
				if (	(std::abs(mat(i,j)-B.mat(i,j))/std::abs(mat(i,j))) > threshold )
				{
					egalite = false;
					std::cout<<" i : "<<i<<" j : "<<j<<std::endl;
				}
			}
		}
	}

	return egalite;
}


// OPERATION BASIQUE //
	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::transpose()
{

#ifdef __COMPILE_TIMERS__
	t_transpose.start();
#endif

	if(isZeros || this->is_identity)
	{
		resize(this->dim2,this->dim1);
#ifdef __COMPILE_TIMERS__
		t_transpose.stop();
#endif
		return;
	}

	mat = mat.transpose().eval();
	faust_unsigned_int dim1_copy = this->dim1;
	this->dim1 = this->dim2;
	this->dim2 = dim1_copy;

#ifdef __COMPILE_TIMERS__
	t_transpose.stop();
#endif

}
	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::conjugate()
{

#ifdef __COMPILE_TIMERS__
	//t_conjugate.start(); //TODO
#endif

	if(isZeros)
	{
#ifdef __COMPILE_TIMERS__
		//			t_conjugate.stop(); //TODO
#endif
		return;
	}
	conjugate(true);
#ifdef __COMPILE_TIMERS__
	//t_conjugate.stop(); //TODO
#endif

}
	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::conjugate(const bool eval)
{

#ifdef __COMPILE_TIMERS__
	//t_conjugate.start(); //TODO
#endif

	if(isZeros)
	{
#ifdef __COMPILE_TIMERS__
		//			t_conjugate.stop(); //TODO
#endif
		return;
	}

	if(eval)
		mat = mat.conjugate().eval();
	else
		mat = mat.conjugate();
#ifdef __COMPILE_TIMERS__
	//t_conjugate.stop(); //TODO
#endif

}

template<typename FPP>
void Faust::MatDense<FPP,Cpu>::real()
{
	mat = mat.real().eval().template cast<FPP>();
}

template<typename FPP>
	void Faust::MatDense<FPP,Cpu>::multiply(Faust::MatSparse<FPP,Cpu> & M, char opThis) const
{
	//TODO: this function should not rely on spgemm which forces us to transconjugate
	bool alt_conj = false;
	if(opThis == 'H')
		//compute (this, M^H)^H
		opThis = 'N';
	else if(opThis == 'T')
	{
		//compute (this->conjugate(), M^H)^H
		const_cast<Faust::MatDense<FPP,Cpu>*>(this)->conjugate();
		opThis = 'N';
		alt_conj = true; //keep track of conj to undo it on the end // dirty
	}
	else
		//compute (this^H, M^H)^H
		opThis = 'H';
	Faust::MatDense<FPP,Cpu> out;
	Faust::spgemm<FPP>(M, *this,out,1.0,0.0,'H', opThis);
	M = out;
	M.makeCompression();
	M.transpose();
	M.conjugate();
	if(alt_conj) const_cast<Faust::MatDense<FPP,Cpu>*>(this)->conjugate();
}

	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::multiplyRight(Faust::MatDense<FPP,Cpu> const& A)
{

#ifdef __COMPILE_TIMERS__
	t_mult_right.start();
#endif

	if (this->dim2 != A.dim1)
	{
		handleError(m_className, "multiplyRight : dimension conflict between matrix");
	}

	if(A.is_id())
	{
#ifdef __COMPILE_TIMERS__
		t_mult_right.stop();
#endif
		return;
	}

	if(isZeros || A.isZeros)
	{
		//std::cout<<"zero"<<std::endl;
		resize(this->dim1,A.dim2);
		FPP *const ptr_data_dst = getData();
		memset(ptr_data_dst, 0, sizeof(FPP) * this->dim1*this->dim2);
		isZeros = true;
	this->is_identity = false;
#ifdef __COMPILE_TIMERS__
		t_mult_right.stop();
#endif
		return;
	}

	if(this->is_identity)
	{
		this->operator=(A);
#ifdef __COMPILE_TIMERS__
		t_mult_right.stop();
#endif
		return;
	}

	MatDense this_copy((*this));
	gemm<FPP>(this_copy,A,(*this),1.0,0.0,'N','N');


#ifdef __COMPILE_TIMERS__
	t_mult_right.stop();
#endif

}

	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::multiplyRight(Faust::MatSparse<FPP,Cpu> const& A)
{

#ifdef __COMPILE_TIMERS__
	t_mult_right.start();
#endif

	if (this->dim2 != A.dim1)
	{
		handleError(m_className, "multiplyRight : dimension conflict between matrix");
	}

	if(isZeros)
	{
		//std::cout<<"zero"<<std::endl;
		resize(this->dim1,A.dim2);
		FPP *const ptr_data_dst = getData();
		memset(ptr_data_dst, 0, sizeof(FPP) * this->dim1*this->dim2);
		isZeros = true;
	this->is_identity = false;
#ifdef __COMPILE_TIMERS__
		t_mult_right.stop();
#endif
		return;
	}

	if(this->is_identity)
	{
		this->operator=(A);
#ifdef __COMPILE_TIMERS__
		t_mult_right.stop();
#endif
		return;
	}

	MatDense this_copy((*this));
	spgemm<FPP>(this_copy,A,(*this),1.0,0.0,'N','N');


#ifdef __COMPILE_TIMERS__
	t_mult_right.stop();
#endif

}

/*
   template<typename FPP>
   void Faust::MatDense<FPP,Cpu>::multiplyLeft(Faust::MatDense<FPP,Cpu> const& A)
   {

#ifdef __COMPILE_TIMERS__
t_mult_left.start();
#endif

if (this->dim1 != A.dim2)
{
handleError(m_className, "multiplyLeft : dimension conflict between matrix");
}

if(A.isIdentity)
{
#ifdef __COMPILE_TIMERS__
t_mult_left.stop();
#endif
return;
}

if(isZeros || A.isZeros)
{
resize(A.dim1,this->dim2);
FPP *const ptr_data_dst = getData();
memset(ptr_data_dst, 0, sizeof(FPP) * this->dim1*this->dim2);
isZeros = true;
isIdentity = false;
#ifdef __COMPILE_TIMERS__
t_mult_left.stop();
#endif
return;
}

if(this->is_identity)
{
this->operator=(A);
#ifdef __COMPILE_TIMERS__
t_mult_left.stop();
#endif
return;
}


MatDense this_copy((*this));
gemm<FPP>(A,this_copy,(*this),1.0,0.0,'N','N');

#ifdef __COMPILE_TIMERS__
t_mult_left.stop();
#endif

}
*/



template<typename FPP>
Real<FPP> Faust::MatDense<FPP,Cpu>::spectralNorm(const faust_unsigned_int nbr_iter_max,FPP threshold, int & flag,Faust::BlasHandle<Cpu> blas_handle/*=Faust::BlasHandle<Cpu>()*/) const
{
#ifdef __COMPILE_FPPIMERS__
	t_spectral_norm2.start();
#endif
	if(isZeros)
	{
		flag = -2;
#ifdef __COMPILE_TIMERS__
		t_spectral_norm2.stop();
#endif
		return 0;
	}

	if(this->is_identity)
	{
		flag = -3;
#ifdef __COMPILE_TIMERS__
		t_spectral_norm2.stop();
#endif
		return 1;
	}

	faust_unsigned_int nb_row = this->getNbRow();
	faust_unsigned_int nb_col = this->getNbCol();


	Faust::MatDense<FPP,Cpu> AtA;
	if (nb_row <= nb_col)
	{
		gemm((*this),(*this),AtA,(FPP)1.,(FPP)0.0,'N', 'H'/*'T'*/);
	}else
	{
		gemm((*this),(*this),AtA,(FPP)1.,(FPP)0.0,'H' /*equiv. to 'T' if FPP real*/,'N');
	}

	Real<FPP> res=std::sqrt(std::abs(power_iteration(AtA,nbr_iter_max,threshold,flag)));


#ifdef __COMPILE_TIMERS__
	t_spectral_norm2.stop();
#endif
	return res;

}







	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::scalarMultiply(FPP const lambda)
{
#ifdef __COMPILE_TIMERS__
	t_scalar_multiply.start();
#endif
	mat = lambda * mat;
#ifdef __COMPILE_TIMERS__
	t_scalar_multiply.stop();
#endif
}


	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::add(Faust::MatDense<FPP,Cpu> const& A)
{
#ifdef __COMPILE_TIMERS__
	t_add.start();
#endif
	if ((this->getNbCol() != A.getNbCol()) || (this->getNbRow() != A.getNbRow()))
	{
		handleError(m_className, "add : matrix dimension not equal");
	}
	mat = mat + A.mat;
	isZeros = false;
	this->is_identity = false;
#ifdef __COMPILE_TIMERS__
	t_add.stop();
#endif
}

	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::sub(Faust::MatDense<FPP,Cpu> const& A)
{
#ifdef __COMPILE_TIMERS__
	t_sub.start();
#endif
	if ((this->getNbCol() != A.getNbCol()) || (this->getNbRow() != A.getNbRow()))
	{
		std::cout<<"sub"<<std::endl;			
		std::cout<<" this dimension ("<<this->getNbRow()<<","<<this->getNbCol()<<")"<<std::endl;		
		std::cout<<" A dimension ("<<A.getNbRow()<<","<<A.getNbCol()<<")"<<std::endl;	
		handleError(m_className, "sub : matrix dimension not equal");
	}
	mat = mat - A.mat;

	isZeros = false;
this->is_identity = false;

#ifdef __COMPILE_TIMERS__
	t_sub.stop();
#endif
}

	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::sub(Faust::MatSparse<FPP,Cpu> const& A)
{
#ifdef __COMPILE_TIMERS__
	t_sub.start();
#endif
	if ((this->getNbCol() != A.getNbCol()) || (this->getNbRow() != A.getNbRow()))
	{
		std::cout<<"sub"<<std::endl;			
		std::cout<<" this dimension ("<<this->getNbRow()<<","<<this->getNbCol()<<")"<<std::endl;		
		std::cout<<" A dimension ("<<A.getNbRow()<<","<<A.getNbCol()<<")"<<std::endl;	
		handleError(m_className, "sub : matrix dimension not equal");
	}
	mat = mat - A.mat;

	isZeros = false;
this->is_identity = false;

#ifdef __COMPILE_TIMERS__
	t_sub.stop();
#endif
}



template<typename FPP>
std::string Faust::MatDense<FPP,Cpu>::to_string(const bool transpose /* set to false by default */, const bool displaying_small_mat_elts /* false by default */) const
{
	//using ostringstream because it's easier for concatenation (of complex and any number)
	std::ostringstream  str;

	str << " (" << (typeid(*getData()) == typeid(std::complex<double>) || typeid(*getData()) == typeid(std::complex<float>)?"complex":"real") << ")";
	str<<" DENSE, "; //number of trailing spaces matters to align with SPARSE to_string()
	str << Faust::MatGeneric<FPP,Cpu>::to_string(transpose);
	if(isZeros)
		str <<"zeros matrix flag" <<std::endl;
	else
	{
		if (displaying_small_mat_elts && this->dim1*this->dim2 < 100)
		{
			for (int i=0 ; i<this->dim1 ; i++)
			{
				for(int j=0 ; j<this->dim2 ; j++)
					str << (*this)(i,j) << " " ;
				str << std::endl;
			}
		}
	}
	return str.str();
}



// Affichage
template<typename FPP>
void Faust::MatDense<FPP,Cpu>::Display() const
{
	std::cout<<this->to_string();
}

template<typename FPP>
void Faust::MatDense<FPP,Cpu>::print_file(const char* filename)const
{
	std::ofstream fichier;
	fichier.open(filename);
	for (int i=0 ; i<this->getNbRow() ; i++)
	{
		for (int j=0 ; j<this->getNbCol() ; j++)
			fichier << std::setprecision(20) <<mat(i,j) << " ";
		fichier << std::endl;
	}
	fichier.close();
}


/// SURCHARGE OPERATEUR ///
// affectation
	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::operator=(Faust::MatDense<FPP,Cpu> const& A)
{
	mat = A.mat;
	this->dim1 = A.dim1;
	this->dim2 = A.dim2;
	isZeros = A.isZeros;
	this->is_identity = A.is_identity;
	this->is_ortho = A.is_ortho;
}

template<typename FPP>
	template<typename FPP1>
void Faust::MatDense<FPP,Cpu>::operator=(Faust::MatDense<FPP1,Cpu> const& A)
{
	resize(A.dim1,A.dim2);
	for (int i=0;i<this->dim1*this->dim2;i++)
		(*this)[i]=(FPP) A(i);
	isZeros = A.isZeros;
	this->is_identity = A.is_identity;
	this->is_ortho = A.is_ortho;
}

	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::operator=(Faust::MatSparse<FPP,Cpu> const& S)
{
	const_cast<Faust::MatSparse<FPP,Cpu>&>(S).update_dim();
	S.check_dim_validity();
	try {
	resize(S.getNbRow(),S.getNbCol());
	setZeros();
	FPP*const ptr_data = getData();



	for(int i=0 ; i< S.mat.outerSize() ; i++)
	{
		for(typename Eigen::SparseMatrix<FPP,Eigen::RowMajor>::InnerIterator it(S.mat,i); it; ++it)
		{
			ptr_data[it.col() * this->dim1 + it.row()] = it.value();

		}
	}		
	}
	catch(std::bad_alloc e){
		std::cerr << "Out of memory." << std::endl;
	}
	isZeros = false;
this->is_identity = false;
	this->is_ortho = S.is_ortho;
}


	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::operator*=(const Faust::MatSparse<FPP,Cpu>& S)
{
	if(this->dim2 != S.dim1)
	{
		handleError(m_className, "operator*= : incorrect matrix dimensions");
	}

	if (this->is_identity)
	{
		this->operator=(S);
	this->is_identity = false;
		isZeros = false;
	}
	else if (isZeros)
	{
		resize(this->dim1, S.dim2);
		setZeros();
	}
	else
	{
		mat = mat * S.mat;
		this->dim2 = S.dim2;
	}

}


	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::operator+=(const Faust::MatSparse<FPP,Cpu>& S)
{
	if(this->dim1!=S.getNbRow() || this->dim2!=S.getNbCol())
	{
		handleError(m_className,"operator+= : incorrect matrix dimensions");
	}
	mat += S.mat;
	this->is_identity = false;
	isZeros = false;
}



	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::scalarMultiply(Faust::MatDense<FPP,Cpu> const& A)
{
	if(this->dim1!=A.dim1 || this->dim2!=A.dim2)
	{
		handleError(m_className,"scalarMultiply : incorrect matrix dimensions\n");
	}
	mat = (mat.array() * A.mat.array()).matrix();
this->is_identity = false;
	isZeros = false;
}


	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::multiplyLeft(const Faust::MatSparse<FPP,Cpu>& S,const char TransS)
{
	faust_unsigned_int nbColOpS,nbRowOpS;	
	S.setOp(TransS,nbRowOpS,nbColOpS);	
	if(nbColOpS != this->dim1)
	{
		//std::cerr << "Error in Faust::MatDense<FPP,Cpu>::operator*= : incorrect matrix dimensions" << std::endl;
		//exit(EXIT_FAILURE);
		handleError(m_className,"multiplyLeft : incorrect matrix dimensions\n");
	}

	if (this->is_identity)
	{
		this->operator=(S);
	this->is_identity = false;
		isZeros = false;

		if (TransS == 'T')
			this->transpose();
	}
	else if (isZeros)
	{
		resize(nbRowOpS, this->dim2);
		setZeros();
	}
	else
	{

		if (TransS == 'N') 			
			mat = S.mat * mat;
		else
			mat = S.mat.transpose() * mat;

		this->dim1 = nbRowOpS;
	}

}

	template<typename FPP>
void Faust::MatDense<FPP,Cpu>::init_from_file(const char* filename)
{
#ifdef __COMPILE_TIMERS__
	t_print_file.start();
#endif
	// la premiere ligne contient 2 entiers : dim1 et dim2
	// chacune des autres lignes contient une valeur par ligne
	// suivant la premiere dimension puis la deuxieme

	std::ifstream* vec_stream;
	vec_stream = new std::ifstream(filename);
	if (!vec_stream->is_open())
		handleError(m_className, "init_from_file : unable to open file");
	std::istream_iterator<FPP> start(*vec_stream), eos;
	std::vector<FPP> vec(start, eos);

	if((vec[0]*vec[1]+2) != vec.size())
	{
		handleError(m_className, "init_from_file : problem with the file");
	}
	resize(vec[0],vec[1]);
	memcpy(getData(), &vec[2], sizeof(FPP) * this->dim1 * this->dim2);

	isZeros = false;
this->is_identity = false;

#ifdef __COMPILE_TIMERS__
	t_print_file.stop();
#endif
}


template<typename FPP>
matvar_t* Faust::MatDense<FPP, Cpu>::toMatIOVar(bool transpose, bool conjugate) const
{
	matvar_t *var = NULL;
	size_t dims[2];
	int opt = typeid(mat(0,0))==typeid(std::complex<double>(1.0,1.0))?MAT_F_COMPLEX:0;
	mat_complex_split_t z = {NULL,NULL};
	//
	if(transpose)
	{
		dims[0] = this->getNbCol();
		dims[1] = this->getNbRow();
		if(!opt){
			Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic> mat_copy(mat.transpose().eval());
			//		mat_copy.transposeInPlace(); //undoing the transposition
			var = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, (FPP*) mat_copy.data() /*mat.transpose().eval().data() //  doesn't work so we copy above */, opt);
		}
		else {
			Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> dst_re(mat.rows(), mat.cols());
			Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> dst_im(mat.rows(), mat.cols());
			dst_re = mat.real()/*.transpose()*/.template cast<double>();
			if(conjugate)
				dst_im = (-(mat.imag()))/*.transpose()*/.template cast<double>();
			else
				dst_im = mat.imag()/*.transpose()*/.template cast<double>();
			dst_re.transposeInPlace();
			dst_im.transposeInPlace();
			z.Re = dst_re.data();
			z.Im = dst_im.data();
			var = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &z /*mat.transpose().eval().data() //  doesn't work so we copy above */, opt);
		}
	}
	else
	{
		dims[0] = this->getNbRow();
		dims[1] = this->getNbCol();
		if(!opt) // we use directly the data pointer (col-major order organized)
			// but without the MAT_F_DONT_COPY_DATA flag, MatIO copies the data internally
			var = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, (FPP*) mat.data(), opt);
		else {
			Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> dst_re(mat.rows(), mat.cols());
			Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> dst_im(mat.rows(), mat.cols());
			dst_re = mat.real().template cast<double>();
			if(conjugate)
				dst_im = (-(mat.imag())).template cast<double>();
			else
				dst_im = mat.imag().template cast<double>();
			z.Re = dst_re.data();
			z.Im = dst_im.data();
			var = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &z /*mat.transpose().eval().data() //  doesn't work so we copy above */, opt);
		}
	}
	return var;
}

template<typename FPP>
Real<FPP> Faust::MatDense<FPP, Cpu>::normL1(faust_unsigned_int& col_id, const bool transpose /* default false */) const
{
	//TODO: refactor this function with MatSparse::normL1()
	faust_unsigned_int i, j, max_j;
	Real<FPP> sum, max_sum;
	Eigen::Matrix<FPP, Eigen::Dynamic,Eigen::Dynamic> vec;
	int dim1, dim2;
	if(transpose){
		dim1 = this->getNbRow();
		dim2 = this->getNbCol();
	}
	{
		dim1 = this->getNbCol();
		dim2 = this->getNbRow();
	}
	for(j=0;j<dim1;j++)
	{
		if(transpose) //get a row
			vec=mat.block(j,0,1,this->getNbCol());
		else //get a column
			vec=mat.block(0,j,this->getNbRow(),1);
		for(i=0,sum=0;i<dim2;i++)
			sum += std::abs(vec.data()[i]);
		if(j==0 || std::abs(sum) > std::abs(max_sum))
		{
			max_sum = sum;
			max_j = j;
		}
	}
	col_id = max_j;
	return max_sum;
}

template<typename FPP>
Real<FPP> Faust::MatDense<FPP, Cpu>::normL1(const bool transpose /* default false */) const
{
	faust_unsigned_int id; //useless but mandatory
	return normL1(id, transpose);
}

template<typename FPP>
Faust::Vect<FPP,Cpu> Faust::MatDense<FPP,Cpu>::get_col(faust_unsigned_int id) const
{
	if(id > this->getNbCol())
		handleError("Faust::MatDense", "Too big column index passed to get_col().");
	Eigen::Matrix<FPP, Eigen::Dynamic,1> vec;
	vec = mat.col(id);
	return Vect<FPP,Cpu>(this->getNbRow(),vec.data());
}

template<typename FPP>
Faust::Vect<FPP,Cpu> Faust::MatDense<FPP,Cpu>::get_row(faust_unsigned_int id) const
{
	if(id > this->getNbRow())
		handleError("Faust::MatDense", "Too big row index passed to get_col().");
	Eigen::Matrix<FPP, Eigen::Dynamic,1> vec;
	vec = mat.row(id);
	return Vect<FPP,Cpu>(this->getNbCol(),vec.data());
}

template<typename FPP>
void Faust::MatDense<FPP,Cpu>::delete_col(int id)
{
	if(id < 0 || id >= this->getNbCol()) throw std::out_of_range(std::string(m_className)+"::delete_col() index out of bounds");
	memcpy(getData()+this->getNbRow()*id, getData()+this->getNbRow()*(id+1), this->getNbRow()*(this->getNbCol()-id-1)*sizeof(FPP));
	this->dim2--;
}

template<typename FPP>
void Faust::MatDense<FPP,Cpu>::delete_row(int id)
{
	if(id < 0 || id >= this->getNbRow()) throw std::out_of_range(std::string(m_className)+"::delete_row() index out of bounds");
	Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic> mat(this->getNbRow()-1, this->getNbCol());
	mat.topRows(id) = this->mat.topRows(id);
	mat.bottomRows(this->getNbRow()-id-1) = this->mat.bottomRows(this->getNbRow()-id-1);
	this->mat = mat;
	this->dim1--;
}

template<typename FPP>
Faust::MatDense<FPP,Cpu>* Faust::MatDense<FPP,Cpu>::get_cols(faust_unsigned_int start_col_id, faust_unsigned_int num_cols) const
{
	//TODO: check args
	FPP *data = new FPP[this->getNbRow()*num_cols];
	memcpy(data, getData()+start_col_id*this->getNbRow(), num_cols*this->getNbRow()*sizeof(FPP));
	Faust::MatDense<FPP, Cpu>* cols = new Faust::MatDense<FPP, Cpu>(data, this->getNbRow(), num_cols);
	delete[] data;
	return cols;
}

template<typename FPP>
Faust::MatDense<FPP,Cpu>* Faust::MatDense<FPP,Cpu>::get_cols(faust_unsigned_int* col_ids, faust_unsigned_int n) const
{
	//TODO: check args
	FPP *data = new FPP[this->getNbRow()*n];
	for(faust_unsigned_int i=0; i<n;i++)
		memcpy(data+i*this->getNbRow(), getData()+col_ids[i]*this->getNbRow(), this->getNbRow()*sizeof(FPP));
	Faust::MatDense<FPP, Cpu>* cols = new Faust::MatDense<FPP, Cpu>(data, this->getNbRow(), n);
	delete[] data;
	return cols;
}

template<typename FPP>
Faust::MatDense<FPP,Cpu>* Faust::MatDense<FPP,Cpu>::get_cols(std::vector<int> col_ids) const
{
	//TODO: check args
	int n = col_ids.size();
	FPP *data = new FPP[this->getNbRow()*n];
	for(int i=0; i<n;i++)
		memcpy(data+i*this->getNbRow(), getData()+col_ids[i]*this->getNbRow(), this->getNbRow()*sizeof(FPP));
	Faust::MatDense<FPP, Cpu>* cols = new Faust::MatDense<FPP, Cpu>(data, this->getNbRow(), n);
	delete[] data;
	return cols;
}

template<typename FPP>
Faust::MatDense<FPP,Cpu>* Faust::MatDense<FPP,Cpu>::get_rows(faust_unsigned_int start_row_id, faust_unsigned_int num_rows) const
{
	//TODO: check args
	FPP *data = new FPP[this->getNbCol()*num_rows];
	for(faust_unsigned_int i = 0; i < this->getNbCol(); i++){
		memcpy(data+i*num_rows, getData()+start_row_id+i*this->getNbRow(), num_rows*sizeof(FPP));
	}
	Faust::MatDense<FPP, Cpu>* rows = new Faust::MatDense<FPP, Cpu>(data, num_rows, this->getNbCol());
	delete[] data;
	return rows;
}

template<typename FPP>
Faust::MatDense<FPP,Cpu>* Faust::MatDense<FPP,Cpu>::get_rows(faust_unsigned_int* row_ids, faust_unsigned_int n) const
{
	//TODO: check args
	FPP *data = new FPP[this->getNbCol()*n];
	for(faust_unsigned_int i = 0; i < n; i++) // row_ids[i]
		for(faust_unsigned_int j = 0; j < this->getNbCol(); j++)
			//copy ele row_ids[i],j
			data[j*n+i] = this->mat(row_ids[i],j);
	Faust::MatDense<FPP, Cpu>* rows = new Faust::MatDense<FPP, Cpu>(data, n, this->getNbCol());
	delete data;
	return rows;
}

template<typename FPP>
FPP Faust::MatDense<FPP,Cpu>::min_coeff() const
{
	return mat.minCoeff();
}

template<typename FPP>
Faust::Vect<FPP,Cpu> Faust::MatDense<FPP,Cpu>::rowwise_min() const
{
	return Faust::Vect<FPP,Cpu>(this->getNbCol(), mat.rowwise().minCoeff().eval().data());
}

template<typename FPP>
Faust::Vect<FPP,Cpu> Faust::MatDense<FPP,Cpu>::rowwise_min(int* col_indices) const
{
	Faust::Vect<FPP,Cpu> vec(this->getNbRow());
	for(int i=0;i<this->getNbRow();i++)
		vec.getData()[i] = mat.row(i).minCoeff(col_indices+i);
	return vec;
}

template<typename FPP>
Faust::Vect<FPP,Cpu> Faust::MatDense<FPP,Cpu>::rowwise_max(int* col_indices) const
{
	Faust::Vect<FPP,Cpu> vec(this->getNbRow());
	for(int i=0;i<this->getNbRow();i++)
	{
//		FPP max = FPP(std::numeric_limits<double>::min()) /* contrary to gcc and clang Visual Studio is not happy with that line: error C2589: '(': illegal token on right side of '::' */, e;
		FPP max, e;
		for(int j=0;j<this->getNbCol();j++)
		{
			e = getData()[j*this->getNbRow()+i];
			if(j == 0 || Faust::fabs(e) > Faust::fabs(max))
			{
				max = e;
				if(col_indices != nullptr) col_indices[i] = j;
			}
		}
		vec.getData()[i] = max;
	}
	return vec;
}

template<typename FPP>
void Faust::MatDense<FPP,Cpu>::abs()
{
	// not able to use cwiseAbs() from Eigen for certain versions of the lib and compilers
	// so doing it by hand, but shouldn't be slower
	faust_unsigned_int col_offset, i, j;
	for(j=0;j<this->getNbCol();j++)
	{
		col_offset = this->getNbRow()*j;
		for(i=0;i<this->getNbRow();i++)
			getData()[col_offset+i] = FPP(Faust::fabs((*this)(i,j)));
	}
}

	template<typename FPP>
Faust::MatDense<FPP, Cpu>* Faust::MatDense<FPP, Cpu>::randMat(faust_unsigned_int num_rows, faust_unsigned_int num_cols)
{
	Faust::MatDense<FPP, Cpu>* mat = NULL;
	try {
		Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic> mat_tmp;
		mat_tmp = Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>::Random(num_rows, num_cols);
		mat = new Faust::MatDense<FPP, Cpu>(mat_tmp.data(), num_rows, num_cols);
	}
	catch(std::bad_alloc e) {
//		std::cerr << "Out of memory." << std::endl;
	}
	return mat;
}

	template<typename FPP>
Faust::MatDense<FPP, Cpu>* Faust::MatDense<FPP, Cpu>::randMat(faust_unsigned_int num_rows, faust_unsigned_int num_cols, float density)
{
	Faust::MatDense<FPP, Cpu>* mat = NULL;
	try {
		//TODO: refactor with above randMat()
		Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic> mat_tmp;
		std::default_random_engine generator(time(NULL));
		std::uniform_real_distribution<double> distribution(0, 1);
		mat_tmp = Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>::Random(num_rows, num_cols);
		for(int i=0;i<num_rows;i++)
			for(int j=0;j<num_cols;j++)
				if(distribution(generator) > density)
					mat_tmp(i,j) = 0;
		mat = new Faust::MatDense<FPP, Cpu>(mat_tmp.data(), num_rows, num_cols);
		mat_tmp.resize(0,0);
	}
	catch(std::bad_alloc e) {
//		std::cerr << "Out of memory." << std::endl;
	}
	return mat;
}

	template<typename FPP>
Faust::MatDense<FPP, Cpu>* Faust::MatDense<FPP, Cpu>::randMat(faust_unsigned_int num_rows, faust_unsigned_int num_cols, float density, bool per_row)
{
	Faust::MatDense<FPP, Cpu>* mat = nullptr;
	try {
		// create a rand MatSparse and convert it to MatDense (because the density accuracy of MatSparse::randMat() is better)
		MatSparse<FPP, Cpu>* spmat = MatSparse<FPP, Cpu>::randMat(num_rows, num_cols, density, per_row);
		if(spmat != nullptr)
		{
			mat = new MatDense<FPP, Cpu>(*spmat);
			delete spmat;
		}
	}
	catch(std::bad_alloc e) {
//		std::cerr << "Out of memory." << std::endl;
	}
	return mat;
}

template<typename FPP>
Faust::Vect<FPP,Cpu> Faust::MatDense<FPP, Cpu>::diagonal(int index)
{
	FPP data[this->getNbRow()-index];
	for(int i=0;i < this->getNbRow()-index; i++)
		data[i] = *(this->getData()+i+(this->getNbRow()*(i+index)));
	Faust::Vect<FPP,Cpu> diag(this->getNbRow()-index, data);
	return diag;
}

#ifdef __COMPILE_TIMERS__
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_constr;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_get_coeff;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_get_coeffs;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_set_coeff;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_set_coeffs;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_set_coeffs2;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_resize;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_check_dim;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_max;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_transpose;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_mult_right;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_mult_left;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_scalar_multiply;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_add;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_sub;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_print_file;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_multiply;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_gemm;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_add_ext;

template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_spectral_norm;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_spectral_norm2;
template<typename FPP>
Faust::Timer Faust::MatDense<FPP,Cpu>::t_power_iteration;

template<typename FPP>
void Faust::MatDense<FPP,Cpu>::print_timers()const
{
	std::cout << "timers in Faust::MatDense :" << std::endl;
	std::cout << "t_constr          = " << t_constr.get_time()          << " s for "<< t_constr.get_nb_call()          << " calls" << std::endl;
	std::cout << "t_get_coeff       = " << t_get_coeff.get_time()       << " s for "<< t_get_coeff.get_nb_call()       << " calls" << std::endl;
	std::cout << "t_get_coeffs      = " << t_get_coeffs.get_time()      << " s for "<< t_get_coeffs.get_nb_call()      << " calls" << std::endl;
	std::cout << "t_set_coeff       = " << t_set_coeff.get_time()       << " s for "<< t_set_coeff.get_nb_call()       << " calls" << std::endl;
	std::cout << "t_set_coeffs      = " << t_set_coeffs.get_time()      << " s for "<< t_set_coeffs.get_nb_call()      << " calls" << std::endl;
	std::cout << "t_set_coeffs2     = " << t_set_coeffs2.get_time()     << " s for "<< t_set_coeffs2.get_nb_call()     << " calls" << std::endl;
	std::cout << "t_resize          = " << t_resize.get_time()          << " s for "<< t_resize.get_nb_call()          << " calls" << std::endl;
	std::cout << "t_check_dim       = " << t_check_dim.get_time()       << " s for "<< t_check_dim.get_nb_call()       << " calls" << std::endl;
	std::cout << "t_max             = " << t_max.get_time()             << " s for "<< t_max.get_nb_call()             << " calls" << std::endl;
	std::cout << "t_transpose       = " << t_transpose.get_time()       << " s for "<< t_transpose.get_nb_call()       << " calls" << std::endl;
	std::cout << "t_mult_right      = " << t_mult_right.get_time()      << " s for "<< t_mult_right.get_nb_call()      << " calls" << std::endl;
	std::cout << "t_mult_left       = " << t_mult_left.get_time()       << " s for "<< t_mult_left.get_nb_call()       << " calls" << std::endl;
	std::cout << "t_scalar_multiply = " << t_scalar_multiply.get_time() << " s for "<< t_scalar_multiply.get_nb_call() << " calls" << std::endl;
	std::cout << "t_add             = " << t_add.get_time()             << " s for "<< t_add.get_nb_call()             << " calls" << std::endl;
	std::cout << "t_sub             = " << t_sub.get_time()             << " s for "<< t_sub.get_nb_call()             << " calls" << std::endl;
	std::cout << "t_print_file      = " << t_print_file.get_time()      << " s for "<< t_print_file.get_nb_call()      << " calls" << std::endl<<std::endl;


	std::cout << "timers in Faust::MatDense / LinearAlgebra :" << std::endl;
	std::cout << "t_multiply        = " << t_multiply.get_time()        << " s for "<< t_multiply.get_nb_call()        << " calls" << std::endl;
	std::cout << "t_gemm            = " << t_gemm.get_time()            << " s for "<< t_gemm.get_nb_call()            << " calls" << std::endl;
	std::cout << "t_add_ext         = " << t_add_ext.get_time()         << " s for "<< t_add_ext.get_nb_call()         << " calls" << std::endl<<std::endl<<std::endl;
}
#endif

#endif
