/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2021):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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



#include <type_traits>
#include "faust_linear_algebra.h"
#include <limits>

//#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/SparseCore>
#include "faust_exception.h"
#include "faust_constant.h"
#include <cassert>
#include <Eigen/SVD>
#include "faust_init_from_matio.h"

namespace Faust
{

	template<typename FPP>
		const char * MatDense<FPP,Cpu>::m_className = "MatDense<FPP,Cpu>::";


	template<typename FPP>
		MatDense<FPP,Cpu>::MatDense(const FPP  *data_,const faust_unsigned_int nbRow, const faust_unsigned_int nbCol ) : MatGeneric<FPP,Cpu>(nbRow,nbCol), mat(nbRow,nbCol), isZeros(false)
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
		MatDense<FPP,Cpu>::MatDense(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, const FPP  *data_): MatDense<FPP,Cpu>(data_, nbRow, nbCol) {}

	template<typename FPP>
		MatGeneric<FPP,Cpu>* MatDense<FPP,Cpu>::Clone(const bool isOptimize /*default value = false*/) const
		{

			if (isOptimize)
			{
				MatSparse<FPP,Cpu> S((*this));
				return optimize((*this),S);
			}
			else
			{
				return new MatDense<FPP,Cpu>((*this));
			}
		}


	template<typename FPP>
		void MatDense<FPP,Cpu>::resize(const faust_unsigned_int nbRow,const faust_unsigned_int nbCol)
		{
#ifdef __COMPILE_TIMERS__
			t_resize.start();
#endif


			if ((this->dim1 != nbRow) || (this->dim2 != nbCol))
			{
				MatGeneric<FPP,Cpu>::resize(nbRow,nbCol);
				mat.resize(nbRow,nbCol);
			}

			isZeros = false;
			this->is_identity = false;

#ifdef __COMPILE_TIMERS__
			t_resize.stop();
#endif
		}


	template<typename FPP>
		faust_unsigned_int MatDense<FPP,Cpu>::getNonZeros()const
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
	size_t MatDense<FPP,Cpu>::getNBytes() const
	{
		return this->getNbCol()*this->getNbRow()*sizeof(FPP);
	}

	template<typename FPP>
		void MatDense<FPP,Cpu>::check_dim_validity()
		{

#ifdef __COMPILE_TIMERS__
			t_check_dim.start();
#endif

			bool verifSize = (this->getNbCol() == mat.cols()) &&  (this->getNbRow() == mat.rows());

			if (!verifSize)
				handleError(m_className, "check_dim_validity : Size incompatibility in the MatDense");
#ifdef __COMPILE_TIMERS__
			t_check_dim.stop();
#endif

		}


	template<typename FPP>
		MatDense<FPP,Cpu> MatDense<FPP,Cpu>::lower_tri(const bool diag) const
		{
			MatDense<FPP,Cpu> tri = MatDense<FPP,Cpu>(this->dim1, this->dim2);
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
		MatDense<FPP,Cpu> MatDense<FPP,Cpu>::upper_tri(const bool diag) const
		{
			MatDense<FPP,Cpu> tri = MatDense<FPP,Cpu>(this->dim1, this->dim2);
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
		std::list<std::pair<int,int>> MatDense<FPP,Cpu>::nonzeros_indices() const
		{
			std::list<std::pair<int,int>> nz_inds;
			if(this->is_identity)
#ifdef _MSC_VER
				for(int i=0;i<min(this->dim1, this->dim2);i++) // VS14 strange issue with std::min //C2589 or C2059 // TODO/ fix rather with -DNOMINMAX
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
		void MatDense<FPP,Cpu>::setOnes()
		{
			FPP* ptr_data = getData();
			for (int i=0 ; i<this->dim1*this->dim2; i++)
				ptr_data[i] = FPP(1.0);
			isZeros = false;
			this->is_identity = false;
		}

	template<typename FPP>
		void MatDense<FPP,Cpu>::setZeros()
		{
			memset(getData(), 0, sizeof(FPP) * this->dim1*this->dim2);
			isZeros = true;
			this->is_identity = false;
		}

	template<typename FPP>
		void MatDense<FPP,Cpu>::setEyes()
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
		void MatDense<FPP, Cpu>::setNZtoOne()
		{
			auto eps = std::numeric_limits<Real<FPP>>::epsilon();
			for(int i=0;i<this->dim1*this->dim2;i++)
			{
				mat.data()[i] = static_cast<FPP>(std::abs((*this)(i)) > eps?1:0);
			}
		}

	template<typename FPP>
		void MatDense<FPP,Cpu>::setData(const FPP* data, int32_t nrows, int32_t ncols)
		{
			if(nrows != this->dim1 || ncols != this->dim2)
				resize(nrows, ncols);
			memcpy(getData(), data, sizeof(FPP)*nrows*ncols);
		}

	template<typename FPP>
		MatDense<FPP,Cpu> MatDense<FPP,Cpu>::eye(faust_unsigned_int nrows, faust_unsigned_int ncols)
		{
			MatDense<FPP,Cpu> mat(nrows, ncols);
			mat.setEyes();
			return mat;
		}

	// EGALITE //
	template<typename FPP>
		bool MatDense<FPP,Cpu>::isEqual(const MatDense<FPP,Cpu> & B) const
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
		bool MatDense<FPP,Cpu>::isEqual(const MatDense<FPP,Cpu> & B, FPP threshold) const
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
		void MatDense<FPP,Cpu>::transpose()
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

			mat = mat.transpose().eval(); // equiv. to transposeInPlace() https://eigen.tuxfamily.org/dox/classEigen_1_1DenseBase.html#ac501bd942994af7a95d95bee7a16ad2a
			faust_unsigned_int dim1_copy = this->dim1;
			this->dim1 = this->dim2;
			this->dim2 = dim1_copy;

#ifdef __COMPILE_TIMERS__
			t_transpose.stop();
#endif

		}
	template<typename FPP>
		void MatDense<FPP,Cpu>::conjugate()
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
		void MatDense<FPP,Cpu>::conjugate(const bool eval)
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
		void MatDense<FPP,Cpu>::adjoint()
		{
			if(isZeros)
			{
				this->resize(this->dim2, this->dim1);
				return;
			}

			mat.adjointInPlace();
			faust_unsigned_int dim1_copy = this->dim1;
			this->dim1 = this->dim2;
			this->dim2 = dim1_copy;
		}

	template<typename FPP>
		void MatDense<FPP,Cpu>::real()
		{
			mat = mat.real().eval().template cast<FPP>();
		}

	template<typename FPP>
		void MatDense<FPP,Cpu>::real(MatDense<Real<FPP>, Cpu> & real_mat) const
		{
			real_mat.resize(this->getNbRow(), this->getNbCol());
			real_mat.mat = mat.real().eval().template cast<Real<FPP>>();
		}

	template<typename FPP>
		template<typename FPP2>
		MatDense<FPP2, Cpu> MatDense<FPP, Cpu>::to_real() const
		{
			Faust::MatDense<FPP2, Cpu> ddmat;
			ddmat.resize(this->getNbRow(), this->getNbCol());
			ddmat.mat = mat.real().eval().template cast<Real<FPP2>>();
			return ddmat;
		}

	template<typename FPP>
		void MatDense<FPP,Cpu>::multiply(MatSparse<FPP,Cpu> & M, char opThis) const
		{
			//TODO: this function should not rely on spgemm which forces us to transconjugate
			bool alt_conj = false;
			if(opThis == 'H')
				//compute (this, M^H)^H
				opThis = 'N';
			else if(opThis == 'T')
			{
				//compute (this->conjugate(), M^H)^H
				const_cast<MatDense<FPP,Cpu>*>(this)->conjugate();
				opThis = 'N';
				alt_conj = true; //keep track of conj to undo it on the end // dirty
			}
			else
				//compute (this^H, M^H)^H
				opThis = 'H';
			MatDense<FPP,Cpu> out;
			spgemm<FPP>(M, *this,out,1.0,0.0,'H', opThis);
			M = out;
			M.makeCompression();
			M.transpose();
			M.conjugate();
			if(alt_conj) const_cast<MatDense<FPP,Cpu>*>(this)->conjugate();
		}

	template<typename FPP>
		void MatDense<FPP,Cpu>::multiplyRight(MatDense<FPP,Cpu> const& A)
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
		void MatDense<FPP,Cpu>::multiplyRight(MatSparse<FPP,Cpu> const& A)
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
	   void MatDense<FPP,Cpu>::multiplyLeft(MatDense<FPP,Cpu> const& A)
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
Real<FPP> MatDense<FPP,Cpu>::spectralNorm(const faust_unsigned_int nbr_iter_max,FPP threshold, int & flag,BlasHandle<Cpu> blas_handle/*=BlasHandle<Cpu>()*/) const
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


	MatDense<FPP,Cpu> AtA;
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
void MatDense<FPP,Cpu>::scalarMultiply(FPP const lambda)
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
void MatDense<FPP,Cpu>::add(MatDense<FPP,Cpu> const& A)
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
void MatDense<FPP,Cpu>::add(MatSparse<FPP,Cpu> const& A)
{
#ifdef __COMPILE_TIMERS__
	t_add.start();
#endif
	if ((this->getNbCol() != A.getNbCol()) || (this->getNbRow() != A.getNbRow()))
	{
		handleError(m_className, "add : matrix dimension not equal");
	}
	mat += A.mat;
	isZeros = false;
	this->is_identity = false;
#ifdef __COMPILE_TIMERS__
	t_add.stop();
#endif
}

	template<typename FPP>
void MatDense<FPP,Cpu>::sub(MatDense<FPP,Cpu> const& A)
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
void MatDense<FPP,Cpu>::sub(MatSparse<FPP,Cpu> const& A)
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
std::string MatDense<FPP,Cpu>::to_string(const bool transpose /* set to false by default */, const bool displaying_small_mat_elts /* false by default */) const
{
	//using ostringstream because it's easier for concatenation (of complex and any number)
	std::ostringstream  str;
	str << MatGeneric<FPP,Cpu>::to_string(Dense, transpose);
	if(isZeros)
		str <<"zeros matrix flag" <<std::endl;
	else
	{
		if (displaying_small_mat_elts && this->dim1*this->dim2 < 1000)
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

template<typename FPP>
std::string Faust::MatDense<FPP,Cpu>::to_string(int32_t nrows, int32_t ncols, bool transpose, Real<FPP> density, int32_t nnz, bool is_identity)
{
	return Faust::MatGeneric<FPP,Cpu>::to_string(nrows, ncols, transpose, density, nnz, is_identity, Dense);
}


// Affichage
template<typename FPP>
void MatDense<FPP,Cpu>::Display() const
{
	std::cout<<this->to_string();
}

template<typename FPP>
void MatDense<FPP,Cpu>::print_file(const char* filename)const
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
void MatDense<FPP,Cpu>::operator=(MatDense<FPP,Cpu> const& A)
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
void MatDense<FPP,Cpu>::operator=(MatDense<FPP1,Cpu> const& A)
{
	resize(A.dim1,A.dim2);
	for (int i=0;i<this->dim1*this->dim2;i++)
		(*this)[i]=(FPP) A(i);
	isZeros = A.isZeros;
	this->is_identity = A.is_identity;
	this->is_ortho = A.is_ortho;
}

	template<typename FPP>
void MatDense<FPP,Cpu>::operator=(MatSparse<FPP,Cpu> const& S)
{
	const_cast<MatSparse<FPP,Cpu>&>(S).update_dim();
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
void MatDense<FPP,Cpu>::operator*=(const MatSparse<FPP,Cpu>& S)
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
void MatDense<FPP,Cpu>::operator+=(const MatSparse<FPP,Cpu>& S)
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
void MatDense<FPP,Cpu>::nonzerosToOnes() //TODO: rename nonZerosToOnes
{
	auto a = (Eigen::abs(mat.array()) > 0 || Eigen::isnan(mat.array()) || Eigen::isinf(mat.array())).select(mat, Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>::Ones(this->dim1, this->dim2));
	this->is_identity = false;
	isZeros = false;
}

	template<typename FPP>
void MatDense<FPP,Cpu>::scalarMultiply(MatDense<FPP,Cpu> const& A)
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
void MatDense<FPP,Cpu>::multiplyLeft(const MatSparse<FPP,Cpu>& S,const char TransS)
{
	faust_unsigned_int nbColOpS,nbRowOpS;	
	S.setOp(TransS,nbRowOpS,nbColOpS);	
	if(nbColOpS != this->dim1)
	{
		//std::cerr << "Error in MatDense<FPP,Cpu>::operator*= : incorrect matrix dimensions" << std::endl;
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
void MatDense<FPP,Cpu>::init_from_file(const char* filename)
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
void MatDense<FPP, Cpu>::from_matio_var(matvar_t *var)
{
	matio_classes matio_class;
	if(std::is_same<FPP, float>::value)
		matio_class = MAT_C_SINGLE;
	else
		matio_class = MAT_C_DOUBLE;

	if(var->class_type != matio_class
			|| var->rank != 2
			|| var->data_size != sizeof(Real<FPP>))
		handleError("MatDense::from_matio_var", "variable isn't of good type.");

	resize(var->dims[0], var->dims[1]);

	if(std::is_same<FPP,Real<FPP>>::value)
		memcpy(getData(), var->data, sizeof(FPP)*this->getNbRow()*this->getNbCol());
	else
	{
		// (it should be) complex data
		mat_complex_split_t c_data = *static_cast<mat_complex_split_t*>(var->data);
		for(int i=0;i<var->dims[0]*var->dims[1];i++)
		{
			((Real<FPP>*)this->getData())[i*2] = ((Real<FPP>*)c_data.Re)[i];
			((Real<FPP>*)this->getData())[i*2+1] = ((Real<FPP>*)c_data.Im)[i];
		}
	}

}

template<typename FPP>
void MatDense<FPP, Cpu>::read_from_mat_file(const char *filepath, const char *variable_name)
{
	matvar_t* matvar = faust_matio_read_variable(filepath, variable_name);
	from_matio_var(matvar);
}

template<typename FPP>
void MatDense<FPP, Cpu>::save_to_mat_file(const char *filepath, const char *var_name)
{
	int ret;
	matvar_t* matvar = toMatIOVar(false, false, var_name);
	mat_t* matfp = Mat_CreateVer(filepath, NULL, MAT_FT_MAT5);
	if(matfp == NULL)
		handleError("Faust::MatDense::save_to_mat_file()", "Failed creating file");
//	Mat_VarPrint(matvar, 1);
	ret = Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE); //TODO: enable compression ?
	if(ret)
		handleError("Faust::MatDense::save_to_mat_file", (std::string("Failed writing the MatDense to a Matlab file error code: ")+std::to_string(ret)).c_str());
	Mat_VarFree(matvar);
	Mat_Close(matfp);
}

template<typename FPP>
matvar_t* MatDense<FPP, Cpu>::toMatIOVar(bool transpose, bool conjugate, const char *var_name/*=nullptr*/) const
{
	matvar_t *var = nullptr;
	size_t dims[2];
	int opt = typeid(mat(0,0))==typeid(std::complex<Real<FPP>>(1.0,1.0))?MAT_F_COMPLEX:0;
	mat_complex_split_t z = {nullptr,nullptr};
	matio_types matio_type;
	matio_classes matio_class;
	if(std::is_same<FPP, float>::value)
	{
		matio_type = MAT_T_SINGLE;
		matio_class = MAT_C_SINGLE;
	}
	else
	{
		matio_type = MAT_T_DOUBLE;
		matio_class = MAT_C_DOUBLE;
	}
	//
	if(transpose)
	{
		dims[0] = this->getNbCol();
		dims[1] = this->getNbRow();
		if(!opt){
			Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic> mat_copy(mat.transpose().eval());
			//		mat_copy.transposeInPlace(); //undoing the transposition
			var = Mat_VarCreate(var_name, matio_class, matio_type, 2, dims, (FPP*) mat_copy.data() /*mat.transpose().eval().data() //  doesn't work so we copy above */, opt);
		}
		else {
			Eigen::Matrix<Real<FPP>,Eigen::Dynamic,Eigen::Dynamic> dst_re(mat.rows(), mat.cols());
			Eigen::Matrix<Real<FPP>,Eigen::Dynamic,Eigen::Dynamic> dst_im(mat.rows(), mat.cols());
			dst_re = mat.real()/*.transpose()*/.template cast<Real<FPP>>();
			if(conjugate)
				dst_im = (-(mat.imag()))/*.transpose()*/.template cast<Real<FPP>>();
			else
				dst_im = mat.imag()/*.transpose()*/.template cast<Real<FPP>>();
			dst_re.transposeInPlace();
			dst_im.transposeInPlace();
			z.Re = dst_re.data();
			z.Im = dst_im.data();
			var = Mat_VarCreate(var_name, matio_class, matio_type, 2, dims, &z /*mat.transpose().eval().data() //  doesn't work so we copy above */, opt);
		}
	}
	else
	{
		dims[0] = this->getNbRow();
		dims[1] = this->getNbCol();
		if(!opt) // we use directly the data pointer (col-major order organized)
			// but without the MAT_F_DONT_COPY_DATA flag, MatIO copies the data internally
			var = Mat_VarCreate(var_name, matio_class, matio_type, 2, dims, (FPP*) mat.data(), opt);
		else {
			Eigen::Matrix<Real<FPP>,Eigen::Dynamic,Eigen::Dynamic> dst_re(mat.rows(), mat.cols());
			Eigen::Matrix<Real<FPP>,Eigen::Dynamic,Eigen::Dynamic> dst_im(mat.rows(), mat.cols());
			dst_re = mat.real().template cast<Real<FPP>>();
			if(conjugate)
				dst_im = (-(mat.imag())).template cast<Real<FPP>>();
			else
				dst_im = mat.imag().template cast<Real<FPP>>();
			z.Re = dst_re.data();
			z.Im = dst_im.data();
			var = Mat_VarCreate(var_name, matio_class, matio_type, 2, dims, &z /*mat.transpose().eval().data() //  doesn't work so we copy above */, opt);
		}
	}
	return var;
}

template<typename FPP>
Real<FPP> MatDense<FPP, Cpu>::normL1(faust_unsigned_int& id, const bool transpose /* default false */) const
{
	Real<FPP> sum, max_sum;
	if(transpose)
		max_sum = mat.cwiseAbs().rowwise().sum().maxCoeff();
	else
		max_sum = mat.cwiseAbs().colwise().sum().maxCoeff();
	return max_sum;
}

template<typename FPP>
Real<FPP> MatDense<FPP, Cpu>::normL1(const bool transpose /* default false */) const
{
	faust_unsigned_int id; //useless but mandatory
	return normL1(id, transpose);
}

template<typename FPP>
Real<FPP> Faust::MatDense<FPP, Cpu>::normInf(const bool transpose/*=false*/) const
{
	return normL1(!transpose);
}

template<typename FPP>
Real<FPP> Faust::MatDense<FPP, Cpu>::normInf(faust_unsigned_int& row_id, const bool transpose/*=false*/) const
{
	return normL1(row_id, !transpose);
}

	template<typename FPP>
void MatDense<FPP,Cpu>::normalize(int norm_type/*=-2*/)
{
	Real<FPP> n;
	switch(norm_type)
	{
		case -2: // frobenius
			n = norm();
			break;
		case 2:
			int flag; // not used
			n = spectralNorm(FAUST_NORM2_MAX_ITER, FAUST_PRECISION, flag);
		case 1:
			n = normL1();
			break;
		case -1:
			n = normInf();
			break;
		default:
			throw std::runtime_error("Unknown kind of norm asked for normalization.");
	}
	if(n != Real<FPP>(0))
		scalarMultiply(FPP(1.0/n));
	else
		throw std::domain_error("the norm is zero, can't normalize.");
}

template<typename FPP>
bool MatDense<FPP,Cpu>::containsNaN() const
{
	return mat.hasNaN();
}

	template<typename FPP>
void MatDense<FPP,Cpu>::copyBuf(FPP* dst_buf) const
{
	memcpy(dst_buf, getData(), sizeof(FPP)*this->getNbCol()*this->getNbRow());
}


template<typename FPP>
Vect<FPP,Cpu> MatDense<FPP,Cpu>::get_col(faust_unsigned_int id) const
{
	if(id > this->getNbCol())
		handleError("MatDense", "Too big column index passed to get_col().");
	Eigen::Matrix<FPP, Eigen::Dynamic,1> vec;
	vec = mat.col(id);
	return Vect<FPP,Cpu>(this->getNbRow(),vec.data());
}

template<typename FPP>
Vect<FPP,Cpu> MatDense<FPP,Cpu>::get_row(faust_unsigned_int id) const
{
	if(id > this->getNbRow())
		handleError("MatDense", "Too big row index passed to get_col().");
	Eigen::Matrix<FPP, Eigen::Dynamic,1> vec;
	vec = mat.row(id);
	return Vect<FPP,Cpu>(this->getNbCol(),vec.data());
}

template<typename FPP>
void MatDense<FPP,Cpu>::swap_cols(const faust_unsigned_int id1, const faust_unsigned_int id2)
{
	Vect<FPP,Cpu> v = this->get_col(id1);
	memcpy(this->getData()+this->getNbRow()*id1, this->getData()+this->getNbRow()*id2, sizeof(FPP)*this->getNbRow());
	memcpy(this->getData()+this->getNbRow()*id2, v.getData(), sizeof(FPP)*this->getNbRow());
}

template<typename FPP>
void MatDense<FPP,Cpu>::swap_rows(const faust_unsigned_int id1, const faust_unsigned_int id2)
{
	faust_unsigned_int ncols = this->getNbCol();
	MatDense<FPP,Cpu> swap_row(1, ncols);
	swap_row.mat.row(0) = mat.row(id1);
	mat.row(id1) = mat.row(id2);
	mat.row(id2) = swap_row.mat.row(0);
}

	template<typename FPP>
void MatDense<FPP,Cpu>::delete_col(int id)
{
	if(id < 0 || id >= this->getNbCol()) throw std::out_of_range(std::string(m_className)+"::delete_col() index out of bounds");
	memcpy(getData()+this->getNbRow()*id, getData()+this->getNbRow()*(id+1), this->getNbRow()*(this->getNbCol()-id-1)*sizeof(FPP));
	this->dim2--;
}

	template<typename FPP>
void MatDense<FPP,Cpu>::delete_row(int id)
{
	if(id < 0 || id >= this->getNbRow()) throw std::out_of_range(std::string(m_className)+"::delete_row() index out of bounds");
	Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic> mat(this->getNbRow()-1, this->getNbCol());
	mat.topRows(id) = this->mat.topRows(id);
	mat.bottomRows(this->getNbRow()-id-1) = this->mat.bottomRows(this->getNbRow()-id-1);
	this->mat = mat;
	this->dim1--;
}

template<typename FPP>
MatDense<FPP,Cpu>* MatDense<FPP,Cpu>::get_cols(faust_unsigned_int start_col_id, faust_unsigned_int num_cols) const
{
	//TODO: check args
	FPP *data = new FPP[this->getNbRow()*num_cols];
	memcpy(data, getData()+start_col_id*this->getNbRow(), num_cols*this->getNbRow()*sizeof(FPP));
	MatDense<FPP, Cpu>* cols = new MatDense<FPP, Cpu>(data, this->getNbRow(), num_cols);
	delete[] data;
	return cols;
}

template<typename FPP>
MatDense<FPP,Cpu>* MatDense<FPP,Cpu>::get_cols(const faust_unsigned_int* col_ids, faust_unsigned_int n) const
{
	FPP *data = new FPP[this->getNbRow()*n];
	faust_unsigned_int j;
	for(faust_unsigned_int i=0; i<n;i++)
	{
		j = col_ids[i];
		if(j >= this->getNbCol())
		{
			delete[] data;
			throw std::runtime_error("Index out of column range.");
		}
		memcpy(data+i*this->getNbRow(), getData()+col_ids[i]*this->getNbRow(), this->getNbRow()*sizeof(FPP));
	}
	MatDense<FPP, Cpu>* cols = new MatDense<FPP, Cpu>(data, this->getNbRow(), n);
	delete[] data;
	return cols;
}

template<typename FPP>
MatDense<FPP,Cpu>* MatDense<FPP,Cpu>::get_cols(std::vector<int> col_ids) const
{
	int n = col_ids.size();
	MatDense<FPP, Cpu>* cols = new MatDense<FPP, Cpu>(this->getNbRow(), n);
	FPP *data = cols->getData();
	int i = 0;
	for(auto j: col_ids)
	{
		if(j >= this->getNbCol() || j < 0)
		{
			throw std::runtime_error("Index out of column range.");
		}
		memcpy(data+i*this->getNbRow(), getData()+j*this->getNbRow(), this->getNbRow()*sizeof(FPP));
		i++;
	}
	return cols;
}

template<typename FPP>
MatDense<FPP,Cpu>* MatDense<FPP,Cpu>::get_rows(faust_unsigned_int start_row_id, faust_unsigned_int num_rows) const
{
	MatDense<FPP, Cpu>* rows = new MatDense<FPP, Cpu>(num_rows, this->getNbCol());
	get_rows(start_row_id, num_rows, *rows);
	return rows;
}

template<typename FPP>
void MatDense<FPP,Cpu>::get_rows(faust_unsigned_int start_row_id, faust_unsigned_int num_rows, MatDense<FPP, Cpu>& out_rows) const
{
	if(start_row_id >= this->getNbRow() || start_row_id+num_rows > this->getNbRow())
		throw std::domain_error("get_rows: arguments out of row indices.");
	out_rows.resize(num_rows, this->getNbCol());
	FPP *data = out_rows.getData();
	for(faust_unsigned_int i = 0; i < this->getNbCol(); i++)
		memcpy(data+i*num_rows, getData()+start_row_id+i*this->getNbRow(), num_rows*sizeof(FPP));
}

template<typename FPP>
MatDense<FPP,Cpu>* MatDense<FPP,Cpu>::get_rows(const faust_unsigned_int* row_ids, faust_unsigned_int n) const
{
	MatDense<FPP, Cpu>* rows = new MatDense<FPP, Cpu>(n, this->getNbCol());
	this->get_rows(row_ids, n, *rows);
	return rows;
}

template<typename FPP>
void MatDense<FPP,Cpu>::get_rows(const faust_unsigned_int* row_ids, faust_unsigned_int n, MatDense<FPP, Cpu>& out_rows) const
{
	out_rows.resize(n, this->getNbCol());
	FPP *data = out_rows.getData();
	for(faust_unsigned_int i = 0; i < n; i++) // row_ids[i]
		for(faust_unsigned_int j = 0; j < this->getNbCol(); j++)
		{
			if(row_ids[i] >= this->getNbRow())
				throw std::runtime_error("Index out of row range.");
			//copy ele row_ids[i],j
			data[j*n+i] = this->mat(row_ids[i],j);
		}
}

	template<typename FPP>
bool MatDense<FPP,Cpu>::eq_cols(const MatDense<FPP, Cpu> & other, faust_unsigned_int id_this, faust_unsigned_int id_other, const Real<FPP>& precision) const
{
	if(this->getNbRow() != other.getNbRow()) return false;
	for(int i=0;i<this->getNbRow();i++)
	{
		if(std::abs((*this)(i, id_this)-other(i, id_other)) > precision)
			return false;
	}
	return true;
}

	template<typename FPP>
bool MatDense<FPP,Cpu>::eq_rows(const MatDense<FPP, Cpu> & other, faust_unsigned_int id_this, faust_unsigned_int id_other, const Real<FPP>& precision) const
{
	if(this->getNbCol() != other.getNbCol()) return false;
	for(int j=0;j<this->getNbCol();j++)
	{
		if(std::abs((*this)(id_this,j)-other(id_other,j)) > precision)
			return false;
	}
	return true;
}

template<typename FPP>
template<typename SVDImpl>
void MatDense<FPP, Cpu>::best_low_rank(const int &r, MatDense<FPP,Cpu> &bestX, MatDense<FPP, Cpu> &bestY, SVDImpl svd) const
{
	if(bestX.getNbRow() != this->getNbRow() || r != bestX.getNbCol())
		bestX.resize(this->getNbRow(), r);
	if(bestY.getNbRow() != this->getNbRow() || r != bestY.getNbCol())
		bestY.resize(r, this->getNbCol());
	auto s = Eigen::sqrt(svd.singularValues().block(0, 0, r, r).array()).matrix().asDiagonal();
	bestX.mat = svd.matrixU().block(0, 0, this->getNbRow(), r) * s;
	bestY.mat = s * svd.matrixV().block(0, 0, this->getNbCol(), r).adjoint();
}

template<typename FPP>
void MatDense<FPP, Cpu>::best_low_rank(const int &r, MatDense<FPP,Cpu> &bestX, MatDense<FPP, Cpu> &bestY) const
{
#if(EIGEN_WORLD_VERSION > 3 || EIGEN_WORLD_VERSION >= 3 && EIGEN_MAJOR_VERSION >= 3)
	Eigen::JacobiSVD<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>> svd(this->mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
#else
	//#ifdef _MSC_VER
	// as far as I tested eigen3.4rc1 doesn't compile with VS 14
	// so use JacobiSVD
	Eigen::JacobiSVD<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>> svd(this->mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
#endif
	best_low_rank(r, bestX, bestY, svd);
}

template<typename FPP>
void MatDense<FPP, Cpu>::approx_rank1(MatDense<FPP,Cpu> &bestX, MatDense<FPP, Cpu> &bestY) const
{
	MatDense<FPP, Cpu> AtA, tAA;
	int flag;
	gemm(*this, *this, AtA, FPP(1.0), FPP(0.0), 'N', 'T');
	gemm(*this, *this, tAA, FPP(1.0), FPP(0.0), 'T', 'N');
	bestX.resize(this->getNbRow(), 1);
	bestY.resize(this->getNbCol(), 1);
	FPP sigma1 = std::sqrt(std::abs(power_iteration(AtA, 2, 1e-2, flag, bestX.getData(), /*rand_init*/ true)));
	FPP sigma2 = std::sqrt(std::abs(power_iteration(tAA, 2, 1e-2, flag, bestY.getData(), true)));
	MatDense<FPP,Cpu> tuav;
	gemm(*this, bestY, tuav, FPP(1.0), FPP(0.0), 'N', 'N');
	gemm(bestX, tuav, tuav, FPP(1.0), FPP(0.0), 'H', 'N');
	if(std::abs(sigma2) < 0 && tuav(0,0) > 0 || std::abs(sigma2) > 0 && tuav(0,0) < 0)
		sigma2 *= -1;
	bestY *= sigma2;
	bestY.adjoint();
	// test the error to determine if the sigma sign is correct (the way to do it above is slightly better)
	//	auto err = bestY;
	//	bestX.multiply(err, 'N');
	//	err -= *this;
	//	//				std::cout << "err:" << err.norm() << std::endl;
	//	if(err.norm() > 1e-6)
	//		bestY *= -1;

}

template<typename FPP>
void MatDense<FPP, Cpu>::initJacobiSVD(Eigen::JacobiSVD<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>>& svd)
{
	svd.compute(this->mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
}

	template<typename FPP>
std::vector<int> MatDense<FPP, Cpu>::col_nonzero_inds(faust_unsigned_int col_id) const
{
	std::vector<int> ids;
	auto data_ptr = this->getData()+col_id*this->getNbRow();
	for(int i=0;i<this->getNbRow();i++)
	{
		if((FPP)data_ptr[i] != FPP(0)) ids.push_back(i);
	}
	return ids;
}

template<typename FPP>
std::vector<int> MatDense<FPP, Cpu>::row_nonzero_inds(faust_unsigned_int row_id) const
{
	std::vector<int> ids;
	auto data_ptr = this->getData();
	for(int j=0;j<this->getNbRow();j++)
		if((FPP)(data_ptr[j*this->getNbRow()+row_id]) != FPP(0)) ids.push_back(j);
	return ids;
}

template<typename FPP>
void MatDense<FPP, Cpu>::submatrix(const std::vector<int> &row_ids, const std::vector<int> &col_ids, MatDense<FPP, Cpu> & submat) const
{
	if(this->dim1 != row_ids.size() || this->dim2 != col_ids.size())
		submat.resize(row_ids.size(), col_ids.size());
#ifdef _MSC_VER
	// as far as I tested eigen3.4rc1 doesn't compile with VS 14
	// so extract the matrix manually in this env.
	for(int i=0;i<row_ids.size();i++)
		for(int j=0;j<col_ids.size();j++)
			submat.mat(i,j) = mat(row_ids[i], col_ids[j]);
#else
	submat.mat = mat(row_ids, col_ids);
#endif
}

template<typename FPP>
void MatDense<FPP, Cpu>::set_col_coeffs(faust_unsigned_int col_id, const std::vector<int> &row_ids, const MatDense<FPP, Cpu> &values, faust_unsigned_int val_col_id)
{
	for(int i=0;i<row_ids.size();i++)
	{
		auto row_id = row_ids[i];
		mat(row_id, col_id) = values(i, val_col_id);
	}
//	this->isZeros = this->getNonZeros() == 0;
	this->isZeros = false; // too costly to verify exhaustively, in doubt say it is not 0
}

template<typename FPP>
void MatDense<FPP, Cpu>::set_row_coeffs(faust_unsigned_int row_id, const std::vector<int> &col_ids, const MatDense<FPP, Cpu> &values, faust_unsigned_int val_row_id)
{
	for(int i=0;i<col_ids.size();i++)
	{
		auto col_id = col_ids[i];
		mat(row_id, col_id) = values(val_row_id, i);
	}
//	this->isZeros = this->getNonZeros() == 0;
	this->isZeros = false; // too costly to verify exhaustively, in doubt say it is not 0
}

template<typename FPP>
FPP MatDense<FPP,Cpu>::sum_col(faust_unsigned_int id) const
{
	FPP sum = (FPP)0;
	for(int i=0;i<this->getNbRow();i++)
		sum += (*this)(i, id);
	return sum;
}

	template<typename FPP>
FPP MatDense<FPP,Cpu>::sum_row(faust_unsigned_int id) const
{
	FPP sum = (FPP)0;
	for(int i=0;i<this->getNbCol();i++)
		sum += (*this)(id, i);
	return sum;
}

	template<typename FPP>
MatDense<FPP,Cpu> MatDense<FPP,Cpu>::get_block(faust_unsigned_int i, faust_unsigned_int j, faust_unsigned_int nrows, faust_unsigned_int ncols)
{
	MatDense<FPP,Cpu> block_mat(nrows, ncols);
	block_mat.mat.block(0,0, nrows, ncols) = mat.block(i, j, nrows, ncols);
	return block_mat;
}

	template<typename FPP>
void MatDense<FPP,Cpu>::set_block(faust_unsigned_int i, faust_unsigned_int j, MatDense<FPP,Cpu> & block_mat)
{
	faust_unsigned_int nrows, ncols;
	nrows = block_mat.getNbRow();
	ncols = block_mat.getNbCol();
	mat.block(i, j, nrows, ncols) = block_mat.mat.block(0, 0, nrows, ncols);
}

template<typename FPP>
FPP MatDense<FPP,Cpu>::min_coeff() const
{
	return mat.minCoeff();
}

template<typename FPP>
Vect<FPP,Cpu> MatDense<FPP,Cpu>::rowwise_min() const
{
	return Vect<FPP,Cpu>(this->getNbCol(), mat.rowwise().minCoeff().eval().data());
}

template<typename FPP>
Vect<FPP,Cpu> MatDense<FPP,Cpu>::rowwise_min(int* col_indices) const
{
	Vect<FPP,Cpu> vec(this->getNbRow());
	for(int i=0;i<this->getNbRow();i++)
		vec.getData()[i] = mat.row(i).minCoeff(col_indices+i);
	return vec;
}

	template<typename FPP>
void kron(const MatDense<FPP,Cpu> &A, const MatDense<FPP, Cpu> &B, MatDense<FPP,Cpu>& out)
{
	auto out_nrows = A.getNbRow()*B.getNbRow();
	auto out_ncols = A.getNbCol()*B.getNbCol();
	out.resize(out_nrows, out_ncols);
	MatDense<FPP,Cpu> tmp;
	//TODO: OpenMP ?
	for(int i=0;i<A.getNbRow();i++)
	{
		for(int j=0;j<A.getNbRow();j++)
		{
			auto s = A(i,j);
			tmp = B;
			tmp *= s;
			// copy column by column in out
			for(int k=0;k<B.getNbCol();k++)
			{
				auto col_offset = out_nrows*j*B.getNbCol()+out_nrows*k+i*B.getNbRow();
				memcpy(out.getData()+col_offset, tmp.getData()+k*B.getNbRow(), sizeof(FPP)*B.getNbRow());
			}
		}
	}
}

template<typename FPP>
Vect<FPP,Cpu> MatDense<FPP,Cpu>::rowwise_max(int* col_indices) const
{
	Vect<FPP,Cpu> vec(this->getNbRow());
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
void MatDense<FPP,Cpu>::abs()
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
MatDense<FPP, Cpu>* MatDense<FPP, Cpu>::randMat(faust_unsigned_int num_rows, faust_unsigned_int num_cols)
{
	MatDense<FPP, Cpu>* mat = nullptr;
	try {
		Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic> mat_tmp;
		mat_tmp = Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>::Random(num_rows, num_cols);
		mat = new MatDense<FPP, Cpu>(mat_tmp.data(), num_rows, num_cols);
	}
	catch(std::bad_alloc e) {
		//		std::cerr << "Out of memory." << std::endl;
	}
	return mat;
}

	template<typename FPP>
MatDense<FPP, Cpu>* MatDense<FPP, Cpu>::randMat(faust_unsigned_int num_rows, faust_unsigned_int num_cols, float density)
{
	MatDense<FPP, Cpu>* mat = nullptr;
	try {
		//TODO: refactor with above randMat()
		Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic> mat_tmp;
		std::default_random_engine generator(time(nullptr));
		std::uniform_real_distribution<double> distribution(0, 1);
		mat_tmp = Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>::Random(num_rows, num_cols);
		for(int i=0;i<num_rows;i++)
			for(int j=0;j<num_cols;j++)
				if(distribution(generator) > density)
					mat_tmp(i,j) = 0;
		mat = new MatDense<FPP, Cpu>(mat_tmp.data(), num_rows, num_cols);
		mat_tmp.resize(0,0);
	}
	catch(std::bad_alloc e) {
		//		std::cerr << "Out of memory." << std::endl;
	}
	return mat;
}

	template<typename FPP>
MatDense<FPP, Cpu>* MatDense<FPP, Cpu>::randMat(faust_unsigned_int num_rows, faust_unsigned_int num_cols, float density, bool per_row)
{
	MatDense<FPP, Cpu>* mat = nullptr;
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
Vect<FPP,Cpu> MatDense<FPP, Cpu>::diagonal(int index)
{
	return gen_diagonal(index, true);
}

	template<typename FPP>
Vect<FPP,Cpu> MatDense<FPP, Cpu>::adiagonal(int index)
{
	return gen_diagonal(index, false);
}

	template<typename FPP>
Vect<FPP,Cpu> MatDense<FPP, Cpu>::gen_diagonal(int index, bool is_diag)
{
	int i, j, ei = 0;
	std::vector<std::pair<int,int>> indices = is_diag?this->get_diag_indices(index):this->get_antidiag_indices(index);
	FPP* data = new FPP[indices.size()];  //use the heap to avoid C2131 (VS14)
	for(auto ind : indices)
	{
		i = ind.first;
		j = ind.second;
		data[ei++] = *(this->getData()+i+this->getNbRow()*j);
	}
	Vect<FPP,Cpu> diag(ei, data);
	delete[] data;
	return diag;
}

	template<typename FPP>
std::vector<std::pair<int,int>> MatDense<FPP, Cpu>::get_diag_indices(int index)
{
	if((index > 0 && index > this->getNbCol()) || (index < 0 && -index > this->getNbRow())) throw std::out_of_range("diagonal index is out of range.");
	std::vector<std::pair<int, int>> indices;
	int i, j;
	if(index >= 0)
	{
		i = 0;
		j = index;
	}
	else
	{
		i = -index;
		j = 0;
	}
	while(i < this->getNbRow() && j < this->getNbCol())
	{
		indices.push_back(std::make_pair(i,j));
		i++;
		j++;
	}
	return indices;
}

	template<typename FPP>
std::vector<std::pair<int,int>> MatDense<FPP, Cpu>::get_antidiag_indices(int index)
{
	if((index > 0 && index > this->getNbRow()) || (index < 0 && -index > this->getNbCol())) throw std::out_of_range("anti-diagonal index is out of range.");
	std::vector<std::pair<int, int>> indices;
	int i, j;
	if(index >= 0)
	{
		i = this->getNbRow()-index-1;
		j = 0;
	}
	else
	{
		i = this->getNbRow()-1;
		j = -index;
	}
	while(i < this->getNbRow() && j < this->getNbCol())
	{
		indices.push_back(std::make_pair(i,j));
		i--;
		j++;
	}
	return indices;
}

	template<typename FPP>
	template<typename MatType1, typename MatType2>
void MatDense<FPP, Cpu>::eigenIndexMul(const faust_unsigned_int* row_ids, const faust_unsigned_int* col_ids, size_t nrows, size_t ncols, const MatType1 &in_mat, MatType2 &out_mat, bool transpose/* = false*/, bool conjugate /*= false*/)
{
	if(row_ids == nullptr && col_ids == nullptr)
	{
		auto tmp = mat(Eigen::all, Eigen::all);
		if(transpose && conjugate)
			out_mat = tmp.adjoint() * in_mat;
		else if(transpose)
			out_mat = tmp.adjoint() * in_mat;
		else if(conjugate)
			out_mat = tmp.conjugate() * in_mat;
		else
			out_mat = tmp * in_mat;
	}
	else if(row_ids == nullptr && col_ids != nullptr)
	{
		std::vector<int> vec_cols(col_ids, col_ids+ncols);
		if(transpose && conjugate)
			out_mat = mat(vec_cols, Eigen::all).adjoint() * in_mat;
		else if(transpose)
			out_mat = mat(vec_cols, Eigen::all).transpose() * in_mat;
		else if(conjugate)
			out_mat = mat(Eigen::all, vec_cols).conjugate() * in_mat;
		else
			out_mat = mat(Eigen::all, vec_cols) * in_mat;
	}
	else if(row_ids != nullptr && col_ids == nullptr)
	{

		std::vector<int> vec_rows(row_ids, row_ids+nrows);
		if(transpose && conjugate)
			out_mat = mat(Eigen::all, vec_rows).adjoint() * in_mat;
		else if(transpose)
			out_mat = mat(Eigen::all, vec_rows).transpose() * in_mat;
		else if(conjugate)
			out_mat = mat(vec_rows, Eigen::all).conjugate() * in_mat;
		else
			out_mat = mat(vec_rows, Eigen::all) * in_mat;
	}
	else // if(row_ids != nullptr && col_ids != nullptr)
	{

		std::vector<int> vec_rows(row_ids, row_ids+nrows);
		std::vector<int> vec_cols(col_ids, col_ids+ncols);
		if(transpose && conjugate)
			out_mat = mat(vec_cols, vec_rows).adjoint() * in_mat;
		else if(transpose)
			out_mat = mat(vec_cols, vec_rows).transpose() * in_mat;
		else if(conjugate)
			out_mat = mat(vec_rows, vec_cols).conjugate() * in_mat;
		else
			out_mat = mat(vec_rows, vec_cols) * in_mat;
	}
}

#ifdef __COMPILE_TIMERS__
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_constr;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_get_coeff;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_get_coeffs;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_set_coeff;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_set_coeffs;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_set_coeffs2;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_resize;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_check_dim;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_max;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_transpose;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_mult_right;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_mult_left;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_scalar_multiply;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_add;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_sub;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_print_file;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_multiply;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_gemm;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_add_ext;

template<typename FPP>
Timer MatDense<FPP,Cpu>::t_spectral_norm;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_spectral_norm2;
template<typename FPP>
Timer MatDense<FPP,Cpu>::t_power_iteration;

template<typename FPP>
void MatDense<FPP,Cpu>::print_timers()const
{
	std::cout << "timers in MatDense :" << std::endl;
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


	std::cout << "timers in MatDense / LinearAlgebra :" << std::endl;
	std::cout << "t_multiply        = " << t_multiply.get_time()        << " s for "<< t_multiply.get_nb_call()        << " calls" << std::endl;
	std::cout << "t_gemm            = " << t_gemm.get_time()            << " s for "<< t_gemm.get_nb_call()            << " calls" << std::endl;
	std::cout << "t_add_ext         = " << t_add_ext.get_time()         << " s for "<< t_add_ext.get_nb_call()         << " calls" << std::endl<<std::endl<<std::endl;
}
#endif


}
#endif
