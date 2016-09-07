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
/****************************************************************************/
#ifndef __FAUST_MAT_HPP__
#define __FAUST_MAT_HPP__



#include "faust_linear_algebra.h"

//#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Sparse>
#include "faust_exception.h"
#include "faust_constant.h"

#ifdef __GEMM_WITH_OPENBLAS__
	#include "faust_cblas_algebra.h"
#endif




template<typename FPP>
const char * Faust::MatDense<FPP,Cpu>::m_className = "Faust::MatDense<FPP,Cpu>::";


template<typename FPP>
Faust::MatDense<FPP,Cpu>::MatDense(const FPP  *data_,const faust_unsigned_int nbRow, const faust_unsigned_int nbCol ) : Faust::MatGeneric<FPP,Cpu>(nbRow,nbCol), mat(nbRow,nbCol), isIdentity(false),isZeros(false)
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
    isIdentity = false;

#ifdef __COMPILE_TIMERS__
t_resize.stop();
#endif
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
void Faust::MatDense<FPP,Cpu>::setZeros()
{
	memset(getData(), 0, sizeof(FPP) * this->dim1*this->dim2);
	isZeros = true;
	isIdentity = false;
}

template<typename FPP>
void Faust::MatDense<FPP,Cpu>::setEyes()
{
	setZeros();
	FPP* ptr_data = getData();
	for (int i=0 ; i<std::min(this->dim1,this->dim2); i++)
		ptr_data[i*this->dim1+i] = 1.0;
	if (this->dim1 == this->dim2)
		isIdentity = true;
	isZeros = false;
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

	if(isZeros || isIdentity)
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
void Faust::MatDense<FPP,Cpu>::multiplyRight(Faust::MatDense<FPP,Cpu> const& A)
{

#ifdef __COMPILE_TIMERS__
t_mult_right.start();
#endif

	if (this->dim2 != A.dim1)
	{
		handleError(m_className, "multiplyRight : dimension conflict between matrix");
	}

	if(A.isIdentity)
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
		isIdentity = false;
		#ifdef __COMPILE_TIMERS__
			t_mult_right.stop();
		#endif
		return;
	}

	if(isIdentity)
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

	if(isIdentity)
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




template<typename FPP>
FPP Faust::MatDense<FPP,Cpu>::spectralNorm(const faust_unsigned_int nbr_iter_max,FPP threshold, faust_int & flag,Faust::BlasHandle<Cpu> blas_handle/*=Faust::BlasHandle<Cpu>()*/) const
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

	if(isIdentity)
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
		gemm((*this),(*this),AtA,(FPP)1.,(FPP)0.0,'N','T');
	}else
	{
		gemm((*this),(*this),AtA,(FPP)1.,(FPP)0.0,'T','N');
	}

	FPP  res=std::sqrt(power_iteration(AtA,nbr_iter_max,threshold,flag));


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
    isIdentity = false;
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
		handleError(m_className, "sub : matrix dimension not equal");
	}
	mat = mat - A.mat;

    isZeros = false;
    isIdentity = false;

#ifdef __COMPILE_TIMERS__
t_sub.stop();
#endif
}




// Affichage
template<typename FPP>
void Faust::MatDense<FPP,Cpu>::Display() const
{
    std::cout << "nb_row=" << this->getNbRow() << endl;
    std::cout << "nb_col=" << this->getNbCol()   <<endl;
    if(isZeros)
        cout << this->dim1 << " by " << this->dim2 << " matrix of zeros" << endl;
    else if (isIdentity)
        cout << this->dim1 << " by " << this->dim2 << " identity matrix" << endl;
    else if (this->dim1*this->dim2==0)
        cout << this->dim1 << " by " << this->dim2 << " empty matrix" << endl;
    else
    {
        for (int i=0 ; i<this->dim1 ; i++)
        {
            for(int j=0 ; j<this->dim2 ; j++)
                cout << (*this)(i,j) << " " ;
            cout << endl;
        }
    }
    cout<<endl;
}

template<typename FPP>
void Faust::MatDense<FPP,Cpu>::print_file(const char* filename)const
{
    ofstream fichier;
	fichier.open(filename);
	for (int i=0 ; i<this->getNbRow() ; i++)
	{
		for (int j=0 ; j<this->getNbCol() ; j++)
			fichier << setprecision(20) <<mat(i,j) << " ";
		fichier << endl;
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
    isIdentity = A.isIdentity;
}

template<typename FPP>
template<typename FPP1>
void Faust::MatDense<FPP,Cpu>::operator=(Faust::MatDense<FPP1,Cpu> const& A)
{
    resize(A.dim1,A.dim2);
	for (int i=0;i<this->dim1*this->dim2;i++)
        (*this)[i]=(FPP) A(i);
    isZeros = A.isZeros;
    isIdentity = A.isIdentity;
}

template<typename FPP>
void Faust::MatDense<FPP,Cpu>::operator=(Faust::MatSparse<FPP,Cpu> const& S)
{
    S.check_dim_validity();
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

	isZeros = false;
	isIdentity = false;
}


template<typename FPP>
void Faust::MatDense<FPP,Cpu>::operator*=(const Faust::MatSparse<FPP,Cpu>& S)
{
	if(this->dim2 != S.dim1)
	{
		handleError(m_className, "operator*= : incorrect matrix dimensions");
	}

	if (isIdentity)
	{
		this->operator=(S);
		isIdentity = false;
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
	if(this->dim1!=S.dim1 || this->dim2!=S.dim2)
	{
		handleError(m_className,"operator+= : incorrect matrix dimensions");
	}
	mat += S.mat;
	isIdentity = false;
	isZeros = false;
}

template<typename FPP>
void Faust::MatDense<FPP,Cpu>::operator-=(const Faust::MatSparse<FPP,Cpu>& S)
{
	if(this->dim1!=S.dim1 || this->dim2!=S.dim2)
	{
		handleError(m_className,"operator-= : incorrect matrix dimensions");
	}
	mat -= S.mat;
	isIdentity = false;
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
	isIdentity = false;
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

	if (isIdentity)
	{
		this->operator=(S);
		isIdentity = false;
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

    ifstream* vec_stream;
    vec_stream = new ifstream(filename);
    if (!vec_stream->is_open())
        handleError(m_className, "init_from_file : unable to open file");
    istream_iterator<FPP> start(*vec_stream), eos;
    vector<FPP> vec(start, eos);

    if((vec[0]*vec[1]+2) != vec.size())
    {
        handleError(m_className, "init_from_file : problem with the file");
    }
    resize(vec[0],vec[1]);
    memcpy(getData(), &vec[2], sizeof(FPP) * this->dim1 * this->dim2);

    isZeros = false;
    isIdentity = false;

#ifdef __COMPILE_TIMERS__
t_print_file.stop();
#endif
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
   cout << "timers in Faust::MatDense :" << endl;
   cout << "t_constr          = " << t_constr.get_time()          << " s for "<< t_constr.get_nb_call()          << " calls" << endl;
   cout << "t_get_coeff       = " << t_get_coeff.get_time()       << " s for "<< t_get_coeff.get_nb_call()       << " calls" << endl;
   cout << "t_get_coeffs      = " << t_get_coeffs.get_time()      << " s for "<< t_get_coeffs.get_nb_call()      << " calls" << endl;
   cout << "t_set_coeff       = " << t_set_coeff.get_time()       << " s for "<< t_set_coeff.get_nb_call()       << " calls" << endl;
   cout << "t_set_coeffs      = " << t_set_coeffs.get_time()      << " s for "<< t_set_coeffs.get_nb_call()      << " calls" << endl;
   cout << "t_set_coeffs2     = " << t_set_coeffs2.get_time()     << " s for "<< t_set_coeffs2.get_nb_call()     << " calls" << endl;
   cout << "t_resize          = " << t_resize.get_time()          << " s for "<< t_resize.get_nb_call()          << " calls" << endl;
   cout << "t_check_dim       = " << t_check_dim.get_time()       << " s for "<< t_check_dim.get_nb_call()       << " calls" << endl;
   cout << "t_max             = " << t_max.get_time()             << " s for "<< t_max.get_nb_call()             << " calls" << endl;
   cout << "t_transpose       = " << t_transpose.get_time()       << " s for "<< t_transpose.get_nb_call()       << " calls" << endl;
   cout << "t_mult_right      = " << t_mult_right.get_time()      << " s for "<< t_mult_right.get_nb_call()      << " calls" << endl;
   cout << "t_mult_left       = " << t_mult_left.get_time()       << " s for "<< t_mult_left.get_nb_call()       << " calls" << endl;
  cout << "t_scalar_multiply = " << t_scalar_multiply.get_time() << " s for "<< t_scalar_multiply.get_nb_call() << " calls" << endl;
   cout << "t_add             = " << t_add.get_time()             << " s for "<< t_add.get_nb_call()             << " calls" << endl;
   cout << "t_sub             = " << t_sub.get_time()             << " s for "<< t_sub.get_nb_call()             << " calls" << endl;
cout << "t_print_file      = " << t_print_file.get_time()      << " s for "<< t_print_file.get_nb_call()      << " calls" << endl<<endl;


   cout << "timers in Faust::MatDense / LinearAlgebra :" << endl;
   cout << "t_multiply        = " << t_multiply.get_time()        << " s for "<< t_multiply.get_nb_call()        << " calls" << endl;
   cout << "t_gemm            = " << t_gemm.get_time()            << " s for "<< t_gemm.get_nb_call()            << " calls" << endl;
   cout << "t_add_ext         = " << t_add_ext.get_time()         << " s for "<< t_add_ext.get_nb_call()         << " calls" << endl<<endl<<endl;
}
#endif

#endif
