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
#ifndef FAUST_MATDENSE_H
#define FAUST_MATDENSE_H


#include <Eigen/Dense>

#include "faust_constant.h"
#include <vector>
#include <iterator>
#include "faust_MatGeneric.h"
#include "faust_exception.h"
#include <iostream>

#ifdef __COMPILE_TIMERS__
    #include "faust_Timer.h"
#endif

// modif AL AL
#include "faust_Vect.h"
#include "faust_Transform.h"

#include "faust_BlasHandle.h"
#include "matio.h"

/*! \class Faust::MatDense MatDenseDense.h
* \brief Class template representing dense matrix <br>
* This class implements basic linear algebra operation (addition, multiplication, frobenius and spectral norm...) <br>
* The matrix format is ColMajor. <br>
* \tparam T scalar numeric type, e.g float or double
*/

//template<typename FPP, Device DEVICE> class MatDense;
//template<typename FPP, Device DEVICE> class MatSparse;
//template<typename FPP, Device DEVICE> class Vect;
//template<typename FPP, Device DEVICE> class Transform;
//
////! \fn add
////! \brief (*this) = (*this) + A
//template<typename FPP>
//void Faust::add(const Faust::MatDense<FPP,Cpu> & A, const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C);
//
////! \fn gemm_core
////! \brief performs ??
//template<typename FPP>
//void Faust::gemm_core(const Faust::MatDense<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB);
//
////! \fn gemm
////! \brief performs ??
//template<typename FPP>
//void Faust::gemm(const Faust::MatDense<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB);
//
////! \fn multiply
////! \brief performs ??
//template<typename FPP>
//void Faust::multiply(const Faust::MatDense<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C);
//template<typename FPP>
//void Faust::multiply(const Faust::Transform<FPP,Cpu> & A, const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, char typeA, char typeMult);
//
////! \fn Faust::gemv
////! \brief performs ??
//template<typename FPP>
//void Faust::gemv(const Faust::MatDense<FPP,Cpu> & A,const Faust::Vect<FPP,Cpu> & x,Faust::Vect<FPP,Cpu> & y,const FPP & alpha, const FPP & beta, char typeA);
//
////! \fn Faust::spgemm
////! \brief performs ??
//template<typename FPP>
//void Faust::spgemm(const Faust::MatSparse<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char typeA, char typeB);
//

//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{


template<typename FPP, Device DEVICE> class MatDense;
template<typename FPP, Device DEVICE> class MatSparse;
template<typename FPP, Device DEVICE> class Vect;
template<typename FPP, Device DEVICE> class Transform;

//! \fn add
//! \brief (*this) = (*this) + A
template<typename FPP>
void add(const Faust::MatDense<FPP,Cpu> & A, const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C);

//! \fn gemm_core
//! \brief performs ??
template<typename FPP>
void gemm_core(const Faust::MatDense<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB);

//! \fn gemm
//! \brief performs ??
template<typename FPP>
void gemm(const Faust::MatDense<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB);

//! \fn multiply
//! \brief performs ??
template<typename FPP>
void multiply(const Faust::MatDense<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C);
template<typename FPP>
void multiply(const Faust::Transform<FPP,Cpu> & A, const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, char typeA, char typeMult);

//! \fn Faust::gemv
//! \brief performs ??
template<typename FPP>
void gemv(const Faust::MatDense<FPP,Cpu> & A,const Faust::Vect<FPP,Cpu> & x,Faust::Vect<FPP,Cpu> & y,const FPP & alpha, const FPP & beta, char typeA);

//! \fn Faust::spgemm
//! \brief performs ??
template<typename FPP>
void spgemm(const Faust::MatSparse<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char typeA, char typeB);







    template<typename FPP, Device DEVICE>
    class MatDense;

    template<typename FPP, Device DEVICE>
    class MatGeneric;

    template<typename FPP, Device DEVICE>
    class Transform;
    //template<Device DEVICE> class BlasHandle;

    template<typename FPP>
    class MatDense<FPP,Cpu> : public Faust::MatGeneric<FPP,Cpu>
    {

	friend class MatSparse<FPP,Cpu>;
        
	/// All derived class template of MatDense are considered as friends
        template<class,Device> friend class MatDense;

        public:
        static const char * name;
        //void gemm(const MatDense<FPP,Cpu> & A,const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char  typeA, char  typeB);

        //! \fn faust_gemm
        //! \brief performs ??
        //!
        void faust_gemm(const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char  typeA, char  typeB)const {gemm<FPP>((*this),B,C,alpha,beta,typeA,typeB);}	/*!
        *  \brief Constructor MatDense
        *	\tparam data : pointer to the data array of the matrix
            \tparam nbRow : number of row of the matrix
            \tparam nbCol : number of column of the matrix
        */
        MatDense(const FPP  *data_,const faust_unsigned_int nbRow, const faust_unsigned_int nbCol );
        MatDense() : MatGeneric<FPP,Cpu>(), mat(0,0), isIdentity(false), isZeros(false) {}
        /*!
        *  \brief Copy Constructor of MatDense
        *  \tparam A : another MatDense
        */
        MatDense(const MatDense<FPP,Cpu> & A) : MatGeneric<FPP,Cpu>(A.dim1,A.dim2), mat(A.mat), isIdentity(A.isIdentity), isZeros(A.isZeros) {}
        template<typename FPP1>
        MatDense(const MatDense<FPP1,Cpu> & A){this->operator=(A);}
        template<typename FPP1>
        MatDense(const Faust::MatSparse<FPP1,Cpu> & A){this->operator=(A);}
        MatDense(const Faust::MatSparse<FPP,Cpu> & A){this->operator=(A);}

        MatDense(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol) : MatGeneric<FPP,Cpu>(nbRow,nbCol), mat(nbRow,nbCol), isIdentity(false), isZeros(false){}
        MatDense(const faust_unsigned_int nbRow) : MatGeneric<FPP,Cpu>(nbRow,nbRow), mat(nbRow,nbRow), isIdentity(false), isZeros(false){}
        /// Destructor of MatDense
        ~MatDense(){resize(0,0);/*std::cout<<"destructor dense mat"<<std::endl;*/}

	
	///**** METHOD INHERITED FROM VIRTUAL MATGENERIC ****///
	
		
	MatType getType() const
	{ return Dense;}
	

	/*!  \brief return a "copy" to the given matrix
	*  \param isOptimize (optionnal) : boolean which the style of copy <br>
			       -True, the return copy is optimized for the product <br>
	                       which means dynamic type of the copy could be different from the original one <br>
	                      -False, the return copy is simple, the dynamic type isn't changed <br>   						(default value is False) <br>					
        //! \return  a pointer of MatGeneric
	//  \warning the dynamic type of the copy can be different from the original object
	*/
	MatGeneric<FPP,Cpu>* Clone(const bool isOptimize=false) const;
	
	
	//\brief : return the number of non-zeros coefficient nnz
	faust_unsigned_int getNonZeros()const;

	
	//! \brief compute MatGeneric-vector multiplication
	//! \param vec : the vector
	//! \param opThis : character	
	//! vec = (*this) * vec if opThis='N'
	// vec = (*this)' * vec if opThis='T'
	void multiply(Faust::Vect<FPP,Cpu> & vec,const char opThis) const
	{Faust::gemv((*this),vec,vec,(FPP) 1.0, (FPP) 0.0,opThis);}




	//! \brief compute MatGeneric-MatDense multiplication
	//! \param M : the dense matrix
	//! \param opThis : character	
	//! M = (*this) * M if opThis='N'
	//  M = (*this)' * M if opThis='T' 
	void multiply(Faust::MatDense<FPP,Cpu> & M, const char opThis) const
	{Faust::gemm<FPP>((*this),M,M,1.0,0.0,opThis,'N');}


 	///********* FIN METHOD INHERITED FROM VIRTUAL MATGENERIC *********///
	

        //! \brief resize the MatDense
        //! \param	   nbRow : new number of row of the matrix
        //! \param	   nbCol : new number of column of the matrix
        //! \warning   nbRow and nbCol must be greater or equal to 0
        void resize(const faust_unsigned_int nbRow,const faust_unsigned_int nbCol);

        //! \brief resize the MatDense
        //! \param	 nbRow : new number of row and column of the matrix
        //! \warning nbRow must be greater or equal to 0
        void resize(const faust_unsigned_int nbRow){resize(nbRow,nbRow);}
	



	
        //! \brief Check if the dimension of the matrix are consistent, if not throws an error
        void check_dim_validity();

        //! \brief Set the matrix to the zero matrix
        void setZeros();

        //! \brief Set the matrix to the one diagonal matrix
        void setEyes();

        //! \brief Access to the ith coefficient of the matrix pointer, Colmajor format and zero indexing
        //! \param i : position
        //! \return : ith coefficient of the matrix
        FPP& operator[](faust_unsigned_int i){isZeros=false; isIdentity=false;return mat.data()[i];}

        //!  \brief Access to the ith coefficient of the matrix pointer, Colmajor format and zero indexing (version with "[]" )
        //! \param i : position
        //! \return  : read-only ith coefficient of the matrix
        const FPP& operator[](faust_unsigned_int i)const{return mat.data()[i];}

        //!  \brief Access to the ith coefficient of the matrix pointer, Colmajor format and zero indexing (version with "()" )
        //! \param i : position
        //! \return : read-only ith coefficient of the matrix
        const FPP& operator()(faust_unsigned_int i)const{return mat.data()[i];}

        //! \brief Access to (i,j) coefficient of the matrix pointer in zero indexing
        //! \param i : row position
        //! \param j : col position
        //! \return : read-only (i,j) coefficient of the matrix
        const FPP& operator()(faust_unsigned_int i, faust_unsigned_int j)const{return mat.data()[j*this->dim1+i];}

        void operator*=(const Faust::MatSparse<FPP,Cpu>& M);
        

	void operator+=(const Faust::MatSparse<FPP,Cpu>& M);

	// modif NB : v1102 comment useless function
	/*        
	void operator-=(const Faust::MatSparse<FPP,Cpu>& M);
	*/
        
	void multiplyLeft(const Faust::MatSparse<FPP,Cpu>& S,const char TransS='N');

        FPP* getData(){isZeros=false; isIdentity=false;return mat.data();}
        const FPP* getData()const{return mat.data();}

        bool isEqual(const MatDense<FPP,Cpu> & B) const;
        bool isEqual(const MatDense<FPP,Cpu> & B, FPP threshold) const;

        //!  \brief Initialize MatDense from text file
        //! \param filename : name of the file
        //! The first line of the file contains 2 integers : the number of row and the number of column. <br>
        //! All the other line contains one coefficient in ColMajor access
        void init_from_file(const char* filename);

        //! Absolute value calculated from the library eigen.
        void abs() {mat=mat.cwiseAbs();}


        //!  \brief Compute the Frobenius norm of the MatDense
        //! \return  the Frobenius norm
        FPP norm() const {return mat.norm();}

        //!  \brief Normalize the matrix according to its Frobenius norm
        void normalize() {scalarMultiply(1.0/norm());}


        //!	\param nbr_iter_max : maximum number of iteration for the power algo
        //! \param threshold : threshold until convergence
        //! \param flag : convergence flag
        //! \return Return the estimated spectral norm (maximum singular value in absolute value) using power iteration algorithm
        //! See also, template<typename FPP> FPP power_iteration(const MatDense<FPP,Cpu> & A, const faust_unsigned_int nbr_iter_max,FPP threshold,faust_int & flag);
        FPP spectralNorm(const faust_unsigned_int nbr_iter_max,FPP threshold, faust_int & flag,Faust::BlasHandle<Cpu>  blas_handle=Faust::BlasHandle<Cpu>()) const;

        //! \brief Compute the trace of the MatDense
        //! \return  the trace
        FPP trace() const {return mat.trace();}

        //! \brief Transpose the MatDense
        void transpose();

		//! \brief Replaces the matrix by its conjugate.
		void conjugate();

        //! \brief Replace this by (this) * A
        void multiplyRight(MatDense<FPP,Cpu> const& A);




        //!  \brief replace this by lambda * (*this)
        void scalarMultiply(FPP const lambda);

        //!  \brief replace this by lambda * (*this) using element by element multiplication
        void scalarMultiply(MatDense<FPP,Cpu> const& A);

        //! \brief (*this) = (*this) + A
        void add(MatDense<FPP,Cpu> const& A);

        //!  \brief (*this) = (*this) - A
        void sub(MatDense<FPP,Cpu> const& A);

        //! \brief Displays the MatDense
        void Display() const;

        //!  \brief Write MatDense into text file
        //! \param filename : name of the file
        //!
        //! The first line of the file contains 2 integer : the number of row and the number of column
        //! All the other line contains one coefficient in ColMajor access of the MatDense
        void print_file(const char* filename)const;

		matvar_t* toMatIOVar(bool transpose, bool conjugate) const;

        void operator=(MatDense<FPP,Cpu> const& A);

        template<typename FPP1>
        void operator=(MatDense<FPP1,Cpu> const& A);
        template<typename FPP1>
        void operator=(Faust::MatSparse<FPP1,Cpu> const& A){Faust::MatSparse<FPP,Cpu> AT(A);this->operator=(AT);};

        void operator=(Faust::MatSparse<FPP,Cpu> const& A);
        void operator-=(MatDense<FPP,Cpu> const& A){sub(A);}
        void operator+=(MatDense<FPP,Cpu> const& A){add(A);}
        void operator*=(MatDense<FPP,Cpu> const& A){multiplyRight(A);}

        void operator*=(FPP lambda){scalarMultiply(lambda);}
        void operator/=(FPP lambda){scalarMultiply(1.0/lambda);}


	


//        friend void gemm_core<>(const MatDense<FPP,Cpu> & A,const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB);
//        friend void multiply<>(const MatDense<FPP,Cpu> & A,const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C);
//        friend void spgemm<>(const Faust::MatSparse<FPP,Cpu> & A,const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char  typeA, char  typeB);
//        friend void multiply<>(const Faust::Transform<FPP,Cpu> & A, const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, char typeA, char typeMult);
//        friend void Faust::gemv<>(const MatDense<FPP,Cpu> & A,const Faust::Vect<FPP,Cpu> & x,Faust::Vect<FPP,Cpu> & y,const FPP & alpha, const FPP & beta, char typeA);
// modif AL AL
        friend void Faust::gemm_core<>(const MatDense<FPP,Cpu> & A,const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB);
        friend void Faust::multiply<>(const MatDense<FPP,Cpu> & A,const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C);
        friend void Faust::spgemm<>(const Faust::MatSparse<FPP,Cpu> & A,const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char  typeA, char  typeB);
        friend void Faust::multiply<>(const Faust::Transform<FPP,Cpu> & A, const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, char typeA, char typeMult);
        friend void Faust::gemv<>(const MatDense<FPP,Cpu> & A,const Faust::Vect<FPP,Cpu> & x,Faust::Vect<FPP,Cpu> & y,const FPP & alpha, const FPP & beta, char typeA);
//	friend void  Faust::MatSparse<FPP,Cpu>::multiply(MatDense<FPP,Cpu> & M,const char opThis) const;

        bool estIdentite()const{return isIdentity;}
        bool estNulle()const{return isZeros;}

        private:
        Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic> mat;
        bool isIdentity;
        bool isZeros;
        static const char * m_className;


    #ifdef __COMPILE_TIMERS__
        public:
        Faust::Timer t_local_muliplyLeft;
        //temporary members
        static Faust::Timer t_constr;
        static Faust::Timer t_get_coeff;
        static Faust::Timer t_get_coeffs;
        static Faust::Timer t_set_coeff;
        static Faust::Timer t_set_coeffs;
        static Faust::Timer t_set_coeffs2;
        static Faust::Timer t_resize;
        static Faust::Timer t_check_dim;
        static Faust::Timer t_max;
        static Faust::Timer t_transpose;
        static Faust::Timer t_mult_right;
        static Faust::Timer t_mult_left;
        static Faust::Timer t_scalar_multiply;
        static Faust::Timer t_add;
        static Faust::Timer t_sub;
        static Faust::Timer t_print_file;
        static Faust::Timer t_spectral_norm;
        static Faust::Timer t_spectral_norm2;

        static Faust::Timer t_power_iteration;
        static Faust::Timer t_multiply;
        static Faust::Timer t_gemm;
        static Faust::Timer t_add_ext;
	
        void print_timers()const;
    #endif
    };
    
     /// TEMPLATE SPECIALIZATION ////	
    //template specialization member function must be done at namespace scope 
    // advice : put the specialization in the same file as the declaration to avoid issue	
    
    /* this specialisation is for raising an error for incompatible Matrix operator=
       for instance a real scalar matrix can't be assigned to a complex scalar matrix*/ 		
    /*
    template<>
    template<>
    void Faust::MatDense<double,Cpu>::operator=(Faust::MatDense<std::complex<double>,Cpu> const& A)
    {
    	handleError("MatDense<double,Cpu>","operator=(MatDense<std::complex<double>,Cpu> impossible to pass from complex to real");
    }


    template<>
    template<>
    void Faust::MatDense<float,Cpu>::operator=(Faust::MatDense<std::complex<float>,Cpu> const& A)
    {
    	handleError("MatDense<float,Cpu>","operator=(MatDense<std::complex<float>,Cpu> impossible to pass from complex to real");
    }

    template<>
    template<>
    void Faust::MatDense<double,Cpu>::operator=(Faust::MatDense<std::complex<float>,Cpu> const& A)
    {
    	handleError("MatDense<double,Cpu>","operator=(MatDense<std::complex<float>,Cpu> impossible to pass from complex to real");
    }				

    template<>
    template<>
    void Faust::MatDense<float,Cpu>::operator=(Faust::MatDense<std::complex<double>,Cpu> const& A)
    {
    	handleError("MatDense<float,Cpu>","operator=(MatDense<std::complex<double>,Cpu> impossible to pass from complex to real");
    }		
    */


}

#include "faust_MatDense.hpp"


#endif
