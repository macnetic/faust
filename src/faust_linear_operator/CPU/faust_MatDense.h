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
#ifndef FAUST_MATDENSE_H
#define FAUST_MATDENSE_H


#include <Eigen/Dense>

#include "faust_constant.h"
#include <vector>
#include <list>
#include <iterator>
#include "faust_MatGeneric.h"
#include "faust_exception.h"
#include <iostream>
#include <utility>
#ifdef __COMPILE_TIMERS__
    #include "faust_Timer.h"
#endif

#include "faust_Vect.h"
#include "faust_Transform.h"

#include "faust_linear_algebra.h"
#ifndef NO_MATIO
#include "matio.h"  // otherwise matvar_t is defined as void from MatGeneric header
#endif
#include <random>
/*! \class Faust::MatDense MatDenseDense.h
* \brief Class template representing dense matrix <br>
* This class implements basic linear algebra operation (addition, multiplication, frobenius and spectral norm...) <br>
* The matrix format is ColMajor. <br>
* \tparam T scalar numeric type, e.g float or double
*/

//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{


	template<typename FPP, FDevice DEVICE> class MatDense;
	template<typename FPP, FDevice DEVICE> class MatSparse;
	template<typename FPP, FDevice DEVICE> class Vect;
	template<typename FPP, FDevice DEVICE> class Transform;
	template<typename FPP, FDevice DEVICE> class TransformHelper;

	//! \fn add
	//! \brief (*this) = (*this) + A
	template<typename FPP>
		void add(const MatDense<FPP,Cpu> & A, const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C);

	//! \fn gemm_core
	//! \brief performs ??
	template<typename FPP>
		void gemm_core(const MatDense<FPP,Cpu> & A,const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB);

	//! \fn gemm
	//! \brief performs ??
	template<typename FPP>
		void gemm(const MatDense<FPP,Cpu> & A,const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB);

	//! \fn multiply
	//! \brief performs ??
	template<typename FPP>
		void multiply(const MatDense<FPP,Cpu> & A,const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C);
	template<typename FPP>
		void multiply(const Transform<FPP,Cpu> & A, const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, char typeA, char typeMult);

	//! \fn gemv
	//! \brief performs ??
	template<typename FPP>
		void gemv(const MatDense<FPP,Cpu> & A,const Vect<FPP,Cpu> & x,Vect<FPP,Cpu> & y,const FPP & alpha, const FPP & beta, char typeA);

	//! \fn spgemm
	//! \brief performs ??
	template<typename FPP>
		void spgemm(const MatSparse<FPP,Cpu> & A,const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char typeA, char typeB);

	template<typename FPP>
		void spgemm(const MatDense<FPP,Cpu> & A,const MatSparse<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char typeA, char typeB);


	template<typename FPP, FDevice DEVICE>
		class MatDense;

	template<typename FPP, FDevice DEVICE>
		class MatGeneric;

	template<typename FPP, FDevice DEVICE>
		class Transform;

	template<typename FPP>
		class MatDiag;

	template<typename FPP,FDevice DEVICE>
		class MatBSR;

	template<typename FPP, FDevice DEV> class MatButterfly;

	template<typename FPP>
		class MatDense<FPP,Cpu> : public MatGeneric<FPP,Cpu>
		{
			using EigDenseMat = Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>;

			friend class MatSparse<FPP,Cpu>;
			friend class MatBSR<FPP, Cpu>;
			friend TransformHelper<FPP,Cpu>; // TODO: limit to needed member functions only
			friend Transform<FPP,Cpu>; //TODO: limit to needed member functions only (multiply)
			friend void  MatDiag<FPP>::multiply(MatDense<FPP,Cpu> & M, char opThis) const;
			friend MatButterfly<FPP, Cpu>;

			/// All derived class template of MatDense are considered as friends
			template<class,FDevice> friend class MatDense;

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
			explicit MatDense(const FPP  *data_,const faust_unsigned_int nbRow, const faust_unsigned_int nbCol );
			explicit MatDense(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, const FPP  *data_);
			MatDense() : MatGeneric<FPP,Cpu>(), mat(0,0), isZeros(false) {}
			/*!
			 *  \brief Copy Constructor of MatDense
			 *  \tparam A : another MatDense
			 */
			MatDense(const MatDense<FPP,Cpu> & A) : MatGeneric<FPP,Cpu>(A.dim1,A.dim2, A.is_ortho, A.is_identity), mat(A.mat), isZeros(A.isZeros) { }
			template<typename FPP1>
				MatDense(const MatDense<FPP1,Cpu> & A){this->operator=(A);}
			template<typename FPP1>
				MatDense(const MatSparse<FPP1,Cpu> & A)  {this->operator=(A);this->is_ortho = A.is_ortho;}
			MatDense(const MatSparse<FPP,Cpu> & A)  {this->operator=(A);this->is_ortho = A.is_ortho;}

			explicit MatDense(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol) : MatGeneric<FPP,Cpu>(nbRow,nbCol), mat(nbRow,nbCol), isZeros(false){}
			explicit MatDense(const faust_unsigned_int nbRow) : MatGeneric<FPP,Cpu>(nbRow,nbRow), mat(nbRow,nbRow), isZeros(false){}
			/// Destructor of MatDense
			~MatDense(){resize(0,0);/*std::cout<<"destructor dense mat"<<std::endl;*/}


			///**** METHOD INHERITED FROM VIRTUAL MATGENERIC ****///


			MatType getType() const
			{ return MatType::Dense;}


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


			//! \brief compute MatDense-vector multiplication
			//! \param vec : the vector
			//! \param opThis : character	
			//! vec = (*this) * vec if opThis='N'
			// vec = (*this)' * vec if opThis='T'
			void multiply(Vect<FPP,Cpu> & vec,const char opThis) const
			{gemv((*this),vec,vec,(FPP) 1.0, (FPP) 0.0,opThis);}

			Vect<FPP,Cpu> multiply(const Vect<FPP,Cpu> & vec) const
			{
				return *this * vec;
			}


			//! \brief compute MatDense multiplication
			//! \param M : the dense matrix
			//! \param opThis : character	
			//! M = (*this) * M if opThis='N'
			//  M = (*this)' * M if opThis='T' 
			void multiply(MatDense<FPP,Cpu> & M, char opThis) const
			{gemm<FPP>((*this),M,M,1.0,0.0,opThis,'N');}

			void multiply(MatSparse<FPP,Cpu> & M, char opThis) const;


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

			//! \brief Set all matrix coefficients to one.
			void setOnes();

			//! \brief Set the matrix to the zero matrix
			void setZeros();

			//! \brief Set the matrix to the one diagonal matrix
			void setEyes();

			// \brief Sets all nonzeros to one.
			void setNZtoOne();

			// \brief Sets the matrix to random values.
			// \note using this function is preferable instead of using randMat functions.
			void setRand();

			//! \brief Returns the identity matrix.
			static MatDense<FPP,Cpu> eye(faust_unsigned_int nrows, faust_unsigned_int ncols);


			//! \brief Access to the ith coefficient of the matrix pointer, Colmajor format and zero indexing
			//! \param i : position
			//! \return : ith coefficient of the matrix
			FPP& operator[](faust_unsigned_int i){isZeros=false;this->is_identity=false;return mat.data()[i];}

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

			void operator*=(const MatSparse<FPP,Cpu>& M);


			void operator+=(const MatSparse<FPP,Cpu>& M);

			// modif NB : v1102 comment useless function
			/*        
				  void operator-=(const MatSparse<FPP,Cpu>& M);
				  */

			void multiplyLeft(const MatSparse<FPP,Cpu>& S,const char TransS='N');

			FPP* getData(){isZeros=false;this->is_identity=false;return mat.data();}
			const FPP* getData()const{return mat.data();}

			void setData(const FPP* data, int32_t nrows, int32_t ncols);

			bool isEqual(const MatDense<FPP,Cpu> & B) const;
			bool isEqual(const MatDense<FPP,Cpu> & B, FPP threshold) const;

			//!  \brief Initialize MatDense from text file
			//! \param filename : name of the file
			//! The first line of the file contains 2 integers : the number of row and the number of column. <br>
			//! All the other line contains one coefficient in ColMajor access
			void init_from_file(const char* filename);

			//! Absolute elementwise value
			void abs();

			//!  \brief Compute the Frobenius norm of the MatDense
			//! \return  the Frobenius norm
			Real<FPP> norm() const {return mat.norm();}

			//!  \brief Normalize the matrix according to its Frobenius norm (norm_type==-2) by default.
			// \param norm_type: the type of norm used to normalize the matrix. -2 Frobenius norm, 1 1-norm, 2 2-norm (spectral norm), -1 inf-norm.
			void normalize(int norm_type=-2);


			//!	\param nbr_iter_max : maximum number of iteration for the power algo
			//! \param threshold : threshold until convergence
			//! \param flag : convergence flag
			//! \return Return the estimated spectral norm (maximum singular value in absolute value) using power iteration algorithm
			//! See also, template<typename FPP> FPP power_iteration(const MatDense<FPP,Cpu> & A, const faust_unsigned_int nbr_iter_max,FPP threshold,faust_int & flag);
			Real<FPP> spectralNorm(const faust_unsigned_int nbr_iter_max,FPP threshold, int & flag) const;

			//! \brief Compute the trace of the MatDense
			//! \return  the trace
			FPP trace() const {return mat.trace();}

			//! \brief Transpose the MatDense
			void transpose();

			//! \brief Replaces the matrix by its conjugate.
			void conjugate();
			void conjugate(const bool eval);

			//! \brief Replaces the matrix by its transconjugate.
			void adjoint();

			//! \brief Replace this by (this) * A
			void multiplyRight(MatDense<FPP,Cpu> const& A);

			void multiplyRight(MatSparse<FPP, Cpu> const& A);


			//!  \brief replace this by lambda * (*this)
			void scalarMultiply(FPP const lambda);

			//!  \brief replace this by lambda * (*this) using element by element multiplication
			void scalarMultiply(MatDense<FPP,Cpu> const& A);

			//! Replaces the matrix by its support (all nonzeros are set to ones).
			void nonzerosToOnes();

			//! \brief (*this) = (*this) + A
			void add(MatDense<FPP,Cpu> const& A);
			//
			//! \brief (*this) = (*this) + A
			void add(MatSparse<FPP,Cpu> const& A);

			//!  \brief (*this) = (*this) - A
			void sub(MatDense<FPP,Cpu> const& A);

			void sub(MatSparse<FPP,Cpu> const& A);

			//! \brief Displays the MatDense
			void Display() const;

			/**
			 * \brief Prints matrix with an optional header/title.
			 */
			void print_mat(std::string header="") const { if (header != "") std::cout << header << std::endl; std::cout << mat << std::endl;}


			//! \brief Returns all the features of the MatDense.
			std::string to_string(const bool transpose=false, const bool displaying_small_mat_elts=false) const;

			std::string to_string(int32_t nrows, int32_t ncols, bool transpose, Real<FPP> density, int32_t nnz, bool is_identity);

			//!  \brief Write MatDense into text file
			//! \param filename : name of the file
			//!
			//! The first line of the file contains 2 integer : the number of row and the number of column
			//! All the other line contains one coefficient in ColMajor access of the MatDense
			void print_file(const char* filename)const;

			matvar_t* toMatIOVar(bool transpose, bool conjugate, const char *var_name=nullptr) const;
			//!  \brief Creates a MatDense from at matio variable
			void from_matio_var(matvar_t* var);
			//!  \brief Creates a MatDense from at .mat file
			void read_from_mat_file(const char *filepath, const char *var_name);
			//!  \brief Saves a MatDense to a .mat file
			void save_to_mat_file(const char *filepath, const char *var_name);
			//! \brief Keeps only real part of the matrix, erase imaginary part.
			//
			//	If the matrix is real, it does nothing.
			void real();
			void real(MatDense<Real<FPP>, Cpu> & real_mat) const;
			/**
			 * Contrary to real() functions it returns a new MatDense of type Real<FPP2> but left this one untouched.
			 */
			template<typename FPP2>
				MatDense<Real<FPP2>, Cpu> to_real() const;
			template<typename FPP2>
				MatDense<FPP2, Cpu> cast() const;
			Real<FPP> normL1(const bool transpose=false) const;
			Real<FPP> normL1(faust_unsigned_int&, const bool transpose=false) const;
			Real<FPP> normInf(const bool transpose=false) const;
			Real<FPP> normInf(faust_unsigned_int&, const bool transpose=false) const;
			Vect<FPP,Cpu> get_col(faust_unsigned_int id) const;
			Vect<FPP,Cpu> get_row(faust_unsigned_int id) const;
			MatDense<FPP,Cpu>* get_cols(faust_unsigned_int start_col_id, faust_unsigned_int num_cols) const;
			MatDense<FPP,Cpu>* get_cols(const faust_unsigned_int* col_ids, faust_unsigned_int n) const;
			MatDense<FPP,Cpu>* get_cols(std::vector<int> col_ids) const;
			/** \brief Returns true if this[:,id_this] == other[:, id_other] at the specified precision. */
			bool eq_cols(const MatDense<FPP, Cpu> & other, faust_unsigned_int id_this, faust_unsigned_int id_other, const Real<FPP>& precision) const;
			/** \brief Returns true if this[id_this, :] == other[id_other, :] at the specified precision. */
			bool eq_rows(const MatDense<FPP, Cpu> & other, faust_unsigned_int id_this, faust_unsigned_int id_other, const Real<FPP>& precision) const;
			/** \brief Returns the sum of the column of index id. */
			FPP sum_col(faust_unsigned_int id) const;
			/** \brief Returns the sum of the row of index id. */
			FPP sum_row(faust_unsigned_int id) const;
			// \brief: sum all entries of this.
			FPP sum() const { return mat.sum();};
			void delete_col(int offset);
			void delete_row(int offset);

			MatDense<FPP,Cpu>* get_rows(faust_unsigned_int start_row_id, faust_unsigned_int num_rows) const;
			void get_rows(faust_unsigned_int start_row_id, faust_unsigned_int num_rows, MatDense<FPP, Cpu>& out_rows) const;
			void get_rows(const faust_unsigned_int* row_ids, faust_unsigned_int n, MatDense<FPP, Cpu>& out_rows) const;
			MatDense<FPP,Cpu>* get_rows(const faust_unsigned_int* row_ids, faust_unsigned_int n) const;

			MatDense<FPP,Cpu> get_block(faust_unsigned_int i, faust_unsigned_int j, faust_unsigned_int nrows, faust_unsigned_int ncols);

			void set_block(faust_unsigned_int i, faust_unsigned_int j, MatDense<FPP,Cpu> & block);

			FPP min_coeff() const;
			Vect<FPP,Cpu> rowwise_min() const;
			Vect<FPP,Cpu> rowwise_min(int* col_indices) const;
			Vect<FPP,Cpu> rowwise_max(int* col_indices) const;

			void operator=(MatDense<FPP,Cpu> const& A);

			template<typename FPP1>
				void operator=(MatDense<FPP1,Cpu> const& A);
			template<typename FPP1>
				void operator=(MatSparse<FPP1,Cpu> const& A){MatSparse<FPP,Cpu> AT(A);this->operator=(AT);};

			void operator=(MatSparse<FPP,Cpu> const& A);
			void operator-=(MatDense<FPP,Cpu> const& A){sub(A);}
			void operator-=(MatSparse<FPP,Cpu> const& A){sub(A);}
			void operator+=(MatDense<FPP,Cpu> const& A){add(A);}
			void operator*=(MatDense<FPP,Cpu> const& A){multiplyRight(A);}

			void operator*=(FPP lambda){scalarMultiply(lambda);}
			void operator/=(FPP lambda){scalarMultiply(1.0/lambda);}


			void swap_cols(const faust_unsigned_int id1, const faust_unsigned_int id2);
			void swap_rows(const faust_unsigned_int id1, const faust_unsigned_int id2);



			//        friend void gemm_core<>(const MatDense<FPP,Cpu> & A,const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB);
			//        friend void multiply<>(const MatDense<FPP,Cpu> & A,const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C);
			//        friend void spgemm<>(const MatSparse<FPP,Cpu> & A,const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char  typeA, char  typeB);
			//        friend void multiply<>(const Transform<FPP,Cpu> & A, const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, char typeA, char typeMult);
			//        friend void gemv<>(const MatDense<FPP,Cpu> & A,const Vect<FPP,Cpu> & x,Vect<FPP,Cpu> & y,const FPP & alpha, const FPP & beta, char typeA);
			// modif AL AL
			friend void gemm_core<>(const MatDense<FPP,Cpu> & A,const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB);
			friend void Faust::multiply<>(const MatDense<FPP,Cpu> & A,const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C);
			friend void spgemm<>(const MatSparse<FPP,Cpu> & A,const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char  typeA, char  typeB);
			friend void spgemm<>(const MatDense<FPP,Cpu> & A,const MatSparse<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char typeA, char typeB);
			friend void Faust::multiply<>(const Transform<FPP,Cpu> & A, const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, char typeA, char typeMult);
			friend void gemv<>(const MatDense<FPP,Cpu> & A,const Vect<FPP,Cpu> & x,Vect<FPP,Cpu> & y,const FPP & alpha, const FPP & beta, char typeA);
			//	friend void  MatSparse<FPP,Cpu>::multiply(MatDense<FPP,Cpu> & M,const char opThis) const;
			friend double Transform<FPP,Cpu>::normL1(const bool transpose, const bool full_array/*=true*/, const int batch_sz/*=1*/) const;

			// cf. faust_linear_algebra.h
			friend void gemm_gen<>(const MatGeneric<FPP, Cpu>& A, const MatGeneric<FPP, Cpu>& B, MatDense<FPP, Cpu>& out, const FPP alpha, const FPP beta, const char opA, const char opB);

			bool estNulle()const{return isZeros;}

			static MatDense<FPP,Cpu>* randMat(faust_unsigned_int num_rows, faust_unsigned_int num_cols);
			static MatDense<FPP,Cpu>* randMat(faust_unsigned_int num_rows, faust_unsigned_int num_cols, float density);
			//\param : per_row means the density applies for each line rather than globally for the matrix
			static MatDense<FPP,Cpu>* randMat(faust_unsigned_int num_rows, faust_unsigned_int num_cols, float density, bool per_row);
			/**
			 * \brief Returns the lower-triangular matrix with or without the diagonal.
			 *
			 * param diag true to include the diagonal, false otherwise.
			 */
			MatDense<FPP,Cpu> lower_tri(const bool diag=true) const;
			/**
			 * \brief Returns the upper-triangular matrix with or without the diagonal.
			 *
			 * param diag true to include the diagonal, false otherwise.
			 */
			MatDense<FPP,Cpu> upper_tri(const bool diag=true) const;

			std::vector<std::pair<int,int>> get_diag_indices(int index);
			std::vector<std::pair<int,int>> get_antidiag_indices(int index);
			Vect<FPP, Cpu> diagonal(int index);
			Vect<FPP, Cpu> adiagonal(int index);
			Vect<FPP, Cpu> gen_diagonal(int index, bool diag /* true for diagonal, false for anti-diagonal*/);

			size_t getNBytes() const;

			void copyBuf(FPP* dst_buf) const;

			/**
			 * \brief Returns the nonzeros indices.
			 */
			std::list<std::pair<int,int>> nonzeros_indices(const double& tol=0) const;

			/**
			 * \brief Returns the best low rank approximation this = bestX * bestY using the svd.
			 * \param svd_impl 0 for BDCSVD, 1 for JacobiSVD.
			 */
			template<typename SVDImpl>
				void best_low_rank(const int &r, MatDense<FPP,Cpu> &bestX, MatDense<FPP, Cpu> &bestY, SVDImpl svd) const;

			void best_low_rank(const int &r, MatDense<FPP,Cpu> &bestX, MatDense<FPP, Cpu> &bestY) const;

			void best_low_rank2(MatDense<FPP,Cpu> &bestX, MatDense<FPP, Cpu> &bestY) const;
			void best_low_rank3(MatDense<FPP,Cpu> &bestX, MatDense<FPP, Cpu> &bestY) const;

			void solveNullspaceQR(MatDense<FPP, Cpu>& X) const;

			void initJacobiSVD(Eigen::JacobiSVD<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>>& svd);
			/**
			 * \brief Returns the best rank-1 approximate this = bestX * bestY using the svd/eigensolver/power iteration.
			 */
			void approx_rank1(MatDense<FPP,Cpu> &bestX, MatDense<FPP, Cpu> &bestY) const;
			/**
			 * Returns the vector of the indices of the nonzeros of the row of index row_id.
			 */
			std::vector<int> row_nonzero_inds(faust_unsigned_int row_id) const;
			/**
			 * Returns the vector of the indices of the nonzeros of the column of index col_id.
			 */
			std::vector<int> col_nonzero_inds(faust_unsigned_int col_id) const;
			/**
			 * Returns in submat the submatrix defined by the row indices row_ids and the column indices col_ids.
			 */
			void submatrix(const std::vector<int> &row_ids, const std::vector<int> &col_ids, MatDense<FPP, Cpu> & submat) const;

			void submatrix(const std::vector<int> &row_ids, const std::vector<int> &col_ids, FPP* submat_data) const;

			/**
			 * Assigns this[row_ids[i], col_id] to values[i, val_col_id] for each i in {0, ..., row_ids.size()}.
			 */
			void set_col_coeffs(faust_unsigned_int col_id, const std::vector<int> &row_ids, const MatDense<FPP, Cpu> &values, faust_unsigned_int val_col_id);
			void set_col_coeffs(faust_unsigned_int col_id, const std::vector<int> &row_ids, const FPP* values, faust_unsigned_int val_col_id, faust_unsigned_int values_nrows);

			/**
			 * Assigns this[row_id, col_ids[j]] to values[val_row_id, j] for each j in {0, ..., col_ids.size()}.
			 */
			void set_row_coeffs(faust_unsigned_int row_id, const std::vector<int> &col_ids, const MatDense<FPP, Cpu> &values, faust_unsigned_int val_row_id);
			void set_row_coeffs(faust_unsigned_int row_id, const std::vector<int> &col_ids, const FPP* values, faust_unsigned_int val_row_id, faust_unsigned_int values_nrows);


			bool containsNaN() const;

			/**
			 * Ensures MatGeneric dims are the same as Eigen matrix.
			 */
			void update_dims() { this->dim1 = mat.rows(); this->dim2 = mat.cols();}

			private:
			Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic> mat;
			bool isZeros;
			static const char * m_className;

#ifdef __COMPILE_TIMERS__
			public:
			Timer t_local_muliplyLeft;
			//temporary members
			static Timer t_constr;
			static Timer t_get_coeff;
			static Timer t_get_coeffs;
			static Timer t_set_coeff;
			static Timer t_set_coeffs;
			static Timer t_set_coeffs2;
			static Timer t_resize;
			static Timer t_check_dim;
			static Timer t_max;
			static Timer t_transpose;
			static Timer t_mult_right;
			static Timer t_mult_left;
			static Timer t_scalar_multiply;
			static Timer t_add;
			static Timer t_sub;
			static Timer t_print_file;
			static Timer t_spectral_norm;
			static Timer t_spectral_norm2;

			static Timer t_power_iteration;
			static Timer t_multiply;
			static Timer t_gemm;
			static Timer t_add_ext;

			void print_timers()const;
#endif
			/**
			 * \brief index this(row_ids, col_ids) Eigen matrix and mutliplies the result by in_mat into out_mat (two Eigen dense matrices of respective types MatType1, MatTyp2)
			 *
			 */
			template<typename MatType1, typename MatType2>
				void eigenIndexMul(const faust_unsigned_int* row_ids, const faust_unsigned_int* col_ids, size_t nrows, size_t ncols, const MatType1 &in_mat, MatType2 &out_mat, bool transpose = false, bool conjugate = false);



			MatDense<FPP, Cpu> to_dense() const;
		};

	// \brief Computes the Kronecker product A \otimes B into out.
	template<typename FPP>
			void kron(const MatDense<FPP,Cpu> &A, const MatDense<FPP, Cpu> &B, MatDense<FPP,Cpu>& out);



}

#include "faust_MatDense.hpp"


#endif
