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
#ifndef __FAUST_MAT_GENERIC_H__
#define __FAUST_MAT_GENERIC_H__

#include "faust_constant.h"
#include "faust_LinearOperator.h"
#include "faust_Vect.h"
#include "faust_MatDense.h"
#include "matio.h"
#include "Eigen/Core"
#include <list>
#include <utility>
/**
 * \class MatGeneric faust_MatGeneric.h
 * \brief This MatGeneric class serves as a base class for the derived class Faust::MatDense and Faust::MatSparse .
 * 	Member variable dim1 and dim2 correspond to the dimension of the matrix
 *
*/


//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{


    template<typename FPP,FDevice DEVICE>
    class LinearOperator;

    template<typename FPP,FDevice DEVICE>
    class MatBSR;

    template<typename FPP,FDevice DEVICE>
    class MatGeneric : public Faust::LinearOperator<FPP,DEVICE>
	{

		public:


			MatGeneric() : dim1(0), dim2(0), is_ortho(false), is_identity(false) {}

			MatGeneric(faust_unsigned_int dim1_, faust_unsigned_int dim2_) : dim1(dim1_), dim2(dim2_), is_ortho(false), is_identity(false) {}

			MatGeneric(faust_unsigned_int dim1_, faust_unsigned_int dim2_, bool is_ortho, bool is_identity) : dim1(dim1_), dim2(dim2_), is_ortho(is_ortho), is_identity(is_identity) {}


			//! \brief get the dimension of op((*this))
			//! \param op : char if op=='N' op(*this)=(*this), if op=='T' op((*this))==transpose((*this))
			//! \param nbRowOp : (in/out) faust_unsigned_int, number of rows of op(*this)
			//! \param nbColOp : (in/out) faust_unsigned_int, number of columns of op(*this)
			void setOp(const char op, faust_unsigned_int& nbRowOp, faust_unsigned_int& nbColOp)const;

			//! \brief return the number of rows of (*this)
			faust_unsigned_int getNbRow() const {return dim1;}

			//! \brief return the number of column of (*this)
			faust_unsigned_int getNbCol() const {return dim2;}

			//! \brief resize (*this)
			//! \param dim1_ : new number of rows
			//! \param dim2_ : new number of columns
			void resize(const faust_unsigned_int dim1_,const faust_unsigned_int dim2_){dim1=dim1_;dim2=dim2_;}

			//purely virtual function : must be redefined in all the descendant class
			/*!  \brief return a "copy" to the given matrix
			 *  \param isOptimize (optionnal) : boolean which the style of copy <br>
			 -True, the return copy is optimized for the product <br>
			 which means dynamic type of the copy could be different from the original one <br>
			 -False, the return copy is simple, the dynamic type isn't changed <br>   						(default value is False) <br>
			//! \return  a pointer of MatGeneric
			//  \warning the dynamic type of the copy can be different from the original object,
			//  \warning dynamic allocation of the return pointer, must be delete by hand
			*/
			virtual MatGeneric<FPP,DEVICE>* Clone(const bool isOptimize=false) const=0;

			//! \brief compute MatGeneric-vector multiplication
			//! \param vec : the vector
			//! \param opThis : character
			//! vec = (*this) * vec if opThis='N'
			// vec = (*this)' * vec if opThis='T'
			virtual void multiply(Faust::Vect<FPP,DEVICE> & vec, char opThis='N') const=0;


			//! \brief compute MatGeneric-MatDense multiplication
			//! \param M : the dense matrix
			//! \param opThis : character
			//! M = (*this) * M if opThis='N'
			// M = (*this)' * M if opThis='T'
			virtual void multiply(Faust::MatDense<FPP,DEVICE> & M, char opThis) const=0;

			virtual void multiply(Faust::MatSparse<FPP, DEVICE>& M, char opThis) const=0;

			//! \brief Replace this by (this) * A
			virtual void multiplyRight(Faust::MatSparse<FPP, DEVICE> const& M) =0;

			//! \brief transpose the matrix
			virtual void transpose()=0;

			//! \brief Replaces the matrix by its conjugate.
			virtual void conjugate(const bool eval=true)=0;

			virtual void adjoint()=0;

			//! \brief return the number of non-zeros element in the matrix
			virtual faust_unsigned_int getNonZeros()const=0;

			virtual size_t getNBytes() const=0;

			//!brief return the percentage of non-zeros coefficient in the matrix,
			//! \return value between 0 and 1
			float density() const{return ((float) this->getNonZeros())/((float)this->getNbCol()*this->getNbRow());}

			//! \brief get the dynamic type of the matrix (SPARSE, DENSE, etc.)
			virtual MatType getType() const=0;

			//! \brief multiply a matrix by the given scalar alpha
			// \param alpha : multplicative scalar
			virtual void operator*=(const FPP alpha)=0;


			//! \brief Display the features of the matrix (type Dense/Sparse, size, nnz, density of nnz ... )
			virtual void Display() const;

			//! \brief Returns the features of the matrix (type Dense/Sparse, size, nnz, density of nnz ... )
			virtual std::string to_string(MatType type, const bool transpose=false, const bool displaying_small_mat_elts=false) const;
			//
			//! \brief Returns the features of the matrix (type Dense/Sparse, size, nnz, density of nnz ... )
			virtual std::string to_string(const bool transpose=false, const bool displaying_small_mat_elts=false) const;

			static std::string to_string(int32_t nrows, int32_t ncols, bool transpose, Real<FPP> density, int32_t nnz, bool is_identity, MatType type);

			static std::string get_scalar_type_str();

			//! \brief Converts the Matrix to a matio variable, especially useful for writing into a file with libmatio.
			// \param transpose: set to true to obtain the matio variable for the transpose Matrix.
			// \param conjugate: set it to true to obtain the matio variable for the conjugate Matrix.
			// \return The matio variable matvar_t if it succeeded or NULL otherwise.
			// \see Faust::Transform::save_mat_file()
			virtual matvar_t* toMatIOVar(bool transpose, bool conjugate) const=0;

			//TODO: only one function norm() with argument for type of norm
			//! \brief Computes the L1-norm of the matrix.
			// \param transpose: to compute the norm of the transpose matrix.
			// \see Faust::normL1(faust_unsigned_int&)
			// \see http://mathworld.wolfram.com/L1-Norm.html
			virtual Real<FPP> normL1(const bool transpose) const=0;

			//! \brief Frobenius norm.
			virtual Real<FPP> norm() const=0;

			//! \brief Computes the L1-norm of the matrix.
			// \param col_id: reference to receive the column index which the L1-norm is equal to the matrix's norm (if several exist, then the greater colummn index is kept).
			// \param transpose: to compute the norm of the transpose matrix.
			// \return The norm (its type is the matrix scalar's).
			// \see Faust::normL1()
			// \see http://mathworld.wolfram.com/L1-Norm.html
			virtual Real<FPP> normL1(faust_unsigned_int& col_id, const bool transpose) const=0;
			//
			//! \brief Returns a column of the matrix as a new Faust::Vect.
			virtual Faust::Vect<FPP,DEVICE> get_col(faust_unsigned_int id) const=0;

			//! \brief Returns a sub-group of columns of this matrix as the same type of matrix
			virtual Faust::MatGeneric<FPP,DEVICE>* get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const=0;
			//! \brief Returns a sub-group of rows of this matrix as the same type of matrix
			virtual Faust::MatGeneric<FPP,DEVICE>* get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const=0;
			//! \brief Returns a sub-group of columns of this matrix as the same type of matrix
			virtual Faust::MatGeneric<FPP,DEVICE>* get_cols(const faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const=0;
			//! \brief Returns a sub-group of rows of this matrix as the same type of matrix
			virtual Faust::MatGeneric<FPP,DEVICE>* get_rows(const faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const=0;
			virtual std::list<std::pair<int,int>> nonzeros_indices() const=0;
			void set_orthogonal(const bool is_ortho) { this->is_ortho = is_ortho; /* TODO: move def in hpp*/}
			void set_id(const bool is_identity) { this->is_identity = is_identity; /* TODO: move def in hpp*/}

			virtual void setZeros()=0;
			virtual bool containsNaN()=0;
			virtual const FPP& operator()(faust_unsigned_int i, faust_unsigned_int j)const =0;

			bool is_orthogonal() { return this->is_ortho; /* TODO: move def in hpp*/}
			bool is_id() const { return this->is_identity; /* TODO: move def in hpp*/}


			//! \brief
			//! \warning : declare a virtual destructor is mandatory for an abstract class
			//! in order to allow child class destructor to clean up in case of pointer to the abstract class
			virtual ~ MatGeneric()=0;




		protected:
			faust_unsigned_int dim1;
			faust_unsigned_int dim2;
			bool is_ortho;
			bool is_identity;
			//TODO: add is_zeros



	};

    //!
    //! \brief compare which format between the sparse matrix and the dense matrix is the quickiest for multiplication with vector and return a pointer to the mother class MatGeneric, with a dynamic type equals to the most efficient format
   //! for multiplication
   //! \tparam M : MatDense
   //! \tparam S : MatSparse
   //  \return a pointer of MatGeneric

   //template <typename FPP, FDevice DEVICE>
   template<typename FPP, FDevice DEV>
   Faust::MatGeneric<FPP,DEV>* optimize(Faust::MatDense<FPP,DEV> const & M,Faust::MatSparse<FPP,DEV> const & S);


}
#include "faust_MatGeneric.hpp"

#endif
