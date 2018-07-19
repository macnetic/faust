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
/****************************************************************************/
#ifndef __FAUST_MAT_GENERIC_H__
#define __FAUST_MAT_GENERIC_H__

#include "faust_constant.h"
#include "faust_LinearOperator.h"
#include "faust_Vect.h"
#include "faust_MatDense.h"
#include "matio.h"
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


    template<typename FPP,Device DEVICE>
    class LinearOperator;

    template<typename FPP,Device DEVICE>
    class MatGeneric : public Faust::LinearOperator<FPP,DEVICE>
    {
        public:
        
	MatGeneric() : dim1(0), dim2(0) {}

        MatGeneric(faust_unsigned_int dim1_, faust_unsigned_int dim2_) : dim1(dim1_), dim2(dim2_){}
    	
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


	
	//! \brief transpose the matrix
	virtual void transpose()=0;

	//! \brief Replaces the matrix by its conjugate.
	virtual void conjugate()=0;
	
	//! \brief return the number of non-zeros element in the matrix
	virtual faust_unsigned_int getNonZeros()const=0;

	//!brief return the percentage of non-zeros coefficient in the matrix, 
	//! \return value between 0 and 1
	float density() const{return ((float) this->getNonZeros())/((float)this->getNbCol()*this->getNbRow());}
	
	//! \brief get the dynamic type of the matrix (SPARSE or DENSE)
	virtual MatType getType() const=0;
	
	//! \brief multiply a matrix by the given scalar alpha
	// \param alpha : multplicative scalar
	virtual void operator*=(const FPP alpha)=0;


	//! \brief Display the features of the matrix (type Dense/Sparse, size, nnz, density of nnz ... )
	virtual void Display() const;

	//! \brief Returns the features of the matrix (type Dense/Sparse, size, nnz, density of nnz ... )
	virtual std::string to_string(const bool transpose=false) const;

	//! \brief Converts the Matrix to a matio variable, especially useful for writing into a file with libmatio.
	// \param transpose: set to true to obtain the matio variable for the transpose Matrix.
	// \param conjugate: set it to true to obtain the matio variable for the conjugate Matrix. 
	// \return The matio variable matvar_t if it succeeded or NULL otherwise.
	// \see Faust::Transform::save_mat_file()
	virtual matvar_t* toMatIOVar(bool transpose, bool conjugate) const=0;

	//! \brief Computes the L1-norm of the matrix.
	// \see Faust::normL1(faust_unsigned_int&)
	// \see http://mathworld.wolfram.com/L1-Norm.html
	virtual FPP normL1() const=0;

	//! \brief Computes the L1-norm of the matrix.
	// \param col_id: reference to receive the column index which the L1-norm is equal to the matrix's norm (if several exist, then the greater colummn index is kept).
	// \return The norm (its type is the matrix scalar's).
	// \see Faust::normL1()
	// \see http://mathworld.wolfram.com/L1-Norm.html
	virtual FPP normL1(faust_unsigned_int& col_id) const=0;
	//
	//! \brief Returns a column of the matrix as a new Faust::Vect.
	virtual Faust::Vect<FPP,DEVICE> get_col(faust_unsigned_int id) const=0;

	//! \brief 
	//! \warning : declare a virtual destructor is mandatory for an abstract class
	//! in order to allow descendant class destructor to clean up in case of pointer to the abstract class 
	virtual ~ MatGeneric()=0;
	

	
	
        protected:
        faust_unsigned_int dim1;
        faust_unsigned_int dim2;


		    
    };

    //! 	
    //! \brief compare which format between the sparse matrix and the dense matrix is the quickiest for multiplication with vector and return a pointer to the mother class MatGeneric, with a dynamic type equals to the most efficient format
   //! for multiplication
   //! \tparam M : MatDense
   //! \tparam S : MatSparse
   //  \return a pointer of MatGeneric
   
   //template <typename FPP, Device DEVICE>
   template<typename FPP>	 
   Faust::MatGeneric<FPP,Cpu>* optimize(Faust::MatDense<FPP,Cpu> const & M,Faust::MatSparse<FPP,Cpu> const & S);
	

}
#include "faust_MatGeneric.hpp"

#endif
