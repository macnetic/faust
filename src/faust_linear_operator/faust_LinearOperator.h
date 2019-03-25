/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2018):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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
#ifndef __FAUST_LINEAR_OPERATOR_H__
#define __FAUST_LINEAR_OPERATOR_H__

#include "faust_constant.h"


// modif AL AL
//template<typename T,Device DEVICE> class MatDense;
//template<typename FPP,Device DEVICE> class MatDense;


/**
 * \class Faust::LinearOperator faust_LinearOperator.h
 * \brief This virtual template class represent all the different type of matrix that are used<br> (Dense matrix (Faust::MatDense), Sparse matrix (Faust::MatSparse),FAUST (Faust::Transform)). <br>
 *
 * \tparam T scalar numeric type, e.g float or double
*/


//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{
    // modif AL AL
    template<typename FPP,Device DEVICE>
    class MatDense;
	template<typename FPP, Device DEVICE>
		class Vect;

    template<typename FPP,Device DEVICE>
    class LinearOperator
    {
        public:

        virtual faust_unsigned_int getNbRow() const=0;
        virtual faust_unsigned_int getNbCol() const=0;
		virtual	void setOp(const char op, faust_unsigned_int& nbRowOp, faust_unsigned_int& nbColOp)const=0;

        /*!
        *  \brief Transpose the Faust::MatDense
        */
        virtual void transpose()=0;

		virtual Faust::Vect<FPP,DEVICE> multiply(const Faust::Vect<FPP,DEVICE> &v) const=0;

        /*! \brief faust_gemm — Performs a matrix calculation of the form alpha . (*this).B +beta.C, with B and C 2 faust_matrix and two real constants you supply. <br>

        * C = alpha .op(this).op(B) + beta.C; <br>
        * op(A) = A if typeA='N', op(A) = transpose(A) if typeA='T'<br>
        * op(B) = B if typeB='N', op(B) = transpose(B) if typeB='T'<br>
        * The Object C must be different of B. <br>
        * \param B : Dense Matrix
        * \param C : Dense Matrix
        * \param alpha : Template scalar coefficient
        * \param beta : Template scalar coefficient
        */
        virtual void faust_gemm(const Faust::MatDense<FPP,DEVICE> & B, Faust::MatDense<FPP,DEVICE> & C,const FPP & alpha, const FPP & beta, char typeA, char typeB)const=0;


	/** \brief 
	* check if the LinearOperator has real scalar or complex scalar 
        * */
	bool isReal() const;



    };
}

#include "faust_LinearOperator.hpp"

#endif
