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
#ifndef __FAUST_MAT_GENERIC_H__
#define __FAUST_MAT_GENERIC_H__

#include "faust_constant.h"
#include "faust_LinearOperator.h"
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

    // modif AL AL
    template<typename FPP,Device DEVICE>
    class LinearOperator;

    template<typename FPP,Device DEVICE>
    class MatGeneric : public Faust::LinearOperator<FPP,DEVICE>
    {
        public:
        MatGeneric() : dim1(0), dim2(0) {}
        MatGeneric(faust_unsigned_int dim1_, faust_unsigned_int dim2_) : dim1(dim1_), dim2(dim2_){}
    	
	void setOp(const char op, faust_unsigned_int& nbRowOp, faust_unsigned_int& nbColOp)const; 
        faust_unsigned_int getNbRow() const {return dim1;}
        faust_unsigned_int getNbCol() const {return dim2;}

        void resize(const faust_unsigned_int dim1_,const faust_unsigned_int dim2_){dim1=dim1_;dim2=dim2_;}

        protected:
        faust_unsigned_int dim1;
        faust_unsigned_int dim2;
    };

}
#include "faust_MatGeneric.hpp"

#endif
