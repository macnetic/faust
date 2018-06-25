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
#ifndef __FAUST_CONSTRAINT_TYPE_H__
#define __FAUST_CONSTRAINT_TYPE_H__


template<typename FPP,Device DEVICE> class ConstraintInt;
template<typename FPP,Device DEVICE> class ConstraintFPP;
template<typename FPP,Device DEVICE> class ConstraintMat;


//! \struct ConstraintType
//! \brief This structure defined the type of constraint. See following table for more precision about the type of constraint. <br>
 //! <img src="../../doc/html/constraint.png" alt="constraint parameters" width=800px />
template<typename FPP,Device DEVICE>
struct ConstraintType
{

	typedef const  Faust::ConstraintFPP<FPP,DEVICE>   ConstraintTypeNormcol;
	typedef const  Faust::ConstraintFPP<FPP,DEVICE>   ConstraintTypeNormlin;
	typedef const  Faust::ConstraintMat<FPP,DEVICE>   ConstraintTypeSupp;
	typedef const  Faust::ConstraintMat<FPP,DEVICE>   ConstraintTypeConst;
	typedef const  Faust::ConstraintInt<FPP,DEVICE>   ConstraintTypeSp;
	typedef const  Faust::ConstraintInt<FPP,DEVICE>   ConstraintTypeSpcol;
	typedef const  Faust::ConstraintInt<FPP,DEVICE>   ConstraintTypeSplin;
	typedef const  Faust::ConstraintInt<FPP,DEVICE>   ConstraintTypeSplincol;
	typedef const  Faust::ConstraintInt<FPP,DEVICE>   ConstraintTypeSpPos;
	typedef const  Faust::ConstraintInt<FPP,DEVICE>   ConstraintTypeBlkdiag;

};


#endif
