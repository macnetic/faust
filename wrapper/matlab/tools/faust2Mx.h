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

#ifndef __FAUST_FAUST2MX_H__
#define __FAUST_FAUST2MX_H__

#include "mex.h"
#include <vector>
#include "faust_constant.h"
#include <complex>


namespace Faust {
	template<typename FPP, Device DEVICE> class ConstraintGeneric;
	template<typename FPP, Device DEVICE> class Vect;
	template<typename FPP, Device DEVICE> class Params;
	template<typename FPP, Device DEVICE> class MatGeneric;
	template<typename FPP, Device DEVICE> class MatDense;
	template<typename FPP, Device DEVICE> class MatSparse;
	template<typename FPP, Device DEVICE> class Transform;
	template<typename FPP, Device DEVICE> class LinearOperator;
}




/*!
*  \brief return a matlab mxArray* dense matrix from a Faust::MatVec<FPP,Cpu>, no shared memory
*  \param[out] mxArray* : pointer to the mxArray* (matlab format) representing a dense matrix  
*  \tparam[in] M : Faust::Vec<FPP,Cpu>  
*/
template<typename FPP>
mxArray*  FaustVec2mxArray(const Faust::Vect<FPP,Cpu>& M);
template<typename FPP>
mxArray*  FaustVec2mxArray(const Faust::Vect<std::complex<FPP>,Cpu>& v);




/*!
*  \brief return a matlab mxArray* dense matrix from a Faust::MatDense<FPP,Cpu>, no shared memory
*  \param[out] mxArray* : pointer to the mxArray* (matlab format) representing a dense matrix  
*  \tparam[in] M : Faust::MatDense<FPP,Cpu>  
*/
template<typename FPP>
mxArray*  FaustMat2mxArray(const Faust::MatDense<FPP,Cpu>& M);
template<typename FPP>
mxArray*  FaustMat2mxArray(const Faust::MatDense<std::complex<FPP>,Cpu>& M);




// split complex value ptr into real value ptr and imaginary value ptr,
//real_ptr and imag_ptr must be allocated for nb_element
template<typename FPP>
void splitComplexPtr(const std::complex<FPP>*  cpx_ptr, int nb_element, FPP* & real_ptr, FPP* & imag_ptr);
template<typename FPP>
void mxArray2Ptr(const mxArray* mxMat, std::complex<FPP>* & ptr_data);


/*!
*  \brief return a matlab mxArray** representing a cell-array of matlab matrix from a std::vector<Faust::MatDense<FPP,Cpu> >, no shared memory
*  \param[out] cellfacts : matlab cell-array of matlab dense matrix 
*  \tparam[in] facts : std::vector<Faust::MatDense<FPP,Cpu> > list of Dense matrix
*/
template<typename FPP>
void setCellFacts(mxArray ** cellFacts,std::vector<Faust::MatDense<FPP,Cpu> > & facts);


template<typename FPP>
void DisplayParams(const Faust::Params<FPP,Cpu> & params);


#include "faust2Mx.hpp"

#endif
