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

#ifndef __FAUST_TOOLS_MEX_H__
#define __FAUST_TOOLS_MEX_H__

#include "mex.h"
#include <vector>
#include "faust_constant.h"

namespace Faust {
	template<typename FPP, Device DEVICE> class ConstraintGeneric;
	template<typename FPP, Device DEVICE> class Vect;
	template<typename FPP, Device DEVICE> class Params;
	template<typename FPP, Device DEVICE> class MatGeneric;
	template<typename FPP, Device DEVICE> class MatDense;
	template<typename FPP, Device DEVICE> class MatSparse;
}

/*!
*  \brief convert the matlab mxArray* into a Faust::Vect<FPP,Cpu>, no shared memory
*  \param vec_array : pointer to the mxArray* (matlab format) representing a dense column Vector 
*  \tparam vec : Faust::Vect<FPP,Cpu>  
        */
template<typename FPP>
void getFaustVec(const mxArray * vec_array,Faust::Vect<FPP,Cpu> & vec);

/*!
*  \brief convert the matlab mxArray* into a Faust::MatDense<FPP,Cpu>, no shared memory
*  \param Mat_array : pointer to the mxArray* (matlab format) representing a dense matrix 
*  \tparam Mat : Faust::MatDense<FPP,Cpu>
        */
template<typename FPP>
void getFaustMat(const mxArray* Mat_array,Faust::MatDense<FPP,Cpu> & Mat);

/*!
*  \brief convert the matlab mxArray* into a Faust::Vect<FPP,Cpu>, no shared memory
*  \param spMat_array : pointer to the mxArray* (matlab format) representing a sparse matrix  
*  \tparam S : Faust::MatSparse<FPP,Cpu>  
        */
template<typename FPP>
void getFaustspMat(const mxArray* spMat_array,Faust::MatSparse<FPP,Cpu> & S);

/*!
*  \brief return a matlab mxArray* dense matrix from a Faust::MatDense<FPP,Cpu>, no shared memory
*  \param[out] mxArray* : pointer to the mxArray* (matlab format) representing a dense matrix  
*  \tparam[in] S : Faust::MatDense<FPP,Cpu>  
*/
template<typename FPP>
mxArray*  FaustMat2mxArray(const Faust::MatDense<FPP,Cpu>& M);

/*!
*  \brief return a matlab mxArray** representing a cell-array of matlab matrix from a std::vector<Faust::MatDense<FPP,Cpu> >, no shared memory
*  \param[out] cellfacts : matlab cell-array of matlab dense matrix 
*  \tparam[in] facts : std::vector<Faust::MatDense<FPP,Cpu> > list of Dense matrix
*/
template<typename FPP>
void setCellFacts(mxArray ** cellFacts,std::vector<Faust::MatDense<FPP,Cpu> > & facts);


/*!
*  \brief convert a matlab mxArray* representing a cell-array of palm4MSA constraint into a list of C++ Faust::ConstraintGeneric<FPP,Cpu>
*  \param[in] mxConsS : matlab cell-array representing the list of constraint of palm4MSA algorithm
*                      each cell C is also 4by1 cell-array : - C{1} is the chararcter defining the type of constraint
*                                                            - C{2} is the number of row of the constraint
*							   *							     - C{3} is the number of columns of the constraint
*                                                            - C{4} is the parameter of the constraint		       
*  \tparam[in] concsS : std::vector of Faust::ConstraintGeneric<FPP,Cpu> C++ equivalent of the matlab palm4MSA constraint
*/
template<typename FPP>
void getConstraint(std::vector<const Faust::ConstraintGeneric<FPP,Cpu>*> & consS,mxArray* mxCons);



/*!
        *  \brief convert the matlab mxMat* into a Faust::MatSparse<FPP,Cpu> and push it to the end of the 			list vec_spmat
        *  \param * mxMat : pointer to the matlab matrix format
	*  \tparam vec_spmat : std::vector of Faust::MatSparse<FPP,Cpu>
        */
template<typename FPP>
void setVectorFaustMat(std::vector<Faust::MatDense<FPP,Cpu> > &vecMat, mxArray *Cells);


/*!
        *  \brief convert the matlab mxMat* into a Faust::MatGeneric<FPP,Cpu>* and push it to the end of list_mat
        *  \param * mxMat : pointer to the matlab matrix format
	*  \tparam list_mat : std::vector of pointer of Faust::MatGeneric<FPP,Cpu>
	   \warning : dynamic allocation is made so , "delete list_mat[list_mat.size()-1]" must be done	
        */
template<typename FPP>
void concatMatGeneric(const mxArray * mxMat,std::vector<Faust::MatGeneric<FPP,Cpu> *> &list_mat);

void testCoherence(const mxArray* params,std::vector<bool> & presentFields);
template<typename FPP>
void DisplayParams(const Faust::Params<FPP,Cpu> & params);

#include "tools_mex.hpp"

#endif
