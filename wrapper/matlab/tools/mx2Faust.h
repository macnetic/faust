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

#ifndef __FAUST_MX2FAUST_H__
#define __FAUST_MX2FAUST_H__

#include "mex.h"
#include <vector>
#include "faust_constant.h"
#include "faust_Params.h"
#include <complex>
#include <string>

namespace Faust {
	class ConstraintGeneric;
	template<typename FPP, Device DEVICE> class Vect;
	template<typename FPP, Device DEVICE, typename FPP2> class Params;
	template<typename FPP, Device DEVICE> class MatGeneric;
	template<typename FPP, Device DEVICE> class MatDense;
	template<typename FPP, Device DEVICE> class MatSparse;
	template<typename FPP, Device DEVICE> class Transform;
	template<typename FPP, Device DEVICE> class LinearOperator;
}

using namespace Faust;

/*!
*  \brief check if the Faust::Transform T has compatible scalar with MATLAB matrix Matlab_Mat (currently real is only compatible with real and complex is only compatible with complex)
*  \param T :  Faust::Transform<FPP,Cpu>
*  \tparam Matlab_Mat : mxArray pointer 
        */
template<typename FPP>
bool isScalarCompatible(const Faust::LinearOperator<FPP,Cpu> & L,const mxArray * Matlab_Mat);






/*!
*  \brief convert the matlab mxArray* into a Faust::MatDense<FPP,Cpu>, no shared memory
*  \param Mat_array : pointer to the mxArray* (matlab format) representing a dense matrix 
*  \tparam Mat : Faust::MatDense<FPP,Cpu>
        */
template<typename FPP>
void mxArray2FaustMat(const mxArray* Mat_array,Faust::MatDense<FPP,Cpu> & Mat);
/*template<typename FPP>
void mxArray2FaustMat(const mxArray* Mat_array,Faust::MatDense<std::complex<FPP>,Cpu> & Mat);*/



/*!
*  \brief convert the matlab mxArray* into a Faust::Vect<FPP,Cpu>, no shared memory
*  \param spMat_array : pointer to the mxArray* (matlab format) representing a sparse matrix  
*  \tparam S : Faust::MatSparse<FPP,Cpu>  
        */
template<typename FPP>
void mxArray2FaustspMat(const mxArray* spMat_array,Faust::MatSparse<FPP,Cpu> & S);


/* !
  \brief convert the data (real part or imaginary part) of a mxArray* into the pointer ptr_data
  \param mxMat : pointer to the mxArray (MATLAB matrix)
  \tparam ptr_data : pointer of template FPP where the data of the mxArray will be stored,
  \tparam mxGetDataFunc : functor to get the data of MATLAB matrix
			  possible value : - mxGetData  (get the real scalar part of the data) 
					                (cf https://fr.mathworks.com/help/matlab/apiref/mxgetdata.html)
			                   - mxGetImagData (get the imaginary scalar part of the data)
							(cf https://fr.mathworks.com/help/matlab/apiref/mxgetimagdata.html)  			  		
*/
template<typename FPP,class FUNCTOR>
void mxArray2PtrBase(const mxArray* mxMat,FPP* & ptr_data,FUNCTOR & mxGetDataFunc);

/* \brief convert the real scalar data  of a mxArray* into the pointer ptr_data 
*/
template<typename FPP>
void mxArray2Ptr(const mxArray* mxMat, FPP* & ptr_data);

/* \brief convert the complex scalar data  of a mxArray* into the pointer ptr_data 
*/
template<typename FPP>
void mxArray2Ptr(const mxArray* mxMat, std::complex<FPP>* & ptr_data);


/*!
        *  \brief convert the matlab mxMat* into a Faust::MatGeneric<FPP,Cpu>* and push it to the end of list_mat
        *  \param * mxMat : pointer to the matlab matrix format
	*  \tparam list_mat : std::vector of pointer of Faust::MatGeneric<FPP,Cpu>
	   \warning : dynamic allocation is made so , "delete list_mat[list_mat.size()-1]" must be done	
        */
template<typename FPP>
void concatMatGeneric(const mxArray * mxMat,std::vector<Faust::MatGeneric<FPP,Cpu> *> &list_mat);


/*!
*  \brief convert a matlab mxArray* representing a cell-array of palm4MSA constraint into a list of C++ Faust::ConstraintGeneric<FPP,Cpu>
*  \param[in] mxConsS : matlab cell-array representing the list of constraint of palm4MSA algorithm
*                      each cell C is also 4by1 cell-array : - C{1} is the chararcter defining the type of constraint
*                                                            - C{2} is the number of row of the constraint
*							   *							     - C{3} is the number of columns of the constraint
*                                                            - C{4} is the parameter of the constraint		       
*  \tparam[in] concsS : std::vector of Faust::ConstraintGeneric<FPP,Cpu> C++ equivalent of the matlab palm4MSA constraint
*/
template<typename FPP, typename FPP2>
void getConstraint(std::vector<const Faust::ConstraintGeneric*> & consS,mxArray* mxCons);



/*!
        *  \brief convert the matlab mxMat* into a Faust::MatSparse<FPP,Cpu> and push it to the end of the 			list vec_spmat
        *  \param * mxMat : pointer to the matlab matrix format
	*  \tparam vec_spmat : std::vector of Faust::MatSparse<FPP,Cpu>
        */
template<typename FPP>
void setVectorFaustMat(std::vector<Faust::MatDense<FPP,Cpu> > &vecMat, mxArray *Cells);


enum MAT_FIELD_TYPE
{
	/** All the matlab fields possible in params struct  (1st arg. in testCoherence() */
	NROW,
	NCOL,
	NFACTS,
	CONS,
	NITER1,
	NITER2,
	VERBOSE,
	FACT_SIDE,
	UPDATE_WAY,
	INIT_LAMBDA,
	SC_IS_CRITERION_ERROR,
	SC_ERROR_TRESHOLD,
	SC_MAX_NUM_ITS,
	SC_IS_CRITERION_ERROR2,
	SC_ERROR_TRESHOLD2,
	SC_MAX_NUM_ITS2,
	INIT_FACTS,
	INIT_D //only for FactHierarchicalF(G)FT
};// if you wan to add a field, dont forget to update mat_field_type2str() and MAT_FIELD_TYPE_LEN

const string mat_field_type2str(MAT_FIELD_TYPE f);
const MAT_FIELD_TYPE mat_field_str2type(const string& fstr);

const unsigned int MAT_FIELD_TYPE_LEN = 18; // must be the number of fields in MAT_FIELD_TYPE

void testCoherence(const mxArray* params,std::vector<bool> & presentFields);

template<typename SCALAR, typename FPP2>
const Params<SCALAR,Cpu,FPP2>* mxArray2FaustParams(const mxArray* matlab_params);

#include "mx2Faust.hpp"


#endif









