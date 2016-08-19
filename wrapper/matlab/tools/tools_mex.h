#ifndef __FAUST_TOOLS_MEX_H__
#define __FAUST_TOOLS_MEX_H__

#include "mex.h"
#include <vector>
#include "faust_constant.h"

namespace Faust {
	template<typename FPP, Device DEVICE> class ConstraintGeneric;
	template<typename FPP, Device DEVICE> class Vect;
	template<typename FPP, Device DEVICE> class Params;
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
        *  \brief convert the matlab mxMat* into a Faust::MatSparse<FPP,Cpu> and push it to the end of the 			list vec_spmat
        *  \param * mxMat : pointer to the matlab matrix format
	*  \tparam vec_spmat : std::vector of Faust::MatSparse<FPP,Cpu>
        */
template<typename FPP>
void addSpmat(const mxArray * mxMat, std::vector<Faust::MatSparse<FPP,Cpu> > &vec_spmat);

void testCoherence(const mxArray* params,std::vector<bool> & presentFields);
template<typename FPP>
void DisplayParams(const Faust::Params<FPP,Cpu> & params);

#include "tools_mex.hpp"

#endif
