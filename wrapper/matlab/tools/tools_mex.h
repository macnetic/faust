#ifndef __FAUST_TOOLS_MEX_H__
#define __FAUST_TOOLS_MEX_H__

#include "mex.h"
#include "faust_MatDense.h"
#include "faust_MatSparse.h"
#include <vector>
#include "faust_Vect.h"

template<typename T> class ConstraintGeneric;
template<typename T> class Vect;
template<typename T> class Params;
template<typename T> class MatDense;
template<typename T> class MatSparse;

template<typename T>
void getFaustVec(const mxArray * vec_array,Faust::Vect<T> & vec);
template<typename T>
void getFaustMat(const mxArray* Mat_array,Faust::MatDense<T> & Mat);
template<typename T>
void getFaustspMat(const mxArray* spMat_array,Faust::MatSparse<T> & S);
template<typename T>
mxArray*  FaustMat2mxArray(const Faust::MatDense<T>& M);
template<typename T>
void setCellFacts(mxArray ** cellFacts,std::vector<Faust::MatDense<T> > & facts);
template<typename T>
void getConstraint(std::vector<const Faust::ConstraintGeneric<T>*> & consS,mxArray* mxCons);
template<typename T>
void setVectorFaustMat(std::vector<Faust::MatDense<T> > &vecMat, mxArray *Cells);
template<typename T>
void addSpmat(const mxArray * mxMat, std::vector<Faust::MatSparse<T> > &vec_spmat);
void testCoherence(const mxArray* params,std::vector<bool> & presentFields);
template<typename T>
void DisplayParams(const Faust::Params<T> & params);

#include "tools_mex.hpp"

#endif
