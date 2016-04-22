#ifndef __FAUST_TOOLS_MEX_H__
#define __FAUST_TOOLS_MEX_H__

#include "mex.h"
#include "faust_mat.h"
#include "faust_spmat.h"
#include <vector>
#include "faust_vec.h"

template<typename T> class faust_constraint_generic;
template<typename T> class faust_vec;
template<typename T> class faust_params;
template<typename T> class faust_mat;
template<typename T> class faust_spmat;

template<typename T>
void getFaustVec(const mxArray * vec_array,faust_vec<T> & vec);
template<typename T>
void getFaustMat(const mxArray* Mat_array,faust_mat<T> & Mat);
template<typename T>
void getFaustspMat(const mxArray* spMat_array,faust_spmat<T> & S);
template<typename T>
mxArray*  FaustMat2mxArray(const faust_mat<T>& M);
template<typename T>
void setCellFacts(mxArray ** cellFacts,std::vector<faust_mat<T> > & facts);
template<typename T>
void getConstraint(std::vector<const faust_constraint_generic<T>*> & consS,mxArray* mxCons);
template<typename T>
void setVectorFaustMat(std::vector<faust_mat<T> > &vecMat, mxArray *Cells);
template<typename T> 
void addSpmat(const mxArray * mxMat, std::vector<faust_spmat<T> > &vec_spmat);
void testCoherence(const mxArray* params,std::vector<bool> & presentFields);
template<typename T>
void DisplayParams(const faust_params<T> & params);

#include "tools_mex.hpp"

#endif
