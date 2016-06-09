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

template<typename FPP>
void getFaustVec(const mxArray * vec_array,Faust::Vect<FPP,Cpu> & vec);
template<typename FPP>
void getFaustMat(const mxArray* Mat_array,Faust::MatDense<FPP,Cpu> & Mat);
template<typename FPP>
void getFaustspMat(const mxArray* spMat_array,Faust::MatSparse<FPP,Cpu> & S);
template<typename FPP>
mxArray*  FaustMat2mxArray(const Faust::MatDense<FPP,Cpu>& M);
template<typename FPP>
void setCellFacts(mxArray ** cellFacts,std::vector<Faust::MatDense<FPP,Cpu> > & facts);
template<typename FPP>
void getConstraint(std::vector<const Faust::ConstraintGeneric<FPP,Cpu>*> & consS,mxArray* mxCons);
template<typename FPP>
void setVectorFaustMat(std::vector<Faust::MatDense<FPP,Cpu> > &vecMat, mxArray *Cells);
template<typename FPP>
void addSpmat(const mxArray * mxMat, std::vector<Faust::MatSparse<FPP,Cpu> > &vec_spmat);
void testCoherence(const mxArray* params,std::vector<bool> & presentFields);
template<typename FPP>
void DisplayParams(const Faust::Params<FPP,Cpu> & params);

#include "tools_mex.hpp"

#endif
