#ifndef __FAUST_TOOLS_MEX_H__
#define __FAUST_TOOLS_MEX_H__

#include "mex.h"
#include "faust_mat.h"
#include "faust_spmat.h"
#include <vector>

class faust_constraint_generic;
class faust_params;


void getFaustMat(const mxArray* Mat_array,faust_mat & Mat);
void getFaustspMat(const mxArray* spMat_array,faust_spmat & S);
mxArray*  FaustMat2mxArray(const faust_mat& M);
void setCellFacts(mxArray ** cellFacts,std::vector<faust_mat> facts);
void getConstraint(std::vector<const faust_constraint_generic*> & consS,mxArray* mxCons);
void setVectorFaustMat(std::vector<faust_mat> &vecMat, mxArray *Cells);
void addSpmat(const mxArray * mxMat,std::vector<faust_spmat> &vec_spmat);
void testCoherence(const mxArray* params,std::vector<bool> & presentFields);
void DisplayParams(const faust_params & params);
#endif
