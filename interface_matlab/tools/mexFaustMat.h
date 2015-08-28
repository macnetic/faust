#include <mex.h>
#include <faust_mat.h>
#include <faust_spmat.h>
#include <vector>
#include <faust_constraint_generic.h>


faust_mat getFaustMat(mxArray* Mat_array);
faust_spmat getFaustspMat(mxArray* spMat_array);
mxArray*  FaustMat2mxArray(faust_mat M);
void setCellFacts(mxArray ** cellFacts,std::vector<faust_mat> facts);
void getConstraint(std::vector<const faust_constraint_generic*> & consS,mxArray* mxCons);
void setVectorFaustMat(std::vector<faust_mat> &vecMat, mxArray *Cells);