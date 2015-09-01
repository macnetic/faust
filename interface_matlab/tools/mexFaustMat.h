#include <mex.h>
#include <faust_mat.h>
#include <faust_spmat.h>
#include <vector>
#include <faust_constraint_generic.h>


void getFaustMat(const mxArray* Mat_array,faust_mat & Mat);
void getFaustspMat(const mxArray* spMat_array,faust_spmat & S);
mxArray*  FaustMat2mxArray(faust_mat M);
void setCellFacts(mxArray ** cellFacts,std::vector<faust_mat> facts);
void getConstraint(std::vector<const faust_constraint_generic*> & consS,mxArray* mxCons);
void setVectorFaustMat(std::vector<faust_mat> &vecMat, mxArray *Cells);
void loadDenseFaust( const mxArray * Cells,std::vector<faust_spmat> &vec_spmat);
void loadSpFaust(const mxArray * Cells,std::vector<faust_spmat> &vec_spmat);
void addSpmat(const mxArray * mxMat,std::vector<faust_spmat> &vec_spmat);