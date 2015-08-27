#include <mex.h>
#include <faust_mat.h>


faust_mat getFaustMat(mxArray* Mat_array);
mxArray*  FaustMat2mxArray(faust_mat M);