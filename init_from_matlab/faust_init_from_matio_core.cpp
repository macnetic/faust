#include "faust_init_from_matio_core.h"
#include "faust_init_from_matio_mat.h"
#include "faust_mat.h"
#include "faust_spmat.h"
#include "faust_core.h"
#include "faust_constant.h"
#include <iostream>
#include <vector>
#include "faust_params.h"
#include "faust_params_palm.h"
#include "stopping_criterion.h"
#include "faust_constraint_int.h"
#include "faust_constraint_generic.h"
#include "faust_constraint_real.h"
#include "faust_constraint_mat.h"

using namespace std;


void init_faust_core_from_matiofile(faust_core& core, const char* fileName, const char* variableName)
{
	
	matvar_t* cell_var = faust_matio_read_variable(fileName, variableName);
   
	matvar_t* current_spmat_var;
	
	char* current_fieldName;
	
	int niter,nfacts,verbose,update_way,dim1,dim2,cons_parameter,cons_dim1,cons_dim2;
	faust_spmat data_spmat;
	
	core.clear();
	
	for (int j=0 ; j<cell_var->dims[1] ; j++)
	{	
		current_spmat_var = Mat_VarGetCell(cell_var, j);
		init_spmat_from_matvar(data_spmat, current_spmat_var);
		core.push_back(data_spmat);	
	}
	
}

