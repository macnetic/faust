// #include "faust_init_from_matio_core.h"
// #include "faust_init_from_matio_mat.h"
// #include "faust_mat.h"
// #include "faust_spmat.h"
// #include "faust_core.h"
// #include "faust_constant.h"
// #include <iostream>
// #include <vector>
// #include "faust_params.h"
// #include "faust_params_palm.h"
// #include "stopping_criterion.h"
// #include "faust_constraint_int.h"
// #include "faust_constraint_generic.h"
// #include "faust_constraint_real.h"
// #include "faust_constraint_mat.h"

// using namespace std;


// void init_faust_core_from_matiofile(faust_core& core, const char* fileName, const char* variableName)
// {
	
	// matvar_t* cell_var = faust_matio_read_variable(fileName, variableName);

	// init_faust_core_from_matvar(core, cell_var);

	// Mat_VarFree(cell_var);
// }


// void init_faust_core_from_matvar(faust_core& core, matvar_t* cell_var )
// {
	// if(cell_var->class_type != MAT_C_CELL
		// || cell_var->rank != 2
		// || cell_var->data_size != sizeof(double))
	// {
		// cout << "Error in init_faust_core_from_matiofile : filename seems not to be a cell" << endl;
		// exit(EXIT_FAILURE);
	
	// }
   
	// matvar_t* current_spmat_var;
	
	
	// faust_spmat data_spmat;
	
	// core.clear();
	// if(cell_var->dims[0] != 1 )
	// {
		// cout << "Error in init_faust_core_from_matiofile :  filename seems not to be a row vector" << endl;
		// exit(EXIT_FAILURE);
	// }	
	// for (int j=0 ; j<cell_var->dims[1] ; j++)
	// {	
		// current_spmat_var = Mat_VarGetCell(cell_var, j);
		// init_spmat_from_matvar(data_spmat, current_spmat_var);
		// core.push_back(data_spmat);
	// }
		
// }


// void init_faust_data_from_matiofile(vector<faust_mat>& full_mat, vector<faust_core>& core, const char* fileName, const char* variableName)
// {

	
	// matvar_t* cell_array_var = faust_matio_read_variable(fileName, variableName);

	// if(cell_array_var->class_type != MAT_C_CELL
		// || cell_array_var->rank != 2
		// || cell_array_var->data_size != sizeof(double))
	// {
		// cout << "Error in init_faust_data_from_matiofile : " << fileName << "seems not to be a cell" << endl;
		// exit(EXIT_FAILURE);
	
	// }

	// matvar_t* current_cell_var;

	// matvar_t* current_mat_var;
	
	
	// faust_spmat data_spmat;
	
	// core.clear();
	// full_mat.clear();

	// if(cell_array_var->dims[0] != 2)
	// {
		// cerr << "Error in init_faust_data_from_matiofile : wrong dimensions of cell aray" << endl;
		// exit(EXIT_FAILURE);
	// }
		
	// //creation des matrices pleines full_mat a partir de la premiere ligne du tableau de cellules (1ere ligne de cell_array_var)
	// for (int j=0 ; j<cell_array_var->dims[1] ; j++)
	// {	
		// current_mat_var = Mat_VarGetCell(cell_array_var, j*cell_array_var->dims[0]);
		// faust_mat mat_tmp;
		// init_mat_from_matvar(mat_tmp, current_mat_var);
		// full_mat.push_back(mat_tmp);	
	// }


	// //creation des objets faust_core core a partir de la deuxieme ligne du tableau de cellules (2eme ligne de cell_array_var)
	// for (int j=0 ; j<cell_array_var->dims[1] ; j++)
	// {	
		// current_cell_var = Mat_VarGetCell(cell_array_var, j*cell_array_var->dims[0]+1);
		// faust_core core_tmp;
		// init_faust_core_from_matvar(core_tmp, current_cell_var );
		// core.push_back(core_tmp);	
	// }
	
// }

