/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
/*                      Luc Le Magoarou, Remi Gribonval                     */
/*                      INRIA Rennes, FRANCE                                */
/*                      http://www.inria.fr/                                */
/*                                                                          */
/*  The FAuST Toolbox is distributed under the terms of the GNU Affero      */
/*  General Public License.                                                 */
/*  This program is free software: you can redistribute it and/or modify    */
/*  it under the terms of the GNU Affero General Public License as          */
/*  published by the Free Software Foundation.                              */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful, but     */
/*  WITHOUT ANY WARRANTY; without even the implied warranty of              */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    */
/*  See the GNU Affero General Public License for more details.             */
/*                                                                          */
/*  You should have received a copy of the GNU Affero General Public        */
/*  License along with this program.                                        */
/*  If not, see <http://www.gnu.org/licenses/>.                             */
/*                                                                          */
/*                             Contacts:                                    */
/*      Nicolas Bellot  : nicolas.bellot@inria.fr                           */
/*      Adrien Leman    : adrien.leman@inria.fr                             */
/*      Thomas Gautrais : thomas.gautrais@inria.fr                          */
/*      Luc Le Magoarou : luc.le-magoarou@inria.fr                          */
/*      Remi Gribonval  : remi.gribonval@inria.fr                           */
/*                                                                          */
/*                              References:                                 */
/*  [1] Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse       */
/*  approximations of matrices and applications", Journal of Selected       */
/*  Topics in Signal Processing, 2016.                                      */
/*  <https://hal.archives-ouvertes.fr/hal-01167948v1>                       */
/****************************************************************************/


#ifndef __FAUST_INIT_FROM_MATIO_CORE_HPP__
#define __FAUST_INIT_FROM_MATIO_CORE_HPP__

//#include "faust_init_from_matio_core.h"
#include "faust_ConstraintGeneric.h"
#include "faust_init_from_matio_mat.h"
#include "faust_MatDense.h"
#include "faust_MatSparse.h"
#include "faust_Transform.h"
#include "faust_constant.h"
#include <iostream>
#include <vector>
#include "faust_Params.h"
#include "faust_ParamsPalm.h"
#include "faust_StoppingCriterion.h"
#include "faust_ConstraintInt.h"
#include "faust_ConstraintFPP.h"
#include "faust_ConstraintMat.h"

#ifdef COMPILE_GPU
	#include   "faust_gpuCore2cpuCore.h"
#endif



#ifdef COMPILE_GPU
template<typename FPP>
void write_faust_core_into_matfile(const Faust::Transform<FPP,Gpu> core, const char* fileName, const char* variableName)
{
	Faust::Transform<FPP,Cpu> fcore;
	faust_gpu2cpu(fcore,core);
	write_faust_core_into_matfile(fcore,fileName,variableName);

}
#endif

using namespace std;

template<typename FPP,Device DEVICE>
void init_faust_core_from_matiofile(Faust::Transform<FPP,DEVICE>& core, const char* fileName, const char* variableName)
{
	matvar_t* cell_var = faust_matio_read_variable(fileName, variableName);

	init_faust_core_from_matvar(core, cell_var);

	Mat_VarFree(cell_var);
}

template<typename FPP,Device DEVICE>
void init_faust_core_from_matvar(Faust::Transform<FPP,DEVICE>& core, matvar_t* cell_var )
{
	if(cell_var->class_type != MAT_C_CELL
		|| cell_var->rank != 2
		|| cell_var->data_size != sizeof(double))
	{
		cout << "Error in init_faust_core_from_matiofile : filename seems not to be a cell" << endl;
		exit(EXIT_FAILURE);

	}

	matvar_t* current_spmat_var;


	Faust::MatSparse<FPP,DEVICE> data_spmat;
	core.clear();
	if(cell_var->dims[0] != 1 )
	{
		cout << "Error in init_faust_core_from_matiofile :  filename seems not to be a row vector" << endl;
		exit(EXIT_FAILURE);
	}
	for (int j=0 ; j<cell_var->dims[1] ; j++)
	{
		current_spmat_var = Mat_VarGetCell(cell_var, j);
		init_spmat_from_matvar(data_spmat, current_spmat_var);
		Faust::MatGeneric<FPP,DEVICE> * ptr_data_spmat = &data_spmat;
		core.push_back(ptr_data_spmat);
	}

}

template<typename FPP,Device DEVICE>
void init_faust_data_from_matiofile(vector<Faust::MatDense<FPP,DEVICE> >& full_mat, vector<Faust::Transform<FPP,DEVICE> >& core, const char* fileName, const char* variableName)
{


	matvar_t* cell_array_var = faust_matio_read_variable(fileName, variableName);

	if(cell_array_var->class_type != MAT_C_CELL
		|| cell_array_var->rank != 2
		|| cell_array_var->data_size != sizeof(double))
	{
		cout << "Error in init_faust_data_from_matiofile : " << fileName << "seems not to be a cell" << endl;
		exit(EXIT_FAILURE);

	}

	matvar_t* current_cell_var;

	matvar_t* current_mat_var;


	Faust::MatSparse<FPP,DEVICE> data_spmat;

	core.clear();
	full_mat.clear();

	if(cell_array_var->dims[0] != 2)
	{
		cerr << "Error in init_faust_data_from_matiofile : wrong dimensions of cell aray" << endl;
		exit(EXIT_FAILURE);
	}

	// creation des matrices pleines full_mat a partir de la premiere ligne du tableau de cellules (1ere ligne de cell_array_var)
	for (int j=0 ; j<cell_array_var->dims[1] ; j++)
	{
		current_mat_var = Mat_VarGetCell(cell_array_var, j*cell_array_var->dims[0]);
		Faust::MatDense<FPP,DEVICE> mat_tmp;
		init_mat_from_matvar(mat_tmp, current_mat_var);
		full_mat.push_back(mat_tmp);
	}


	// creation des objets Faust::Transform core a partir de la deuxieme ligne du tableau de cellules (2eme ligne de cell_array_var)
	for (int j=0 ; j<cell_array_var->dims[1] ; j++)
	{
		current_cell_var = Mat_VarGetCell(cell_array_var, j*cell_array_var->dims[0]+1);
		Faust::Transform<FPP,DEVICE> core_tmp;
		init_faust_core_from_matvar(core_tmp, current_cell_var );
		core.push_back(core_tmp);
	}

}



template<typename FPP>
void write_faust_core_into_matfile(const Faust::Transform<FPP,Cpu> core, const char* fileName, const char* variableName)
{
	mat_t* matfp = Mat_Open(fileName,MAT_ACC_RDWR);
   matvar_t *matvar;




   if(matfp == NULL)
   {
		matfp = Mat_CreateVer(fileName,NULL,MAT_FT_DEFAULT);
		if ( NULL == matfp ) {
			cerr << "error in write_faust_mat<FPP,DEVICE>_into_matfile : unable to create "<< fileName << endl;
			 exit(EXIT_FAILURE);
		}
	}


	while ( (matvar = Mat_VarReadNextInfo(matfp)) != NULL ) {
		if (strcmp(matvar->name,variableName) == 0)
		{
			Mat_VarDelete(matfp,matvar->name);
		}
        matvar = NULL;
    }

    write_faust_core_into_matvar(core,&matvar,variableName);
		Mat_VarWrite(matfp,matvar,MAT_COMPRESSION_NONE);

        Mat_VarFree(matvar);


	Mat_Close(matfp);

}


template<typename FPP>
void write_faust_core_into_matvar(const Faust::Transform<FPP,Cpu> core, matvar_t** matvar, const char* variableName)
{

	
	// TODO : not  compatible with faust mat generic 
	cerr << "error in write_faust_core_into_matvar(...) : TODO NOT COMPATIBLE with Faust_Transform since they used MatGeneric as factor "<< endl;
	exit(EXIT_FAILURE);
	std::vector<Faust::MatSparse<FPP,Cpu> > sparse_facts;	
	//core.get_facts(sparse_facts);
	write_faust_spmat_list_into_matvar(sparse_facts,matvar,variableName);

}




#endif
