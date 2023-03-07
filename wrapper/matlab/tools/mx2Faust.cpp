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






//why this macro here ? it's not a header
#ifndef __FAUST_MX2FAUST_CPP__
#define __FAUST_MX2FAUST_CPP__




#include <string>
#include <stdexcept>
#include <complex>
#include "mx2Faust.h"
#include "faust_ConstraintGeneric.h"
#include "faust_ConstraintFPP.h"
#include "faust_ConstraintMat.h"
#include "faust_ConstraintInt.h"
#include "faust_Params.h"

using namespace Faust;
using namespace std;


const string mat_field_type2str(MAT_FIELD_TYPE f)
{
	// better map than an array because the order of fields doesn't matter
	switch(f)
	{
		case NROW:
			return "nrow";
		case NCOL:
			return "ncol";
		case NFACTS:
			return "nfacts";
		case CONS:
			return "cons";
		case NITER1:
			return "niter1";
		case NITER2:
			return "niter2";
		case VERBOSE:
			return "verbose";
		case FACT_SIDE:
			return "fact_side";
		case UPDATE_WAY:
			return "update_way";
		case INIT_LAMBDA:
			return "init_lambda";
		case SC_IS_CRITERION_ERROR:
			return "sc_is_criterion_error";
		case SC_ERROR_TRESHOLD:
			return "sc_error_treshold";
		case SC_MAX_NUM_ITS:
			return "sc_max_num_its";
		case SC_IS_CRITERION_ERROR2:
			return "sc_is_criterion_error2";
		case SC_ERROR_TRESHOLD2:
			return "sc_error_treshold2";
		case SC_MAX_NUM_ITS2:
			return "sc_max_num_its2";
		case INIT_FACTS:
			return "init_facts";
		case INIT_D:
			return "init_D";
		case PACKING_RL:
			return "packing_RL";
		case NO_NORMALIZATION:
			return "no_normalization";
		case NO_LAMBDA:
			return "no_lambda";
		case FACTOR_FORMAT:
			return "factor_format";
		case NORM2_MAX_ITER:
			return "norm2_max_iter";
		case NORM2_THRESHOLD:
			return "norm2_threshold";
		case CONSTANT_STEP_SIZE:
			return "constant_step_size";
		case STEP_SIZE:
			return "step_size";
		case GRAD_CALC_OPT_MODE:
			return "grad_calc_opt_mode";
		case SC_ERREPS:
			return "sc_erreps";
		case SC_ERREPS2:
			return "sc_erreps2";
	}
}


const MAT_FIELD_TYPE mat_field_str2type(const string& fstr)
{
	for(int i=0;i<MAT_FIELD_TYPE_LEN;i++)
	{
		MAT_FIELD_TYPE ft = static_cast<MAT_FIELD_TYPE>(i);
		if(fstr == mat_field_type2str(ft))
			return ft;
	}
	throw invalid_argument("Invalid matlab struct field name: "+fstr);
}



void testCoherence(const mxArray* params,std::vector<bool> & presentFields)
{
  int nbr_field=mxGetNumberOfFields(params);
  presentFields.resize(MAT_FIELD_TYPE_LEN);
  presentFields.assign(MAT_FIELD_TYPE_LEN,false);

  if(nbr_field < 3)
      mexErrMsgTxt("The number of fields in params must be at least 3 ");

  for(int i=0;i<nbr_field;i++)
  {
	  auto field_name = string(mxGetFieldNameByNumber(params,i));
	  // dirty bypass to tolerate fields starting with mhtp (MHTPParams) parsed by mxArray2FaustMHTPParams
	  if (field_name.rfind("mhtp_", 0) != 0)
		  presentFields[mat_field_str2type(field_name)] = true;
  }

}

void testCoherencePALM4MSA(const mxArray* params,std::vector<bool> & presentFields)
{
	////TODO: this function should be modified to be more reliable/simple as the function testCoherence() in mx2Faust.cpp has been modified
	int nbr_field=mxGetNumberOfFields(params);
	presentFields.resize(20);
	presentFields.assign(20,false);
	if(nbr_field < 3)
	{
		mexErrMsgTxt("The number of field of params must be at least 3 ");
	}


	for (int i=0;i<nbr_field;i++)
	{
		const char * fieldName;
		fieldName = mxGetFieldNameByNumber(params,i);
		//        mexPrintf("fieldname %d : %s\n",i,fieldName);

		if (strcmp(fieldName,"data") == 0)
		{
			presentFields[0] = true;
		}

		else if (strcmp(fieldName,"nfacts") == 0)
		{
			presentFields[1] = true;
		}
		else if (strcmp(fieldName,"cons") == 0)
		{
			presentFields[2] = true;
		}
		else if (strcmp(fieldName,"niter") == 0)
		{
			presentFields[3] = true;
		}
		//sc stands for StoppingCriterion
		else if(!strcmp(fieldName, "sc_is_criterion_error"))
		{
			presentFields[8] = true;
		}
		else if(!strcmp(fieldName, "sc_error_treshold"))
		{
			presentFields[9] = true;
		}
		else if(!strcmp(fieldName, "sc_max_num_its"))
		{
			presentFields[10] = true;
		}
		else if (strcmp(fieldName,"init_facts") == 0)
		{
			presentFields[4] = true;
		}
		else if (strcmp(fieldName,"verbose") == 0)
		{
			presentFields[5] = true;
		}
		else if (strcmp(fieldName,"init_lambda") == 0)
		{
			presentFields[6] = true;
		}
		else if (strcmp(fieldName,"update_way") == 0)
		{
			presentFields[7] = true;
		}
		else if(strcmp(fieldName, "grad_calc_opt_mode") == 0)
		{
			presentFields[11] = true;
		}
		else if(strcmp(fieldName, "constant_step_size") == 0)
		{
			presentFields[12] = true;
		}
		else if(strcmp(fieldName, "step_size") == 0)
		{
			presentFields[13] = true;
		}
		else if(strcmp(fieldName, "norm2_max_iter") == 0)
			presentFields[14] = true;
		else if(strcmp(fieldName, "norm2_threshold") == 0)
			presentFields[15] = true;
		else if(strcmp(fieldName, "factor_format") == 0)
			presentFields[16] = true;
		else if(strcmp(fieldName, "packing_RL") == 0)
			presentFields[17] = true;
		else if(strcmp(fieldName, "no_normalization") == 0)
			presentFields[18] = true;
		else if(strcmp(fieldName, "no_lambda") == 0)
			presentFields[19] = true;
		else if(strcmp(fieldName, "sc_erreps") == 0)
			presentFields[20] = true;
	}

}


// void DisplayParams(const Faust::Params & params)
// {
    // mexPrintf("/////// PARAMS //////\n");
    // Faust::MatDense data = params.data;

    // int nbRow = data.getNbRow();
    // int nbCol = data.getNbCol();
    // mexPrintf("DATA DIM : %d %d\n",nbRow,nbCol);
    // if (nbRow*nbCol < 1000)
    // {
        // for (int i = 0;i<data.getNbRow();i++)
        // {
            // for (int j = 0;j<data.getNbCol();j++)
            // {
                // mexPrintf("%f ",data(i,j));
            // }
            // mexPrintf("\n");
        // }
    // }
    // mexPrintf("\n\n");

    // mexPrintf("NB_FACTS : %d\n",params.m_nbFact);
    // mexPrintf("CONSTRAINTS : nbr %d\n",params.cons[0].size());
    // for (unsigned int L=0;L<params.cons[0].size();L++)
	// {
		// for (unsigned int jl=0;jl<params.cons.size();jl++)
		// {
			// mexPrintf(" %s ", params.cons[jl][L]->get_constraint_name());
			// mexPrintf(" DIMS (%d,%d) : ",(*params.cons[jl][L]).get_rows(),(*params.cons[jl][L]).get_cols());


			// if (params.cons[jl][L]->is_constraint_parameter_int())
			// {
				// Faust::ConstraintInt* const_int = (Faust::ConstraintInt*)(params.cons[jl][L]);
				// mexPrintf(" parameter : %d",(*const_int).get_parameter());
			// }

			// if (params.cons[jl][L]->is_constraint_parameter_real())
			// {
				// Faust::ConstraintFPP* const_real = (Faust::ConstraintFPP*)(params.cons[jl][L]);
				// mexPrintf(" parameter : %f",(*const_real).get_parameter());
			// }
			// mexPrintf("\n");


		// }

	// }
// }

#endif


