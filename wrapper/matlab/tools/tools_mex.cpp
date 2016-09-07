/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.gforge.inria.fr>                         */
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













#include "tools_mex.h"
#include "faust_ConstraintGeneric.h"
#include "faust_ConstraintFPP.h"
#include "faust_ConstraintMat.h"
#include "faust_ConstraintInt.h"
#include "faust_Params.h"










void testCoherence(const mxArray* params,std::vector<bool> & presentFields)
{
  int nbr_field=mxGetNumberOfFields(params);
  presentFields.resize(9);
  presentFields.assign(9,false);
  if(nbr_field < 3)
  {
      mexErrMsgTxt("The number of field of params must be at least 3 ");
  }


  for (int i=0;i<nbr_field;i++)
  {
        const char * fieldName;
        fieldName = mxGetFieldNameByNumber(params,i);
        //mexPrintf("fieldname %d : %s\n",i,fieldName);

        if (strcmp(fieldName,"data") == 0)
        {
            presentFields[0] = true;
        }

        if (strcmp(fieldName,"nfacts") == 0)
        {
            presentFields[1] = true;
        }
        if (strcmp(fieldName,"cons") == 0)
        {
            presentFields[2] = true;
        }
        if (strcmp(fieldName,"niter1") == 0)
        {
            presentFields[3] = true;
        }
        if (strcmp(fieldName,"niter2") == 0)
        {
            presentFields[4] = true;
        }
        if (strcmp(fieldName,"verbose") == 0)
        {
            presentFields[5] = true;
        }
        if (strcmp(fieldName,"fact_side") == 0)
        {
            presentFields[6] = true;
        }
        if (strcmp(fieldName,"update_way") == 0)
        {
            presentFields[7] = true;
        }
        if (strcmp(fieldName,"init_lambda") == 0)
        {
            presentFields[8] = true;
        }
        // if (strcmp(fieldName,"compute_lambda") == 0)
        // {
            // presentFields[9] = true;
        // }
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


