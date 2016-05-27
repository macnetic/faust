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
        mexPrintf("fieldname %d : %s\n",i,fieldName);

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

    // mexPrintf("NB_FACTS : %d\n",params.nb_fact);
    // mexPrintf("CONSTRAINTS : nbr %d\n",params.cons[0].size());
    // for (unsigned int L=0;L<params.cons[0].size();L++)
	// {
		// for (unsigned int jl=0;jl<params.cons.size();jl++)
		// {
			// mexPrintf(" %s ", params.cons[jl][L]->get_constraint_name());
			// mexPrintf(" DIMS (%d,%d) : ",(*params.cons[jl][L]).getRows(),(*params.cons[jl][L]).getCols());


			// if (params.cons[jl][L]->isConstraintParameterInt())
			// {
				// Faust::ConstraintInt* const_int = (Faust::ConstraintInt*)(params.cons[jl][L]);
				// mexPrintf(" parameter : %d",(*const_int).getParameter());
			// }

			// if (params.cons[jl][L]->isConstraintParameterReal())
			// {
				// Faust::ConstraintFPP* const_real = (Faust::ConstraintFPP*)(params.cons[jl][L]);
				// mexPrintf(" parameter : %f",(*const_real).getParameter());
			// }
			// mexPrintf("\n");


		// }

	// }
// }


