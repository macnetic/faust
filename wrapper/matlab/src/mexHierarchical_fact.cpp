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


#include "mex.h"
#include "faust_HierarchicalFact.h"
#include <vector>
#include <string>
#include <algorithm>
#include "tools_mex.h"
#include <stdexcept>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	#ifdef FAUST_VERBOSE
	if (typeid(FFPP) == typeid(float))
	{
		std::cout<<"FFPP == float"<<std::endl;

	}

	if (typeid(FFPP) == typeid(double))
	{
		std::cout<<"FFPP == double"<<std::endl;
	}
	system("sleep 7");
	#endif

	if (nrhs != 2)
	{
		mexErrMsgTxt("Bad Number of inputs arguments");
	}
    const mxArray* matlab_matrix = prhs[0];	
    const mxArray* matlab_params = prhs[1];  	
    if(!mxIsStruct(matlab_params))
    {
        mexErrMsgTxt("Input must be a structure.");
    }

    std::vector<bool> presentFields;
    // test if all the needed parameter are present
    testCoherence(matlab_params,presentFields);


     // mexPrintf(" NUMBER FIELDS %d\n",presentFields.size());
//     for (int i=0;i<presentFields.size();i++)
//     {
//
//         mexPrintf(" FIELDS n %d present : ",i);
//         if (presentFields[i])
//         {
//             mexPrintf(" true\n");
//         }else
//         {
//             mexPrintf(" false\n");
//         }
//     }
//

    mxArray    *mxCurrentField,*mxCurrentCons;
	
    // initialization of the matrix that will be factorized
    Faust::MatDense<FFPP,Cpu> matrix;
    mxArray2FaustMat(matlab_matrix,matrix);
	 
    // nrow of the matrix that will be factorized	
    int nb_row = -1;
    if (presentFields[0])
    {
        mxCurrentField = mxGetField(matlab_params,0,"nrow");
        nb_row =(int)  mxGetScalar(mxCurrentField);
    }else
    {
		mexErrMsgTxt("params.nrow must be specified");
    }

    //nb column of the matrix that will be factorized
    int nb_col = -1;	
    if (presentFields[1])
    {
        mxCurrentField = mxGetField(matlab_params,0,"ncol");
        nb_col =(int)  mxGetScalar(mxCurrentField);
    }else
    {
		mexErrMsgTxt("params.nrow must be specified");
    }	


   //nbFact initialisation
   int nbFact = -1;
   if (presentFields[2])
   {

        mxCurrentField = mxGetField(matlab_params,0,"nfacts");
        nbFact =(int)  mxGetScalar(mxCurrentField);

   }else
   {
        mexErrMsgTxt("params.nfacts must be specified");
   }


    //constraints
    std::vector<std::vector<const Faust::ConstraintGeneric<FFPP,Cpu>*> > consSS;
    if (presentFields[3])
    {
        mwSize nbRowCons,nbColCons;
        mxCurrentField = mxGetField(matlab_params,0,"cons");

        if(!mxIsCell(mxCurrentField))
        {
            mexErrMsgTxt("cons must be a cell-array");
        }
        nbRowCons = mxGetM(mxCurrentField);
        nbColCons = mxGetN(mxCurrentField);

        /*if(nbRowCons !=2)
        {
            mexPrintf("\n cons has %d rows \n",nbRowCons);
            mexErrMsgTxt("cons must have 2 rows");
        }*/
        /*if(nbColCons != (nbFact-1))
        {
            mexPrintf("\n cons has %d cols and nbFact = %d\n",nbColCons,nbFact);
            mexErrMsgTxt("incoherence between the number of columns of cons and nfacts ");
        }*/
        //mexPrintf("\n cons has %d rows and %d cols \n",nbRowCons,nbColCons);
        //Faust::ConstraintGeneric * consToAdd;
        std::vector<const Faust::ConstraintGeneric<FFPP,Cpu>*> consS;

        for (mwSize i=0;i<nbRowCons;i++)
        {
            //mexPrintf("%d / %d\n",i,nbRowCons);
            for (mwSize j=0;j<nbColCons;j++)
            {
                //mexPrintf("cons(%d , %d)\n",i,j);
                mxCurrentCons=mxGetCell(mxCurrentField,i+(j*nbRowCons));
                getConstraint<FFPP>(consS,mxCurrentCons);
                //consS.push_back(consToAdd);
            }
            consSS.push_back(consS);
            consS.resize(0);
        }

    }else
    {
        mexErrMsgTxt("params.cons must be specified");
    }

    //niter1
    Faust::StoppingCriterion<FFPP> crit1;
    if (presentFields[4])
    {
         mxCurrentField = mxGetField(matlab_params,0,"niter1");
        int nb_iter1 =(int)  mxGetScalar(mxCurrentField);
        Faust::StoppingCriterion<FFPP> newCrit1(nb_iter1);
        crit1 = newCrit1;
    }
    //mexPrintf("\n crit1 nb_it = %d\n",crit1.get_crit());

    //niter2
    Faust::StoppingCriterion<FFPP> crit2;
    if (presentFields[5])
    {
         mxCurrentField = mxGetField(matlab_params,0,"niter2");
        int nb_iter2 =(int)  mxGetScalar(mxCurrentField);
        Faust::StoppingCriterion<FFPP> newCrit2(nb_iter2);
        crit2 = newCrit2;
    }
    //mexPrintf("\n crit2 nb_it = %d\n",crit2.get_crit());

    //verbosity
    bool isVerbose = false;
    if (presentFields[6])
    {
       mxCurrentField = mxGetField(matlab_params,0,"verbose");
        isVerbose =(bool)  mxGetScalar(mxCurrentField);
    }

    //fact_side
    bool factside = false;
    if (presentFields[7])
    {
      mxCurrentField = mxGetField(matlab_params,0,"fact_side");
      factside =(bool)  mxGetScalar(mxCurrentField);
    }

    //update_way
    bool updateway = false;
    if (presentFields[8])
    {
        mxCurrentField = mxGetField(matlab_params,0,"update_way");
        updateway =(bool)  mxGetScalar(mxCurrentField);
    }


   //init_lambda
   FFPP init_lambda = (FFPP) 1.0;
   if (presentFields[9])
   {
       mxCurrentField = mxGetField(matlab_params,0,"init_lambda");
       init_lambda = (FFPP) mxGetScalar(mxCurrentField);
   }

   //compute_lambda
   // bool compute_lambda = true;
   // if (presentFields[9])
   // {
       // mxCurrentField = mxGetField(matlab_params,0,"compute_lambda");
       // compute_lambda = (bool) mxGetScalar(mxCurrentField);
   // }


      ///////////// HIERARCHICAL LAUNCH ///////////////
     // creation des parametres
	 try{
		std::cout<<"nb_row : "<<nb_row<<std::endl;
		std::cout<<"nb_col : "<<nb_col<<std::endl;
		Faust::Params<FFPP,Cpu> params(nb_row,nb_col,nbFact,consSS,std::vector<Faust::MatDense<FFPP,Cpu> >(),crit1,crit2,isVerbose,updateway,factside,init_lambda);

	 //DisplayParams(params);
     //creation de hierarchical fact
     Faust::BlasHandle<Cpu> blas_handle;
     Faust::SpBlasHandle<Cpu> spblas_handle;

     Faust::HierarchicalFact<FFPP,Cpu> hier_fact(matrix,params,blas_handle,spblas_handle);
     hier_fact.compute_facts();


     //extraction des resultats
     FFPP lambda = hier_fact.get_lambda();

     plhs[0]=mxCreateDoubleScalar((double) lambda);

     std::vector<Faust::MatDense<FFPP,Cpu> > facts;
     hier_fact.get_facts(facts);


     Faust::MatDense<FFPP,Cpu> current_fact = facts[0];
     mxArray * cellFacts;
     setCellFacts(&cellFacts,facts);


     plhs[1]=cellFacts;

	}
	 catch (const std::exception& e)
	{
		 mexErrMsgTxt(e.what());
	}

}

