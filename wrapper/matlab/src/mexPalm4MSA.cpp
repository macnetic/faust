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


#include "mex.h"
//#include "mexutils.h"
#include "faust_MatDense.h"
#include <vector>
#include <string>
#include <algorithm>
#include "faust_constant.h"
#include "faust_Palm4MSA.h"
#include <stdexcept>
#include "mx2Faust.h"
#include "faust2Mx.h"
#include "faust_TransformHelper.h"
#include "class_handle.hpp"
using namespace Faust;

void testCoherencePalm4MSA(const mxArray* params,std::vector<bool> & presentFields);



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


	if (nrhs != 1)
	{
		mexErrMsgTxt("Bad Number of inputs arguments");
	}

    if(!mxIsStruct(prhs[0]))
    {
        mexErrMsgTxt("Input must be a structure.");
    }

    std::vector<bool> presentFields;
    testCoherencePalm4MSA(prhs[0],presentFields);
//     mexPrintf(" NUMBER FIELDS %d\n",presentFields.size());




    mxArray    *mxCurrentField,*mxCurrentCons;


    // data initialisation
    Faust::MatDense<FFPP,Cpu> data;
    if (presentFields[0])
    {
        mxCurrentField = mxGetField(prhs[0],0,"data");

        mxArray2FaustMat(  mxCurrentField,data ) ;
        /*mexPrintf("DATA");
        for (int i = 0;i<data.getNbRow();i++)
        {
            for (int j = 0;j<data.getNbCol();j++)
            {
                //bidon = std::snprintf(coeff,10,"%d",A(i,j));
                mexPrintf("%f ",data(i,j));
            }
            mexPrintf("\n")
        };*/
    }else
    {
         mexErrMsgTxt("params.data must be specified");
    }

   //nbFact initialisation
   int nbFact=0;
   if (presentFields[1])
   {

        mxCurrentField = mxGetField(prhs[0],0,"nfacts");
        nbFact =(int)  mxGetScalar(mxCurrentField);
//        mexPrintf("NB FACT : %d\n",nbFact);
   }else
   {
        mexErrMsgTxt("params.nfacts must be specified");
   }



    //constraints
    std::vector<const Faust::ConstraintGeneric<FFPP,Cpu>*> consS;
    if (presentFields[2])
    {
        mwSize nbRowCons,nbColCons;
        mxCurrentField = mxGetField(prhs[0],0,"cons");
        if(!mxIsCell(mxCurrentField))
        {
            mexErrMsgTxt("cons must be a cell-array");
        }
        nbRowCons = mxGetM(mxCurrentField);
        nbColCons = mxGetN(mxCurrentField);

        if(nbRowCons !=1)
        {

            mexErrMsgTxt("cons must have 1 rows");
        }
//		mexPrintf("cons nbColCons=%d\n", nbColCons);
        if(nbColCons != (nbFact))
        {
            //mexPrintf("\n cons has %d cols and nbFact = %d\n",nbColCons,nbFact);
            //mexErrMsgTxt("incoherence between the number of columns of cons and nfacts ");
        }
        //mexPrintf("\n cons has %d rows and %d cols \n",nbRowCons,nbColCons);
        //Faust::ConstraintGeneric * consToAdd;



            for (mwSize j=0;j<nbColCons;j++)
            {
//                mexPrintf("cons(%d)\n",j);
                mxCurrentCons=mxGetCell(mxCurrentField,j);
                getConstraint<FFPP>(consS,mxCurrentCons);
                //consS.push_back(consToAdd);
            }

    }else
    {
        mexErrMsgTxt("params.cons must be specified");
    }
//    std::cout<<"FINI_CONS"<<std::endl;
    //niter1
//    Faust::StoppingCriterion<FFPP> crit1;
//    if (presentFields[3])
//    {
//         mxCurrentField = mxGetField(prhs[0],0,"niter");
//        int nb_iter1 =(int)  mxGetScalar(mxCurrentField);
//        Faust::StoppingCriterion<FFPP> newCrit1(nb_iter1);
//        crit1 = newCrit1;
//    }
//    mexPrintf("\n crit1 nb_it = %d\n",crit1.get_crit());

	//TODO: replace by default values as constants from StoppingCriterion class
	bool is_criterion_error = false;
	int num_its = 500;
	FFPP error_treshold = 0.3;
	int max_num_its = 10000;
	if(presentFields[3])
	{
		mxCurrentField = mxGetField(prhs[0], 0, "niter");
		num_its = (int) mxGetScalar(mxCurrentField);
	}
	if(presentFields[8]){
		mxCurrentField = mxGetField(prhs[0], 0, "sc_is_criterion_error");
		is_criterion_error =  (bool) mxGetScalar(mxCurrentField);
	}
	if(presentFields[9])
	{
		mxCurrentField = mxGetField(prhs[0], 0, "sc_error_treshold");
		error_treshold = (FFPP) mxGetScalar(mxCurrentField);
	}
	if(presentFields[10])
	{
		mxCurrentField = mxGetField(prhs[0], 0, "sc_max_num_its");
		max_num_its = (int) mxGetScalar(mxCurrentField);
	}
	Faust::StoppingCriterion<FFPP> crit1(num_its, is_criterion_error, error_treshold, max_num_its);
//	crit1.Display();
    //init_facts
    std::vector<Faust::MatDense<FFPP,Cpu> > init_facts;
    if (presentFields[4])
    {
         mxCurrentField = mxGetField(prhs[0],0,"init_facts");
//		 std::cout<<"PASSERbeforeInitFact"<<std::endl;
		 setVectorFaustMat(init_facts,mxCurrentField);
//		 std::cout<<"PASSERafterInitFact"<<std::endl;

    }else
	{
		mexErrMsgTxt("init_facts must be must be specified");
	}

//	std::cout<<"PASSER1"<<std::endl;
    //verbosity
    bool isVerbose = false;
    if (presentFields[5])
    {
       mxCurrentField = mxGetField(prhs[0],0,"verbose");
        isVerbose =(bool)  mxGetScalar(mxCurrentField);
    }
    //update_way
    bool updateway = false;
    if (presentFields[7])
    {
        mxCurrentField = mxGetField(prhs[0],0,"update_way");
        updateway =(bool)  mxGetScalar(mxCurrentField);
    }


   //init_lambda
   FFPP init_lambda = 1.0;
   if (presentFields[6])
   {
       mxCurrentField = mxGetField(prhs[0],0,"init_lambda");
       init_lambda = (FFPP) mxGetScalar(mxCurrentField);
   }

   //compute_lambda
   // bool compute_lambda = true;
   // if (presentFields[8])
   // {
       // mxCurrentField = mxGetField(prhs[0],0,"compute_lambda");
       // compute_lambda = (bool) mxGetScalar(mxCurrentField);
   // }


      ///////////// HIERARCHICAL LAUNCH ///////////////
     // creation des parametres
	try{






     //creation de hierarchical fact
    Faust::ParamsPalm<FFPP,Cpu> params(data,nbFact,consS,init_facts,crit1,isVerbose,updateway,init_lambda);
//	params.Display();
	Faust::BlasHandle<Cpu> blas_handle;
	Faust::Palm4MSA<FFPP,Cpu> palm(params,blas_handle,false);
	palm.compute_facts();


     //extraction des resultats
     FFPP lambda = palm.get_lambda();

     plhs[0]=mxCreateDoubleScalar((double) lambda);

     std::vector<Faust::MatDense<FFPP,Cpu> > facts;
	 std::vector<Faust::MatGeneric<FFPP, Cpu>*> sp_facts;
     facts=palm.get_facts();

	 for(typename std::vector<Faust::MatDense<FFPP, Cpu>>::iterator it = facts.begin(); it != facts.end(); it++)
	 {
		 Faust::MatSparse<FFPP, Cpu> * M = new Faust::MatSparse<FFPP, Cpu>(*it);
		 sp_facts.push_back(M);
	 }


	 TransformHelper<FFPP, Cpu> *F = new TransformHelper<FFPP, Cpu>(sp_facts, lambda, false);
	 plhs[1] = convertPtr2Mat<Faust::TransformHelper<FFPP, Cpu>>(F);

	 for(typename std::vector<Faust::MatGeneric<FFPP,Cpu>*>::iterator it = sp_facts.begin(); it != sp_facts.end(); it++)
	 {
		 delete *it;
	 }
	}
	catch (const std::exception& e)
	{
		 mexErrMsgTxt(e.what());
	}

}









void testCoherencePalm4MSA(const mxArray* params,std::vector<bool> & presentFields)
{
  int nbr_field=mxGetNumberOfFields(params);
  presentFields.resize(8);
  presentFields.assign(8,false);
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

  }

}













