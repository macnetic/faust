#include "mex.h"
//#include "mexutils.h"
#include "HierarchicalFact.h"
#include "faust_MatDense.h"
#include "faust_Tools.h"
#include "faust_ConstraintInt.h"
#include "faust_ConstraintFPP.h"
#include "faust_ConstraintMat.h"
#include "faust_ConstraintGeneric.h"
#include <vector>
#include <string>
#include <algorithm>
#include "faust_Params.h"
#include "faust_constant.h"
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

	if (nrhs != 1)
	{
		mexErrMsgTxt("Bad Number of inputs arguments");
	}

    if(!mxIsStruct(prhs[0]))
    {
        mexErrMsgTxt("Input must be a structure.");
    }

    std::vector<bool> presentFields;
    testCoherence(prhs[0],presentFields);
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

    // data initialisation
    Faust::MatDense<FFPP> data;
    if (presentFields[0])
    {
        mxCurrentField = mxGetField(prhs[0],0,"data");

        getFaustMat(  mxCurrentField,data ) ;
        mexPrintf("DATA (%d,%d)",data.getNbRow(),data.getNbCol());
		if ((data.getNbRow() < 10) && (data.getNbCol()))
		{
			for (int i = 0;i<data.getNbRow();i++)
			{
				for (int j = 0;j<data.getNbCol();j++)
				{
                //bidon = std::snprintf(coeff,10,"%d",A(i,j));
					mexPrintf("%f ",data(i,j));
				}
				mexPrintf("\n");
			}
		}
	}else
	{
		mexErrMsgTxt("params.data must be specified");
	}

   //nb_fact initialisation
   int nb_fact = -1;
   if (presentFields[1])
   {

        mxCurrentField = mxGetField(prhs[0],0,"nfacts");
        mexPrintf("a\n");
        nb_fact =(int)  mxGetScalar(mxCurrentField);
        mexPrintf("NB FACT : %d",nb_fact);
   }else
   {
        mexErrMsgTxt("params.nfacts must be specified");
   }


    //constraints
    std::vector<std::vector<const Faust::ConstraintGeneric<FFPP>*> > consSS;
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

        /*if(nbRowCons !=2)
        {
            mexPrintf("\n cons has %d rows \n",nbRowCons);
            mexErrMsgTxt("cons must have 2 rows");
        }*/
        /*if(nbColCons != (nb_fact-1))
        {
            mexPrintf("\n cons has %d cols and nb_fact = %d\n",nbColCons,nb_fact);
            mexErrMsgTxt("incoherence between the number of columns of cons and nfacts ");
        }*/
        //mexPrintf("\n cons has %d rows and %d cols \n",nbRowCons,nbColCons);
        //Faust::ConstraintGeneric * consToAdd;
        std::vector<const Faust::ConstraintGeneric<FFPP>*> consS;

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
    StoppingCriterion<FFPP> crit1;
    if (presentFields[3])
    {
         mxCurrentField = mxGetField(prhs[0],0,"niter1");
        int nb_iter1 =(int)  mxGetScalar(mxCurrentField);
        StoppingCriterion<FFPP> newCrit1(nb_iter1);
        crit1 = newCrit1;
    }
    //mexPrintf("\n crit1 nb_it = %d\n",crit1.get_crit());

    //niter2
    StoppingCriterion<FFPP> crit2;
    if (presentFields[4])
    {
         mxCurrentField = mxGetField(prhs[0],0,"niter2");
        int nb_iter2 =(int)  mxGetScalar(mxCurrentField);
        StoppingCriterion<FFPP> newCrit2(nb_iter2);
        crit2 = newCrit2;
    }
    //mexPrintf("\n crit2 nb_it = %d\n",crit2.get_crit());

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
    //fact_side
    bool factside = false;
    if (presentFields[6])
    {
      mxCurrentField = mxGetField(prhs[0],0,"fact_side");
      factside =(bool)  mxGetScalar(mxCurrentField);
    }

   //init_lambda
   FFPP init_lambda = (FFPP) 1.0;
   if (presentFields[8])
   {
       mxCurrentField = mxGetField(prhs[0],0,"init_lambda");
       init_lambda = (FFPP) mxGetScalar(mxCurrentField);
   }

   //compute_lambda
   // bool compute_lambda = true;
   // if (presentFields[9])
   // {
       // mxCurrentField = mxGetField(prhs[0],0,"compute_lambda");
       // compute_lambda = (bool) mxGetScalar(mxCurrentField);
   // }


      ///////////// HIERARCHICAL LAUNCH ///////////////
     // creation des parametres
	 try{
		std::cout<<"avant "<<std::endl;
		Faust::Params<FFPP> params(data,nb_fact,consSS,std::vector<Faust::MatDense<FFPP> >(),crit1,crit2,isVerbose,updateway,factside,init_lambda);

	 //DisplayParams(params);
     //creation de hierarchical fact
     std::cout<<"youpi"<<std::endl;
	 HierarchicalFact<FFPP> hier_fact(params);
     hier_fact.compute_facts();


     //extraction des resultats
     FFPP lambda = hier_fact.get_lambda();

     plhs[0]=mxCreateDoubleScalar((double) lambda);

     std::vector<Faust::MatDense<FFPP> > facts;
     hier_fact.get_facts(facts);


     Faust::MatDense<FFPP> current_fact = facts[0];
     mxArray * cellFacts;
     setCellFacts(&cellFacts,facts);


     plhs[1]=cellFacts;

	}
	 catch (const std::exception& e)
	{
		 mexErrMsgTxt(e.what());
	}

}

