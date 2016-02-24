#include "mex.h"
//#include "mexutils.h"
#include "faust_mat.h"
#include "faust_constraint_int.h"
#include "faust_constraint_generic.h"
#include "faust_constraint_real.h"
#include "faust_constraint_mat.h"
#include <vector>
#include <string>
#include <algorithm>
#include "faust_params_palm.h"
#include "faust_constant.h"
#include "palm4MSA.h"
#include <stdexcept>
#include "tools_mex.h"


void testCoherencePalm4MSA(const mxArray* params,std::vector<bool> & presentFields);



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	#ifdef FAUST_VERBOSE
		if (typeid(faust_real) == typeid(float))
		{
			std::cout<<"faust_real == float"<<std::endl;
		}
	
		if (typeid(faust_real) == typeid(double))
		{
			std::cout<<"faust_real == double"<<std::endl;
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
     mexPrintf(" NUMBER FIELDS %d\n",presentFields.size());
    
    
    
    
    mxArray    *mxCurrentField,*mxCurrentCons;
    
    
    // data initialisation
    faust_mat<faust_real> data;
    if (presentFields[0])
    {    
        mxCurrentField = mxGetField(prhs[0],0,"data");  
        
        getFaustMat(  mxCurrentField,data ) ;   
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
    
   //nb_fact initialisation
   int nb_fact=0; 
   if (presentFields[1])
   {    
        
        mxCurrentField = mxGetField(prhs[0],0,"nfacts");
        mexPrintf("a\n");
        nb_fact =(int)  mxGetScalar(mxCurrentField);
        mexPrintf("NB FACT : %d\n",nb_fact);
   }else
   {
        mexErrMsgTxt("params.nfacts must be specified");
   }
   
   

    //constraints
    std::vector<const faust_constraint_generic*> consS;
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
        if(nbColCons != (nb_fact))
        {
            //mexPrintf("\n cons has %d cols and nb_fact = %d\n",nbColCons,nb_fact);
            //mexErrMsgTxt("incoherence between the number of columns of cons and nfacts ");
        }
        //mexPrintf("\n cons has %d rows and %d cols \n",nbRowCons,nbColCons);
        //faust_constraint_generic * consToAdd;
        
        
        
            for (mwSize j=0;j<nbColCons;j++)
            {
                mexPrintf("cons(%d)\n",j);
                mxCurrentCons=mxGetCell(mxCurrentField,j);
                getConstraint<faust_real>(consS,mxCurrentCons);
                //consS.push_back(consToAdd);
            }

    }else
    {
        mexErrMsgTxt("params.cons must be specified");
    } 
    std::cout<<"FINI_CONS"<<std::endl;
    //niter1
    stopping_criterion<faust_real> crit1;
    if (presentFields[3])
    {   
         mxCurrentField = mxGetField(prhs[0],0,"niter");
        int nb_iter1 =(int)  mxGetScalar(mxCurrentField);
        stopping_criterion<faust_real> newCrit1(nb_iter1);
        crit1 = newCrit1;
    }
    mexPrintf("\n crit1 nb_it = %d\n",crit1.get_crit());
    
    //init_facts
    std::vector<faust_mat<faust_real> > init_facts;
    if (presentFields[4])
    {   
         mxCurrentField = mxGetField(prhs[0],0,"init_facts");
		 std::cout<<"PASSERbeforeInitFact"<<std::endl;
		setVectorFaustMat(init_facts,mxCurrentField);
		 std::cout<<"PASSERafterInitFact"<<std::endl;
		
    }else
	{
		mexErrMsgTxt("init_facts must be must be specified");
	}

	std::cout<<"PASSER1"<<std::endl;
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
   faust_real init_lambda = 1.0;
   if (presentFields[6])
   {
       mxCurrentField = mxGetField(prhs[0],0,"init_lambda");
       init_lambda = (faust_real) mxGetScalar(mxCurrentField);
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
    faust_params_palm<faust_real> params(data,nb_fact,consS,init_facts,crit1,isVerbose,updateway,init_lambda); 
	std::cout<<"PASSER3"<<std::endl;
	palm4MSA<faust_real> palm(params,false);	
	palm.compute_facts();	
     
     
     //extraction des resultats
     faust_real lambda = palm.get_lambda();
     
     plhs[0]=mxCreateDoubleScalar((double) lambda);
     
     std::vector<faust_mat<faust_real> > facts;
     facts=palm.get_facts();

     
     faust_mat<faust_real> current_fact = facts[0];
     mxArray * cellFacts;
     setCellFacts(&cellFacts,facts);
     

     plhs[1]=cellFacts;
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
        if (strcmp(fieldName,"niter") == 0)
        {
            presentFields[3] = true;
        }
        if (strcmp(fieldName,"init_facts") == 0)
        {
            presentFields[4] = true;
        }
        if (strcmp(fieldName,"verbose") == 0)
        {
            presentFields[5] = true;
        }
        if (strcmp(fieldName,"init_lambda") == 0)
        {
            presentFields[6] = true;
        }
        if (strcmp(fieldName,"update_way") == 0)
        {
            presentFields[7] = true;
        }

  }
  
}













