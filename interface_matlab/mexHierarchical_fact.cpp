#include <mex.h>
//#include <mexutils.h>
#include <hierarchical_fact.h>
#include <faust_mat.h>
#include <faust_constraint_int.h>
#include <faust_constraint_generic.h>
#include <faust_constraint_real.h>
#include <faust_constraint_mat.h>
#include <vector>
#include <string>
#include <algorithm>
#include <faust_params.h>
#include <faust_constant.h>

#include <mexFaustMat.h>

void getConstraint(std::vector<const faust_constraint_generic*> & consS,mxArray* mxCons);
void setCellFacts(mxArray ** cellFacts,std::vector<faust_mat> facts);
void testCoherence(const mxArray* params,std::vector<bool> & presentFields);
void DisplayParams(const faust_params & params);



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
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
     mexPrintf(" NUMBER FIELDS %d\n",presentFields.size());
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
    faust_mat data;
    if (presentFields[0])
    {    
        mxCurrentField = mxGetField(prhs[0],0,"data");  
        
        data = getFaustMat(  mxCurrentField ) ;   
        mexPrintf("DATA");
        for (int i = 0;i<data.getNbRow();i++)
        {
            for (int j = 0;j<data.getNbCol();j++)
            {
                //bidon = std::snprintf(coeff,10,"%d",A(i,j));	
                mexPrintf("%f ",data(i,j));
            }
            mexPrintf("\n");
        }
    }else
    {
         mexErrMsgTxt("params.data must be specified");
    }
    
   //nb_fact initialisation
   int nb_fact; 
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
    std::vector<std::vector<const faust_constraint_generic*> > consSS;
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
        
        if(nbRowCons !=2)
        {
            mexPrintf("\n cons has %d rows \n",nbRowCons);
            mexErrMsgTxt("cons must have 2 rows");
        }
        if(nbColCons != (nb_fact-1))
        {
            mexPrintf("\n cons has %d cols and nb_fact = %d\n",nbColCons,nb_fact);
            mexErrMsgTxt("incoherence between the number of columns of cons and nfacts ");
        }
        mexPrintf("\n cons has %d rows and %d cols \n",nbRowCons,nbColCons);
        faust_constraint_generic * consToAdd;
        std::vector<const faust_constraint_generic*> consS;
        
        for (int i=0;i<nbRowCons;i++)
        {
            mexPrintf("%d / %d\n",i,nbRowCons);
            for (int j=0;j<nbColCons;j++)
            {
                mexPrintf("cons(%d , %d)\n",i,j);
                mxCurrentCons=mxGetCell(mxCurrentField,i+(j*nbRowCons));
                getConstraint(consS,mxCurrentCons);
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
    stopping_criterion crit1;
    if (presentFields[3])
    {   
         mxCurrentField = mxGetField(prhs[0],0,"niter1");
        int nb_iter1 =(int)  mxGetScalar(mxCurrentField);
        stopping_criterion newCrit1(nb_iter1);
        crit1 = newCrit1;
    }
    mexPrintf("\n crit1 nb_it = %d\n",crit1.get_crit());
    
    //niter2
    stopping_criterion crit2;
    if (presentFields[4])
    {   
         mxCurrentField = mxGetField(prhs[0],0,"niter2");
        int nb_iter2 =(int)  mxGetScalar(mxCurrentField);
        stopping_criterion newCrit2(nb_iter2);
        crit2 = newCrit2;
    }
    mexPrintf("\n crit2 nb_it = %d\n",crit2.get_crit());
    
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
   double init_lambda = 1.0;
   if (presentFields[8])
   {
       mxCurrentField = mxGetField(prhs[0],0,"init_lambda");
       init_lambda = (faust_real) mxGetScalar(mxCurrentField);
   }
   
   //compute_lambda
   bool compute_lambda = true;
   if (presentFields[9])
   {
       mxCurrentField = mxGetField(prhs[0],0,"compute_lambda");
       compute_lambda = (bool) mxGetScalar(mxCurrentField);
   }
    
    
      ///////////// HIERARCHICAL LAUNCH ///////////////  
     // creation des parametres   
     faust_params params(data,nb_fact,consSS,std::vector<faust_mat>(),crit1,crit2,isVerbose,updateway,factside,init_lambda,compute_lambda);   
     DisplayParams(params);
     //creation de hierarchical fact
     hierarchical_fact hier_fact(params);
     hier_fact.compute_facts();	
     
     
     //extraction des resultats
     faust_real lambda = hier_fact.get_lambda();
     
     plhs[0]=mxCreateDoubleScalar((double) lambda);
     
     std::vector<faust_mat> facts;
     hier_fact.get_facts(facts);

     
     faust_mat current_fact = facts[0];
     mxArray * cellFacts;
     setCellFacts(&cellFacts,facts);
     

     plhs[1]=cellFacts;

     
   
   
   
   
   

}


void testCoherence(const mxArray* params,std::vector<bool> & presentFields)
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
        if (strcmp(fieldName,"compute_lambda") == 0)
        {
            presentFields[9] = true;
        }
  }
  
}



void setCellFacts(mxArray **  cellFacts,std::vector<faust_mat> facts)
{   
    int rowFact,colFact;
    int nb_fact = facts.size();
    faust_mat mat;
    //mexPrintf("size CellFacts nb fact %d ",nb_fact);
    (*cellFacts) = mxCreateCellMatrix(1,nb_fact);
    mxArray * mxMat;
    double* mat_ptr;
    mxMat = mxCreateDoubleMatrix(0,0,mxREAL);
    
    for (size_t k = 0; k < nb_fact; k++)
    {
        mat = facts[k];
        rowFact = mat.getNbRow();
        colFact = mat.getNbCol();
        
        mat_ptr = (double *) mxCalloc(rowFact*colFact,sizeof(double));
        //mat_ptr = mxGetPr(mxMat);
        memcpy(mat_ptr,mat.getData(),rowFact*colFact*sizeof(double));
        mxSetM(mxMat, rowFact);
        mxSetN(mxMat, colFact);
        mxSetPr(mxMat, mat_ptr);
        
//          mexPrintf("DIM : %d %d\n",rowFact,colFact);
//         for (int i=1;i<rowFact;i++)
//         {
//             for (int j=1;j<colFact;j++)
//             {
//                 mexPrintf("%f ",mat_ptr[i+j*rowFact]);
//             }
//             mexPrintf("\n");
//         }
//          mexPrintf("\n");
//          mexPrintf("\n");
         //rhs[0]=mxMat;
         //mexCallMATLAB(0,NULL,1,rhs, "disp");
        mxSetCell((*cellFacts),k,mxDuplicateArray(mxMat));
    }
    
//     mexPrintf("Display Cell \n");
//     for (size_t l = 0; l < mxGetNumberOfElements((*cellFacts)); l++)
//     {    
//         mexPrintf("facts %d \n",l);
//          mexPrintf("\n");
//          mxArray * rhs[1]; 
//          rhs[0]=mxGetCell((*cellFacts),l);
//          mexCallMATLAB(0,NULL,1,rhs, "disp");
//     }
    
}





void getConstraint(std::vector<const faust_constraint_generic*> & consS,mxArray* mxCons)
{     
    mwSize bufCharLen,nbRowCons,nbColCons;
    int status;        
    char * consName;
    double paramCons;
    mxArray * mxConsParams; 
    
    mxConsParams=mxGetCell(mxCons,0);
    bufCharLen = mxGetNumberOfElements(mxConsParams)+1;
    consName = (char *) mxCalloc(bufCharLen,sizeof(char));
    status = mxGetString(mxConsParams,consName,bufCharLen);
    paramCons = mxGetScalar(mxConsParams);
    mxConsParams=mxGetCell(mxCons,2);
    nbRowCons   = (int) mxGetScalar(mxConsParams);
    mxConsParams=mxGetCell(mxCons,3);
    nbColCons = (int) mxGetScalar(mxConsParams);
    
    bool is_const_int =((strcmp(consName,"sp") == 0) || (strcmp(consName,"sppos")==0));
	is_const_int = ((is_const_int) || ((strcmp(consName,"spcol") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"splin") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"splincol") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"lOpen") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"l1pen") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"splin") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"wav") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"blkdiag") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"splin_test") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"supp") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(consName,"normlin") == 0)));
    
    int const_type = -1;
	if (is_const_int)
	{
		const_type = 0;
	}
	/*if (is_const_real)
	{
		const_type = 1;
	}
	if (is_const_mat)
	{
		const_type = 2;
	}
     */
    int intParameter;
    switch(const_type)
	{
        case 0:
            mxConsParams=mxGetCell(mxCons,1);
             intParameter = (int) std::floor(mxGetScalar(mxConsParams)+0.5);
             mexPrintf("NAME  %s PARAMS %d DIMS : (%d,%d)\n",consName,intParameter,nbRowCons,nbColCons);
            faust_constraint_name consNameType;
            if (strcmp(consName,"sp") == 0)
            {
                consNameType = CONSTRAINT_NAME_SP;
            }
            if (strcmp(consName,"spcol") == 0)
            {
                consNameType = CONSTRAINT_NAME_SPCOL;
            }
            if (strcmp(consName,"splin") == 0)
            {
                consNameType = CONSTRAINT_NAME_SPLIN;
            }
            
            consS.push_back((new faust_constraint_int(consNameType,intParameter,nbRowCons,nbColCons)));
            break;
            
            
            default :
                mexErrMsgTxt("Unknown constraint name ");
			break;
    }

   
   
    
}


void DisplayParams(const faust_params & params)
{   
    mexPrintf("/////// PARAMS //////\n");
    faust_mat data = params.data; 
   
    int nbRow = data.getNbRow();
    int nbCol = data.getNbCol();
    mexPrintf("DATA DIM : %d %d\n",nbRow,nbCol);
    if (nbRow*nbCol < 1000)
    {
        for (int i = 0;i<data.getNbRow();i++)
        {
            for (int j = 0;j<data.getNbCol();j++)
            {
                mexPrintf("%f ",data(i,j));
            }
            mexPrintf("\n");
        }
    }
    mexPrintf("\n\n");
    
    mexPrintf("NB_FACTS : %d\n",params.nb_fact);
    mexPrintf("CONSTRAINTS : nbr %d\n",params.cons[0].size());
    for (int L=0;L<params.cons[0].size();L++)
	{
		for (int jl=0;jl<params.cons.size();jl++)
		{	
			//char * type_cons =  getType((*params.cons[jl][L]).get_constraint_name());
            char * type_cons =  (*params.cons[jl][L]).getType();
			mexPrintf(" %s ",(*params.cons[jl][L]).get_constraint_name());
			mexPrintf(" DIMS (%d,%d) : ",(*params.cons[jl][L]).getRows(),(*params.cons[jl][L]).getCols());
			
			
			if (strcmp(type_cons,"INT") == 0)
			{	
				faust_constraint_int* const_int = (faust_constraint_int*)(params.cons[jl][L]);
				mexPrintf(" parameter : %d",(*const_int).getParameter());
			}
			
			if (strcmp(type_cons,"FAUST_REAL") == 0)
			{	
				faust_constraint_real* const_real = (faust_constraint_real*)(params.cons[jl][L]);
				mexPrintf(" parameter : %f",(*const_real).getParameter());
			}
			mexPrintf("\n"); 

			
		}
		
	}
}