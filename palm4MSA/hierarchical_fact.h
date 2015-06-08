#ifndef __FAUST_HIERARCHICAL_FACT_H__
#define __FAUST_HIERARCHICAL_FACT_H__


class hierarchical_fact
{
   public:
      hierarchical_fact() // voir avec Luc les parametres par defaut
      hierarchical_fact(const faust_params& params_);

      void init();
      void next_step();


   private:
      const vector<vector<faust_constraint> > cons;
      bool isUpdateWayR2L;
      bool isFactSideLeft; 
      bool isVerbose;
      int ind_fact ; //indice de factorisation (!= palm4MSA::ind_fact : indice de facteur)
      palm4MSA palm_2;
      palm4MsA palm_global;
      faust_real lambda;
      vector<faust_mat> S;
      
}

#endif
