#ifndef __FAUST_HIERARCHICAL_FACT_H__
#define __FAUST_HIERARCHICAL_FACT_H__

class faust_params;
class faust_constraint;
class palm4MSA;
class faust_real;


class hierarchical_fact
{
   public:
      hierarchical_fact() // voir avec Luc les parametres par defaut
      hierarchical_fact(const faust_params& params_);

      void init();
      void next_step();


   private:
      const std::vector<const std::vector<const faust_constraint_generic*> > cons;
      bool isUpdateWayR2L;
      bool isFactSideLeft; 
      bool isVerbose;
      int ind_fact ; //indice de factorisation (!= palm4MSA::ind_fact : indice de facteur)
      palm4MSA palm_2;
      palm4MsA palm_global;
      faust_real lambda;
      std::vector<faust_mat> S;
      
};

#endif
