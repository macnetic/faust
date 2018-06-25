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
/****************************************************************************/
#ifndef __HIERARCHICAL_FACT_CU_HPP__
#define __HIERARCHICAL_FACT_CU_HPP__

//#include "faust_HierarchicalFact.h"
#ifdef __COMPILE_TIMERS__
    #include "faust_Timer.h"
#endif

#ifdef __COMPILE_GPU__
	#include "faust_MatSparse_gpu.h"
	#include "faust_Transform_gpu.h"
#else
	#include "faust_MatSparse.h"
	#include "faust_Transform.h"
#endif

#include "faust_exception.h"
using namespace std;

//Faust::HierarchicalFact::Faust::HierarchicalFact(){} // voir avec Luc les parametres par defaut

template<typename FPP,Device DEVICE>
const char * Faust::HierarchicalFact<FPP,DEVICE>::m_className="Faust::HierarchicalFact";

template<typename FPP,Device DEVICE>
Faust::HierarchicalFact<FPP,DEVICE>::HierarchicalFact(const Faust::MatDense<FPP,DEVICE>& M,const Faust::Params<FPP,DEVICE>& params_, Faust::BlasHandle<DEVICE> cublasHandle, Faust::SpBlasHandle<DEVICE> cusparseHandle):
   m_indFact(0),
   cons(params_.cons),
   m_isUpdateWayR2L(params_.isUpdateWayR2L),
   m_isFactSideLeft(params_.isFactSideLeft),
   m_isVerbose(params_.isVerbose),
   nbFact(params_.m_nbFact-1),
   palm_2(Palm4MSA<FPP,DEVICE>(M,params_, cublasHandle, false)),
   palm_global(Palm4MSA<FPP,DEVICE>(M,params_, cublasHandle, true)),
   cons_tmp_global(vector<const Faust::ConstraintGeneric<FPP,DEVICE>*>()),
   default_lambda(params_.init_lambda),
   isFactorizationComputed(false),
   errors(std::vector<std::vector<FPP> >(2,std::vector<FPP >(params_.m_nbFact-1,0.0))),
   cublas_handle(cublasHandle),
   cusparse_handle(cusparseHandle)
{
	   // check if the params and M are compatible
   if ((M.getNbRow() != params_.m_nbRow) |  (M.getNbCol() != params_.m_nbCol))
		handleError(m_className,"constructor : params and matrix haven't compatible size"); 
}


template<typename FPP,Device DEVICE>
void Faust::HierarchicalFact<FPP,DEVICE>::init()
{
#ifdef __COMPILE_TIMERS__
t_init.start();
#endif
   
  

   cons_tmp_global.clear();
   if(m_isFactSideLeft)
      cons_tmp_global.push_back(cons[0][m_indFact]);
   else
      cons_tmp_global.push_back(cons[1][m_indFact]);


   palm_global.set_constraint(cons_tmp_global);
   palm_global.init_fact(1);


#ifdef __COMPILE_TIMERS__
t_init.stop();
#endif

}





void check_validity()
{















}















template<typename FPP,Device DEVICE>
void Faust::HierarchicalFact<FPP,DEVICE>::next_step()
{
#ifdef __COMPILE_TIMERS__
t_next_step.start();
#endif

   if(isFactorizationComputed)
   {
      handleError(m_className,"next_step : factorization has already been computed");
   }

   vector<const Faust::ConstraintGeneric<FPP,DEVICE>*> cons_tmp_2(2);
   cons_tmp_2[0]=cons[0][m_indFact];
   cons_tmp_2[1]=cons[1][m_indFact];


   palm_2.set_constraint(cons_tmp_2);

   palm_2.init_fact(2);

   palm_2.set_lambda(default_lambda);

#ifdef __COMPILE_TIMERS__
palm_2.init_local_timers();
#endif
   //while(palm_2.do_continue())
    //  palm_2.next_step();
	palm_2.compute_facts();


#ifdef __COMPILE_TIMERS__
palm_2.print_local_timers();
#endif
   palm_global.update_lambda_from_palm(palm_2);


   if (m_isFactSideLeft)
   {
      cons_tmp_global[0]=cons[0][m_indFact];
      typename vector<const Faust::ConstraintGeneric<FPP,DEVICE>*>::iterator it;
      it = cons_tmp_global.begin();
      cons_tmp_global.insert(it+1,cons[1][m_indFact]);
   }
   else
   {
      typename vector<const Faust::ConstraintGeneric<FPP,DEVICE>*>::iterator it;
      it = cons_tmp_global.begin();
      cons_tmp_global.insert(it+m_indFact,cons[0][m_indFact]);
      cons_tmp_global[m_indFact+1]=cons[1][m_indFact];
   }

   palm_global.set_constraint(cons_tmp_global);


   palm_global.init_fact_from_palm(palm_2, m_isFactSideLeft);

#ifdef __COMPILE_TIMERS__
palm_global.init_local_timers();
#endif
   //while(palm_global.do_continue())
    //  palm_global.next_step();
	palm_global.compute_facts();
#ifdef __COMPILE_TIMERS__
palm_global.print_local_timers();
#endif

   palm_2.set_data(palm_global.get_res(m_isFactSideLeft, m_indFact));


   compute_errors();


   m_indFact++;


#ifdef __COMPILE_TIMERS__
t_next_step.stop();
palm_2.print_prox_timers();
#endif
}

template<typename FPP,Device DEVICE>
void Faust::HierarchicalFact<FPP,DEVICE>::get_facts(Faust::Transform<FPP,DEVICE> & fact)const
{
	std::vector<Faust::MatSparse<FPP,DEVICE> > spfacts;
	get_facts(spfacts);
	Faust::Transform<FPP,DEVICE> res(spfacts);
	fact = res;
}

template<typename FPP,Device DEVICE>
void Faust::HierarchicalFact<FPP,DEVICE>::get_facts(std::vector<Faust::MatSparse<FPP,DEVICE> >& sparse_facts)const
{
   /*if(!isFactorizationComputed)
   {
      cerr << "Error in Faust::HierarchicalFact<FPP,DEVICE>::get_facts : factorization has not been computed" << endl;
      exit(EXIT_FAILURE);
   }*/

   const std::vector<Faust::MatDense<FPP,DEVICE> >& full_facts = palm_global.get_facts();
   sparse_facts.resize(full_facts.size());
   for (int i=0 ; i<sparse_facts.size() ; i++)
{
      sparse_facts[i].init(full_facts[i],cusparse_handle);
}
}



template<typename FPP,Device DEVICE>
void Faust::HierarchicalFact<FPP,DEVICE>::compute_facts()
{
   if(isFactorizationComputed)
   {
      handleError(m_className,"compute_facts : factorization has already been computed");
   }

  init();
  for (int i=0 ; i<=nbFact-1 ; i++)
  {
     cout << "Faust::HierarchicalFact<FPP,DEVICE>::compute_facts : factorisation "<<i+1<<"/"<<nbFact <<endl;
     next_step();
  }

  isFactorizationComputed = true;

}


template<typename FPP,Device DEVICE>
const std::vector<std::vector< FPP> >& Faust::HierarchicalFact<FPP,DEVICE>::get_errors()const
{
    if(!isFactorizationComputed)
    {
        handleError(m_className,"get_errors() : Factorization has not been computed");
    }
    return errors;
}

template<typename FPP,Device DEVICE>
void Faust::HierarchicalFact<FPP,DEVICE>::compute_errors()
{
   vector<Faust::MatSparse<FPP,DEVICE> > sp_facts;
   get_facts(sp_facts);
   int nb_factor = sp_facts.size();	
   vector< Faust::MatGeneric<FPP,DEVICE> *> Transform_facts;
   Transform_facts.resize(sp_facts.size());

   for (int i=0; i < nb_factor;i++)
   	Transform_facts[i] = sp_facts[i].Clone();   	

   Faust::Transform<FPP,DEVICE> faust_Transform_tmp(Transform_facts, get_lambda());
   
   // delete all the dynamic memory allocated by the Clone method
   for (int i=0; i < nb_factor;i++)
	delete 	Transform_facts[i];
  
   const Faust::MatDense<FPP,DEVICE> estimate_mat = faust_Transform_tmp.get_product(cublas_handle, cusparse_handle);
   Faust::MatDense<FPP,DEVICE> data(palm_global.get_data());
   	
   FPP data_norm = data.norm();

   
   data -= estimate_mat;
   	
   errors[0][m_indFact] =  estimate_mat.norm()/data_norm;
   errors[1][m_indFact] =  faust_Transform_tmp.get_total_nnz()/data.getNbRow()/data.getNbCol();

}


#ifdef __COMPILE_TIMERS__
template<typename FPP,Device DEVICE> Faust::Timer Faust::HierarchicalFact<FPP,DEVICE>::t_init;
template<typename FPP,Device DEVICE> Faust::Timer Faust::HierarchicalFact<FPP,DEVICE>::t_next_step;

template<typename FPP,Device DEVICE>
void Faust::HierarchicalFact<FPP,DEVICE>::print_timers()const
{
   palm_global.print_global_timers();
   cout << "timers in Faust::HierarchicalFact :" << endl;
   cout << "t_init      = " << t_init.get_time()      << " s for "<< t_init.get_nb_call()      << " calls" << endl;
   cout << "t_next_step = " << t_next_step.get_time() << " s for "<< t_next_step.get_nb_call() << " calls" << endl<<endl;



}
#endif

#endif
