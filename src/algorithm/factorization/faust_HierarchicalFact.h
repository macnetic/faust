/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2018):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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
/*      Hakim H. hakim.hadj-djilani@inria.fr                                */
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
#ifndef __FAUST_HIERARCHICAL_FACT_H__
#define __FAUST_HIERARCHICAL_FACT_H__

#include "faust_constant.h"
#include <vector>
#include "faust_Palm4MSA.h"
#include "faust_Params.h"

#ifdef __COMPILE_FPPIMERS__
  #include "faust_Timer.h"
#endif

/*! \class HierarchicalFact
   * \brief template class implementing hierarchical factorisation algorithm
    : <br>  factorization of a data matrix into multiple factors using Faust::Palm4MSA algorithm mixed with a hierarchical approach.
   *\tparam FPP scalar numeric type, e.g float or double
   */


namespace Faust
{

    class ConstraintGeneric;
    template<typename FPP,Device DEVICE,typename FPP2> class Palm4MSA;
    template<typename FPP,Device DEVICE> class MatDense;
    template<typename FPP,Device DEVICE> class MatSparse;
    template<typename FPP,Device DEVICE> class Transform;

    template<typename FPP2> class StoppingCriterion;
//    template<Device DEVICE> class BlasHandle;
//    template<Device DEVICE> class SpBlasHandle;

    template<typename FPP,Device DEVICE,typename FPP2 = double>
    class HierarchicalFact
    {

       public:

          HierarchicalFact(const Faust::MatDense<FPP,DEVICE>& M, const Faust::Params<FPP,DEVICE,FPP2>& params_, Faust::BlasHandle<DEVICE> cublasHandle, SpBlasHandle<DEVICE> cusparseHandle);
          void get_facts(Faust::Transform<FPP,DEVICE> &)const;
          void get_facts(std::vector<Faust::MatSparse<FPP,DEVICE> >&)const;
          void get_facts(std::vector<Faust::MatDense<FPP,DEVICE> >& fact)const{fact = palm_global.get_facts();}
          void compute_facts();
          FPP get_lambda()const{return palm_global.get_lambda();}
          const std::vector<std::vector< FPP> >& get_errors()const;


        private:
          void init();
          void next_step();
		protected:
          void compute_errors();



        protected:
          const std::vector< std::vector<const Faust::ConstraintGeneric*>> cons;
          bool m_isUpdateWayR2L;
          bool m_isFactSideLeft;
          bool m_isVerbose;
          int m_indFact ; //indice de factorisation (!= Faust::Palm4MSA::m_indFact : indice de facteur)
          int nbFact; // nombre de factorisations (!= Faust::Palm4MSA::nbFact : nombre de facteurs)
          Faust::Palm4MSA<FPP,DEVICE,FPP2> palm_2;
          Faust::Palm4MSA<FPP,DEVICE,FPP2> palm_global;
          const FPP default_lambda; // initial value of lambda for factorization into two factors
          //std::vector<Faust::MatDense<FPP,DEVICE> > S;
		  std::vector<const Faust::ConstraintGeneric*> cons_tmp_global;
          bool isFactorizationComputed;
          std::vector<std::vector<FPP> > errors;
          static const char * m_className;
          Faust::BlasHandle<DEVICE> cublas_handle;
          Faust::SpBlasHandle<DEVICE> cusparse_handle;


    #ifdef __COMPILE_TIMERS__
       public:
          static Faust::Timer t_init;
          static Faust::Timer t_next_step;

        void print_timers()const;
        //void print_prox_timers()const;
    #endif



    };

}

#include "faust_HierarchicalFact.hpp"

#endif
