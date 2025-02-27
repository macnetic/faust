/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2019):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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
#ifndef __FAUST_PARAMS_PALM_H__
#define __FAUST_PARAMS_PALM_H__

#include "faust_constant.h"
#include <vector>
#ifdef __COMPILE_GPU__
   #include "faust_MatDense_gpu.h"
#else
   #include "faust_MatDense.h"
#endif
#include "faust_StoppingCriterion.h"
#include "faust_ConstraintGeneric.h"
#include "faust_Params.h"


/*! \class Faust::ParamsPalm
  * \brief template class representing the parameters for building HierarchicalFact object
  * \param FPP scalar numeric type, e.g float or double
*/


//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{
    template<typename FPP,FDevice DEVICE> class MatDense;

    template<typename FPP,FDevice DEVICE,typename FPP2 = double>
    class ParamsPalm
    {


       public:

          ParamsPalm(
             const MatDense<FPP,DEVICE>& data_,
             const int nbFact_,
             const std::vector<const ConstraintGeneric*>& cons_,
             const std::vector<MatDense<FPP,DEVICE> >& init_fact_,
             const StoppingCriterion<FPP2> & stop_crit_ = StoppingCriterion<FPP2>(defaultNiter),
             const bool isVerbose_ = defaultVerbosity ,
             const bool isUpdateWayR2L_ = defaultUpdateWayR2L ,
             const FPP2 init_lambda_ = defaultLambda,
             const bool constant_step_size_ = defaultConstantStepSize,
             const FPP2 step_size_ = defaultStepSize,
			 const GradientCalcOptMode gradCalcOptMode = Params<FPP,DEVICE,FPP2>::defaultGradCalcOptMode,
			 const bool use_MHTP = Params<FPP,DEVICE, FPP2>::defaultUseMHTP,
			 const bool no_normalization = Params<FPP, DEVICE, FPP2>::defaultNoNormalization,
			 const bool no_lambda = Params<FPP,DEVICE, FPP2>::defaultNoLambda);

          void check_constraint_validity();
           ParamsPalm();
          ~ParamsPalm(){}


       public:
          // Required members
          MatDense<FPP,DEVICE> data;
          int nbFact; // number of factors
          std::vector<const ConstraintGeneric*> cons; // vector of constraints

          // Optional members (set to default values if not defined)
          std::vector<MatDense<FPP,DEVICE> > init_fact;
          StoppingCriterion<FPP2> stop_crit;
          bool isVerbose;
          bool isUpdateWayR2L;
          bool isConstantStepSize;
          FPP2 step_size;
          FPP2 init_lambda;
		  GradientCalcOptMode gradCalcOptMode;
          Real<FPP> norm2_threshold;
          unsigned int norm2_max_iter;
		  FactorsFormat factors_format;
		  bool packing_RL;
		  bool no_normalization;
		  bool no_lambda;
		  bool use_MHTP;
          StoppingCriterion<Real<FPP>> stop_crit_MHTP;

          void Display() const;
          void init_factors();
          static const int defaultNiter;
          static const bool defaultVerbosity;
          static const bool defaultUpdateWayR2L;
          static const FPP2 defaultLambda;
          static const bool defaultConstantStepSize;
          static const FPP2 defaultStepSize;
          private :
          static const char *  m_className;


          /*const int nb_it;   // number of iterations
          // if isStoppingCriterionError then criterion is error else criterion is number of iteration
          bool  isStoppingCriterionError;
          const FPP errorThreshold;
          // only used as stopping criterion, if isStoppingCriterionError, when error is still greater than
          int maxIteration;*/
    };

}

#include "faust_ParamsPalm.hpp"


#endif
