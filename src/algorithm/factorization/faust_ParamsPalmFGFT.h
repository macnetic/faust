#ifndef __FAUST_PARAMS_PALM_FGFT_H__
#define __FAUST_PARAMS_PALM_FGFT_H__

#include "faust_ParamsPalm.h"

using namespace Faust;

namespace Faust
{


	template<typename FPP, FDevice DEVICE, typename FPP2 = double>
		class ParamsPalmFGFT : public ParamsPalm<FPP,DEVICE,FPP2>
	{

		public:

			//ctor definitions in header because it consists mainly to call parent ctor
			ParamsPalmFGFT(const MatDense<FPP,DEVICE>& data_,
					const int nbFact_,
					const std::vector<const ConstraintGeneric*>& cons_,
					const std::vector<MatDense<FPP,DEVICE> >& init_fact_,
					const MatDense<FPP,DEVICE>& init_D,
					const StoppingCriterion<FPP2> & stop_crit_ = StoppingCriterion<FPP2>(ParamsPalm<FPP,DEVICE,FPP2>::defaultNiter),
					const bool isVerbose_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultVerbosity ,
					const bool isUpdateWayR2L_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultUpdateWayR2L ,
					const FPP2 init_lambda_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultLambda,
					const FPP2 step_size_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultStepSize, const GradientCalcOptMode gradCalcOptMode = Params<FPP,DEVICE,FPP2>::defaultGradCalcOptMode) : ParamsPalm<FPP, DEVICE, FPP2>(data_, nbFact_, cons_, init_fact_, stop_crit_, isVerbose_, isUpdateWayR2L_, init_lambda_, true /*constant_step_size is always true for Palm4MSAFGFT */, step_size_, gradCalcOptMode), init_D(init_D) {}

			ParamsPalmFGFT() : ParamsPalm<FPP,DEVICE,FPP2>(), init_D(0,0) {}

			ParamsPalmFGFT(const MatDense<FPP,DEVICE>& data_,
					const int nbFact_,
					const std::vector<const ConstraintGeneric*>& cons_,
					const std::vector<MatDense<FPP,DEVICE> >& init_fact_,
					const Vect<FPP,DEVICE>& init_D_diag,
					const StoppingCriterion<FPP2> & stop_crit_ = StoppingCriterion<FPP2>(ParamsPalm<FPP,DEVICE,FPP2>::defaultNiter),
					const bool isVerbose_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultVerbosity ,
					const bool isUpdateWayR2L_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultUpdateWayR2L ,
					const FPP2 init_lambda_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultLambda,
					const FPP2 step_size_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultStepSize,const GradientCalcOptMode gradCalcOptMode = Params<FPP,DEVICE,FPP2>::defaultGradCalcOptMode);

			ParamsPalmFGFT(const MatDense<FPP,DEVICE>& data_,
					const int nbFact_,
					const std::vector<const ConstraintGeneric*>& cons_,
					const std::vector<MatDense<FPP,DEVICE> >& init_fact_,
					const FPP* init_D_diag,
					const StoppingCriterion<FPP2> & stop_crit_ = StoppingCriterion<FPP2>(ParamsPalm<FPP,DEVICE,FPP2>::defaultNiter),
					const bool isVerbose_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultVerbosity ,
					const bool isUpdateWayR2L_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultUpdateWayR2L ,
					const FPP2 init_lambda_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultLambda,
					const FPP2 step_size_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultStepSize, const GradientCalcOptMode gradCalcOptMode = Params<FPP,DEVICE,FPP2>::defaultGradCalcOptMode) : ParamsPalmFGFT<FPP,DEVICE,FPP2>(data_, nbFact_, cons_, init_fact_, Vect<FPP,DEVICE>(data_.getNbRow(), init_D_diag), stop_crit_, isVerbose_, isUpdateWayR2L_, init_lambda_, step_size_, gradCalcOptMode) {}


			MatDense<FPP,DEVICE> init_D;


	};

#include "faust_ParamsPalmFGFT.hpp"

}
#endif
