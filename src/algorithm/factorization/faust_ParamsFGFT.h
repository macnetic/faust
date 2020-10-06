
#ifndef __FAUST_PARAMSFGFT_H__
#define __FAUST_PARAMSFGFT_H__

#include "faust_Params.h"

using namespace Faust;

namespace Faust
{
	template<typename FPP, FDevice DEVICE, typename FPP2 = double>
	class ParamsFGFT : public Params<FPP, DEVICE, FPP2>
	{

		public:
		MatDense<FPP, DEVICE> init_D; //TODO: convert to Sparse or Diag repres. and set private or protected
		//TODO: does it really need to be public 
		//TODO: move the ctor def into .hpp
		ParamsFGFT(
				const faust_unsigned_int nbRow,
				const faust_unsigned_int nbCol,
				const unsigned int nbFact,
				const std::vector<std::vector<const ConstraintGeneric*>> & cons,
				const std::vector<MatDense<FPP,DEVICE> >& init_fact,
				const MatDense<FPP, DEVICE>& init_D,
				const StoppingCriterion<FPP2>& stop_crit_2facts = StoppingCriterion<FPP2>(Params<FPP,DEVICE,FPP2>::defaultNiter1),
				const StoppingCriterion<FPP2>& stop_crit_global = StoppingCriterion<FPP2>(Params<FPP,DEVICE,FPP2>::defaultNiter2),
				const bool isVerbose = Params<FPP,DEVICE,FPP2>::defaultVerbosity,
				const bool isUpdateWayR2L = Params<FPP,DEVICE,FPP2>::defaultUpdateWayR2L,
				const bool isFactSideLeft = Params<FPP,DEVICE,FPP2>::defaultFactSideLeft,
				const FPP2 init_lambda = Params<FPP,DEVICE,FPP2>::defaultLambda,
				const bool constant_step_size = Params<FPP,DEVICE,FPP2>::defaultConstantStepSize,
				const FPP2 step_size = Params<FPP,DEVICE,FPP2>::defaultStepSize, const GradientCalcOptMode gradCalcOptMode = Params<FPP,DEVICE,FPP2>::defaultGradCalcOptMode): Params<FPP, DEVICE, FPP2>(nbRow, nbCol, nbFact, cons, init_fact, stop_crit_2facts, stop_crit_global, isVerbose, isUpdateWayR2L, isFactSideLeft, init_lambda, constant_step_size, step_size, gradCalcOptMode), init_D(init_D)
		{

		}

		ParamsFGFT(
				const faust_unsigned_int nbRow,
				const faust_unsigned_int nbCol,
				const unsigned int nbFact,
				const std::vector<std::vector<const ConstraintGeneric*>> & cons,
				const std::vector<MatDense<FPP,DEVICE> >& init_fact,
				const Vect<FPP, DEVICE>& init_D_diag,
				const StoppingCriterion<FPP2>& stop_crit_2facts = StoppingCriterion<FPP2>(Params<FPP,DEVICE,FPP2>::defaultNiter1),
				const StoppingCriterion<FPP2>& stop_crit_global = StoppingCriterion<FPP2>(Params<FPP,DEVICE,FPP2>::defaultNiter2),
				const bool isVerbose = Params<FPP,DEVICE,FPP2>::defaultVerbosity,
				const bool isUpdateWayR2L = Params<FPP,DEVICE,FPP2>::defaultUpdateWayR2L,
				const bool isFactSideLeft = Params<FPP,DEVICE,FPP2>::defaultFactSideLeft,
				const FPP2 init_lambda = Params<FPP,DEVICE,FPP2>::defaultLambda,
				const bool constant_step_size = Params<FPP,DEVICE,FPP2>::defaultConstantStepSize,
				const FPP2 step_size = Params<FPP,DEVICE,FPP2>::defaultStepSize, const GradientCalcOptMode gradCalcOptMode = Params<FPP,DEVICE,FPP2>::defaultGradCalcOptMode): Params<FPP, DEVICE, FPP2>(nbRow, nbCol, nbFact, cons, init_fact, stop_crit_2facts, stop_crit_global, isVerbose, isUpdateWayR2L, isFactSideLeft, init_lambda, constant_step_size, step_size, gradCalcOptMode), init_D(nbRow, nbCol)
		{
			init_D.setZeros();
			// set init_D from diagonal vector init_D_diag
			for(int i=0;i<nbRow;i++)
				init_D.getData()[i*nbRow+i] = init_D_diag.getData()[i];
		}

		ParamsFGFT(
				const faust_unsigned_int nbRow,
				const faust_unsigned_int nbCol,
				const unsigned int nbFact,
				const std::vector<std::vector<const ConstraintGeneric*>> & cons,
				const std::vector<MatDense<FPP,DEVICE> >& init_fact,
				const FPP* init_D_diag,
				const StoppingCriterion<FPP2>& stop_crit_2facts = StoppingCriterion<FPP2>(Params<FPP,DEVICE,FPP2>::defaultNiter1),
				const StoppingCriterion<FPP2>& stop_crit_global = StoppingCriterion<FPP2>(Params<FPP,DEVICE,FPP2>::defaultNiter2),
				const bool isVerbose = Params<FPP,DEVICE,FPP2>::defaultVerbosity,
				const bool isUpdateWayR2L = Params<FPP,DEVICE,FPP2>::defaultUpdateWayR2L,
				const bool isFactSideLeft = Params<FPP,DEVICE,FPP2>::defaultFactSideLeft,
				const FPP2 init_lambda = Params<FPP,DEVICE,FPP2>::defaultLambda,
				const bool constant_step_size = Params<FPP,DEVICE,FPP2>::defaultConstantStepSize,
				const FPP2 step_size = Params<FPP,DEVICE,FPP2>::defaultStepSize, const GradientCalcOptMode gradCalcOptMode = Params<FPP,DEVICE,FPP2>::defaultGradCalcOptMode): ParamsFGFT<FPP, DEVICE, FPP2>(nbRow, nbCol, nbFact, cons, init_fact, Vect<FPP,DEVICE>(nbRow, init_D_diag), stop_crit_2facts, stop_crit_global, isVerbose, isUpdateWayR2L, isFactSideLeft, init_lambda, constant_step_size, step_size, gradCalcOptMode)
		{

		}



		ParamsFGFT() {}

		void Display() const;


	};
}

#include "faust_ParamsFGFT.hpp"

#endif
