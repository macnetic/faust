
#ifndef __FAUST_PARAMSFFT_H__
#define __FAUST_PARAMSFFT_H__

#include "faust_Params.h"

using namespace Faust;

namespace Faust
{
	template<typename FPP, Device DEVICE, typename FPP2 = double>
	class ParamsFFT : public Params<FPP, DEVICE, FPP2>
	{

		public:
		MatDense<FPP, DEVICE> init_D; //TODO: convert to Sparse or Diag repres.
		//TODO: does it really need to be public 
		//TODO: move the ctor def into .hpp
		ParamsFFT(
				const faust_unsigned_int nbRow,
				const faust_unsigned_int nbCol,
				const unsigned int nbFact,
				const std::vector<std::vector<const Faust::ConstraintGeneric*>> & cons,
				const std::vector<Faust::MatDense<FPP,DEVICE> >& init_fact,
				const MatDense<FPP, DEVICE>& init_D,
				const Faust::StoppingCriterion<FPP2>& stop_crit_2facts = StoppingCriterion<FPP2>(Params<FPP,DEVICE,FPP2>::defaultNiter1),
				const Faust::StoppingCriterion<FPP2>& stop_crit_global = StoppingCriterion<FPP2>(Params<FPP,DEVICE,FPP2>::defaultNiter2),
				const bool isVerbose = Params<FPP,DEVICE,FPP2>::defaultVerbosity,
				const bool isUpdateWayR2L = Params<FPP,DEVICE,FPP2>::defaultUpdateWayR2L,
				const bool isFactSideLeft = Params<FPP,DEVICE,FPP2>::defaultFactSideLeft,
				const FPP init_lambda = Params<FPP,DEVICE,FPP2>::defaultLambda,
				const bool constant_step_size = Params<FPP,DEVICE,FPP2>::defaultConstantStepSize,
				const FPP step_size = Params<FPP,DEVICE,FPP2>::defaultStepSize): Params<FPP, DEVICE, FPP2>(nbRow, nbCol, nbFact, cons, init_fact, stop_crit_2facts, stop_crit_global, isVerbose, isUpdateWayR2L, isFactSideLeft, init_lambda, constant_step_size, step_size), init_D(init_D)
		{

		}

		ParamsFFT(
				const faust_unsigned_int nbRow,
				const faust_unsigned_int nbCol,
				const unsigned int nbFact,
				const std::vector<std::vector<const Faust::ConstraintGeneric*>> & cons,
				const std::vector<Faust::MatDense<FPP,DEVICE> >& init_fact,
				const Faust::Vect<FPP, DEVICE>& init_D_diag,
				const Faust::StoppingCriterion<FPP2>& stop_crit_2facts = StoppingCriterion<FPP2>(Params<FPP,DEVICE,FPP2>::defaultNiter1),
				const Faust::StoppingCriterion<FPP2>& stop_crit_global = StoppingCriterion<FPP2>(Params<FPP,DEVICE,FPP2>::defaultNiter2),
				const bool isVerbose = Params<FPP,DEVICE,FPP2>::defaultVerbosity,
				const bool isUpdateWayR2L = Params<FPP,DEVICE,FPP2>::defaultUpdateWayR2L,
				const bool isFactSideLeft = Params<FPP,DEVICE,FPP2>::defaultFactSideLeft,
				const FPP init_lambda = Params<FPP,DEVICE,FPP2>::defaultLambda,
				const bool constant_step_size = Params<FPP,DEVICE,FPP2>::defaultConstantStepSize,
				const FPP step_size = Params<FPP,DEVICE,FPP2>::defaultStepSize): Params<FPP, DEVICE, FPP2>(nbRow, nbCol, nbFact, cons, init_fact, stop_crit_2facts, stop_crit_global, isVerbose, isUpdateWayR2L, isFactSideLeft, init_lambda, constant_step_size, step_size), init_D(nbRow, nbCol)
		{
			init_D.setZeros();
			// set init_D from diagonal vector init_D_diag
			for(int i=0;i<nbRow;i++)
				init_D.getData()[i*nbRow+i] = init_D_diag.getData()[i];

		}

		ParamsFFT(
				const faust_unsigned_int nbRow,
				const faust_unsigned_int nbCol,
				const unsigned int nbFact,
				const std::vector<std::vector<const Faust::ConstraintGeneric*>> & cons,
				const std::vector<Faust::MatDense<FPP,DEVICE> >& init_fact,
				const FPP* init_D_diag,
				const Faust::StoppingCriterion<FPP2>& stop_crit_2facts = StoppingCriterion<FPP2>(Params<FPP,DEVICE,FPP2>::defaultNiter1),
				const Faust::StoppingCriterion<FPP2>& stop_crit_global = StoppingCriterion<FPP2>(Params<FPP,DEVICE,FPP2>::defaultNiter2),
				const bool isVerbose = Params<FPP,DEVICE,FPP2>::defaultVerbosity,
				const bool isUpdateWayR2L = Params<FPP,DEVICE,FPP2>::defaultUpdateWayR2L,
				const bool isFactSideLeft = Params<FPP,DEVICE,FPP2>::defaultFactSideLeft,
				const FPP init_lambda = Params<FPP,DEVICE,FPP2>::defaultLambda,
				const bool constant_step_size = Params<FPP,DEVICE,FPP2>::defaultConstantStepSize,
				const FPP step_size = Params<FPP,DEVICE,FPP2>::defaultStepSize): ParamsFFT<FPP, DEVICE, FPP2>(nbRow, nbCol, nbFact, cons, init_fact, Faust::Vect<FPP,DEVICE>(nbRow, init_D_diag), stop_crit_2facts, stop_crit_global, isVerbose, isUpdateWayR2L, isFactSideLeft, init_lambda, constant_step_size, step_size)
		{

		}



		ParamsFFT() {}

	};
}

#include "faust_ParamsFFT.hpp"

#endif
