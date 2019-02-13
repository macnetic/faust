#ifndef __FAUST_PARAMS_PALM_FFT_H__
#define __FAUST_PARAMS_PALM_FFT_H__

#include "faust_ParamsPalm.h"

using namespace Faust;

namespace Faust
{


	template<typename FPP, Device DEVICE, typename FPP2 = double>
		class ParamsPalmFFT : public Faust::ParamsPalm<FPP,DEVICE,FPP2>
	{

		public:

			//ctor definitions in header because it consists mainly to call parent ctor
			ParamsPalmFFT(const Faust::MatDense<FPP,DEVICE>& data_,
					const int nbFact_,
					const std::vector<const Faust::ConstraintGeneric*>& cons_,
					const std::vector<Faust::MatDense<FPP,DEVICE> >& init_fact_,
					const Faust::MatDense<FPP,DEVICE>& init_D,
					const Faust::StoppingCriterion<FPP2> & stop_crit_ = StoppingCriterion<FPP2>(ParamsPalm<FPP,DEVICE,FPP2>::defaultNiter),
					const bool isVerbose_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultVerbosity ,
					const bool isUpdateWayR2L_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultUpdateWayR2L ,
					const FPP init_lambda_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultLambda,
					const FPP step_size_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultStepSize) : ParamsPalm<FPP, DEVICE, FPP2>(data_, nbFact_, cons_, init_fact_, stop_crit_, isVerbose_, isUpdateWayR2L_, init_lambda_, true /*constant_step_size is always true for Palm4MSAFFT */, step_size_), init_D(init_D) {}

			ParamsPalmFFT() : ParamsPalm<FPP,DEVICE,FPP2>(), init_D(0,0) {}

			ParamsPalmFFT(const Faust::MatDense<FPP,DEVICE>& data_,
					const int nbFact_,
					const std::vector<const Faust::ConstraintGeneric*>& cons_,
					const std::vector<Faust::MatDense<FPP,DEVICE> >& init_fact_,
					const Faust::Vect<FPP,DEVICE>& init_D_diag,
					const Faust::StoppingCriterion<FPP2> & stop_crit_ = StoppingCriterion<FPP2>(ParamsPalm<FPP,DEVICE,FPP2>::defaultNiter),
					const bool isVerbose_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultVerbosity ,
					const bool isUpdateWayR2L_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultUpdateWayR2L ,
					const FPP init_lambda_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultLambda,
					const FPP step_size_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultStepSize);

			ParamsPalmFFT(const Faust::MatDense<FPP,DEVICE>& data_,
					const int nbFact_,
					const std::vector<const Faust::ConstraintGeneric*>& cons_,
					const std::vector<Faust::MatDense<FPP,DEVICE> >& init_fact_,
					const FPP* init_D_diag,
					const Faust::StoppingCriterion<FPP2> & stop_crit_ = StoppingCriterion<FPP2>(ParamsPalm<FPP,DEVICE,FPP2>::defaultNiter),
					const bool isVerbose_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultVerbosity ,
					const bool isUpdateWayR2L_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultUpdateWayR2L ,
					const FPP init_lambda_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultLambda,
					const FPP step_size_ = ParamsPalm<FPP,DEVICE,FPP2>::defaultStepSize) : ParamsPalmFFT<FPP,DEVICE,FPP2>(data_, nbFact_, cons_, init_fact_, Faust::Vect<FPP,DEVICE>(data_.getNbRow(), init_D_diag), stop_crit_, isVerbose_, isUpdateWayR2L_, init_lambda_, step_size_) {}


			MatDense<FPP,DEVICE> init_D;


	};

#include "faust_ParamsPalmFFT.hpp"

}
#endif
