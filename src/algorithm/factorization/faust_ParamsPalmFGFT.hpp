template<typename FPP,FDevice DEVICE,typename FPP2>
 ParamsPalmFGFT<FPP,DEVICE,FPP2>::ParamsPalmFGFT(const Faust::MatDense<FPP,DEVICE>& data_,
					const int nbFact_,
					const std::vector<const Faust::ConstraintGeneric*>& cons_,
					const std::vector<Faust::MatDense<FPP,DEVICE> >& init_fact_,
					const Faust::Vect<FPP,DEVICE>& init_D_diag,
					const Faust::StoppingCriterion<FPP2> & stop_crit_ ,
					const bool isVerbose_ ,
					const bool isUpdateWayR2L_ ,
					const FPP2 init_lambda_ ,
					const FPP2 step_size_, const Faust::GradientCalcOptMode gradCalcOptMode) : ParamsPalm<FPP, DEVICE, FPP2>(data_, nbFact_, cons_, init_fact_, stop_crit_, isVerbose_, isUpdateWayR2L_, init_lambda_, true /*constant_step_size is always true for Palm4MSAFGFT */, step_size_),  init_D(data_.getNbRow(), data_.getNbCol())
{
	init_D.setZeros();
	// set init_D from diagonal vector init_D_diag
	for(int i=0;i<data_.getNbRow();i++)
		init_D.getData()[i*init_D.getNbRow()+i] = init_D_diag.getData()[i];
}
