template<typename FPP, FDevice DEVICE>
Faust::TransformHelper<FPP,DEVICE>* Faust::hierarchical(const Faust::MatDense<FPP,DEVICE>&  A,
		Params<FPP,DEVICE, Real<FPP>> & p,
		Real<FPP>& lambda, const bool compute_2norm_on_array,
		const bool on_gpu)
{
	auto S = new Faust::TransformHelper<FPP,DEVICE>(); // A is copied
	S->push_back(&A);
	if(on_gpu) S->enable_gpu_meth_for_mul();
	Faust::MatGeneric<FPP,DEVICE> *Si;
	const Faust::ConstraintGeneric *fac_cons, *res_cons;
	int zero_dims[2], id_dims[2];
	Faust::MatGeneric<FPP,DEVICE> *zero_mat = nullptr, *id_mat = nullptr;
	Faust::MatDense<FPP,DEVICE> * tmp_dense = nullptr;
	Faust::MatSparse<FPP,DEVICE> * tmp_sparse = nullptr;
	std::vector<Faust::MatGeneric<FPP,DEVICE>*> Si_vec;
	std::vector<Faust::ConstraintGeneric*> Si_cons;
	Real<FPP> lambda_ = p.init_lambda;
	Real<FPP> glo_lambda = 1;
	//TODO: remove these local variables and use directly p.
	const bool is_update_way_R2L = p.isUpdateWayR2L;
	const bool is_fact_side_left = p.isFactSideLeft;
	const bool use_csr = p.use_csr;
	const bool packing_RL = p.packing_RL;
	const Real<FPP> norm2_threshold = p.norm2_threshold;
	const unsigned int norm2_max_iter = p.norm2_max_iter;
	const double step_size = p.step_size;
	const bool constant_step_size = p.isConstantStepSize;
	std::vector<const Faust::ConstraintGeneric*> & fac_constraints = is_fact_side_left? p.cons[1]: p.cons[0];
	std::vector<const Faust::ConstraintGeneric*> & res_constraints = is_fact_side_left? p.cons[0]: p.cons[1];;


	if(p.isVerbose) p.Display();
	for(int i=0;i < fac_constraints.size();i++)
	{
		cout << "Faust::hierarchical: " << i+1 << endl;
		if(is_fact_side_left)
		{
			Si = S->get_gen_fact_nonconst(0);

		}
		else
			Si = S->get_gen_fact_nonconst(i);
		fac_cons = fac_constraints[i];
		res_cons = res_constraints[i];
		// init factors for the local optimization (factorization of Si)
		//TODO: refactor into a separate function init_zero_id
		if(! is_update_way_R2L && ! is_fact_side_left || is_update_way_R2L && is_fact_side_left)
		{
			zero_dims[0] = fac_cons->get_rows();
			zero_dims[1] = fac_cons->get_cols();
			id_dims[0] = res_cons->get_rows();
			id_dims[1] = res_cons->get_cols();
		}
		else //if(! is_update_way_R2L && is_fact_side_left || is_update_way_R2L && ! is_fact_side_left)
		{
			zero_dims[0] = res_cons->get_rows();
			zero_dims[1] = res_cons->get_cols();
			id_dims[0] = fac_cons->get_rows();
			id_dims[1] = fac_cons->get_cols();
		}
		if(use_csr)
		{
			tmp_sparse = new Faust::MatSparse<FPP,DEVICE>((faust_unsigned_int)zero_dims[0], (faust_unsigned_int)zero_dims[1]);
			tmp_sparse->setZeros();
			zero_mat = tmp_sparse;
			tmp_sparse = new Faust::MatSparse<FPP,DEVICE>((faust_unsigned_int)id_dims[0], (faust_unsigned_int)id_dims[1]);
			tmp_sparse->setEyes();
			id_mat = tmp_sparse;
		}
		else
		{

			tmp_dense = new Faust::MatDense<FPP,DEVICE>((faust_unsigned_int)zero_dims[0], (faust_unsigned_int)zero_dims[1]);
			tmp_dense->setZeros();
			zero_mat = tmp_dense;
			tmp_dense = new Faust::MatDense<FPP,DEVICE>((faust_unsigned_int)id_dims[0], (faust_unsigned_int)id_dims[1]);
			tmp_dense->setEyes();
			id_mat = tmp_dense;
		}

		if(is_update_way_R2L)
			Si_vec = {id_mat, zero_mat};
		else
			Si_vec = {zero_mat, id_mat};

		if(is_fact_side_left)
			Si_cons = { const_cast<Faust::ConstraintGeneric*>(res_cons), const_cast<Faust::ConstraintGeneric*>(fac_cons) };
		else
			Si_cons = { const_cast<Faust::ConstraintGeneric*>(fac_cons), const_cast<Faust::ConstraintGeneric*>(res_cons) };


		Faust::TransformHelper<FPP,DEVICE> Si_th(Si_vec, 1.0, false, false, true);
		//		Si_th.display();
		//		if(on_gpu) Si_th.enable_gpu_meth_for_mul();
		lambda_ = 1;
		tmp_dense = dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(Si);
		if(tmp_dense == nullptr)
		{
			tmp_sparse = dynamic_cast<Faust::MatSparse<FPP,DEVICE>*>(Si);
			tmp_dense = new MatDense<FPP,DEVICE>(*tmp_sparse);
		}
		else tmp_sparse = nullptr;
		Faust::palm4msa2(*tmp_dense, Si_cons, Si_th, lambda_, p.stop_crit_2facts, is_update_way_R2L , use_csr, packing_RL, compute_2norm_on_array,
				norm2_threshold, norm2_max_iter, constant_step_size, step_size, on_gpu);
		if(tmp_sparse != nullptr)
			// the Si factor has been converted into a MatDense in the memory
			// storage
			// delete it // TODO: palm4msa2 should handle this on its own
			delete tmp_dense;

		//prepare global optimization
		glo_lambda *= lambda_;
		if(is_fact_side_left)
		{
			S->pop_front();
			S->push_first(*(Si_th.begin()+1), false, false);
			S->push_first(*(Si_th.begin()), false, false);
		}
		else
		{
			S->pop_back();
			S->push_back(*(Si_th.begin()), false, false);
			S->push_back(*(Si_th.begin()+1), false, false);
		}
		//TODO: verify if the constraints order doesn't depend on
		//is_fact_side_left
		std::vector<Faust::ConstraintGeneric*> glo_cons;
		for(auto ite_cons=fac_constraints.begin(); ite_cons != fac_constraints.begin()+i+1;ite_cons++)
			if(is_fact_side_left)
				glo_cons.insert(glo_cons.begin(), const_cast<Faust::ConstraintGeneric*>(*ite_cons));
			else
				glo_cons.push_back(const_cast<Faust::ConstraintGeneric*>(*ite_cons));
			if(is_fact_side_left)
				glo_cons.insert(glo_cons.begin(), const_cast<Faust::ConstraintGeneric*>(res_constraints[i]));
			else
				glo_cons.push_back(const_cast<Faust::ConstraintGeneric*>(res_constraints[i]));

		// global optimization
		Faust::palm4msa2(A, glo_cons, *S, glo_lambda, p.stop_crit_global ,is_update_way_R2L, use_csr, packing_RL, compute_2norm_on_array,
				norm2_threshold, norm2_max_iter, constant_step_size, step_size, on_gpu);
	}
	lambda = glo_lambda;
	return S;
}

template<typename FPP, FDevice DEVICE>
Faust::TransformHelper<FPP,DEVICE>* Faust::hierarchical(const Faust::MatDense<FPP,DEVICE>& A,
//        const int nites,
		std::vector<StoppingCriterion<Real<FPP>>>& sc,
        std::vector<const Faust::ConstraintGeneric*> & fac_constraints,
        std::vector<const Faust::ConstraintGeneric*> & res_constraints,
        Real<FPP>& lambda,
        const bool is_update_way_R2L, const bool is_fact_side_left,
        const bool use_csr, const bool packing_RL,
        const bool compute_2norm_on_array,
        const Real<FPP> norm2_threshold,
        const unsigned int norm2_max_iter, const bool is_verbose,
		const bool constant_step_size, const Real<FPP> step_size,
		const bool on_gpu)
{
    Faust::Params<FPP,DEVICE,Real<FPP>> p(A.getNbRow(), A.getNbCol(), fac_constraints.size()+1, {fac_constraints, res_constraints}, {}, sc[0], sc[1], is_verbose, is_update_way_R2L, is_fact_side_left, lambda, constant_step_size, step_size);
	p.use_csr = use_csr;
	p.packing_RL = packing_RL;
	p.norm2_threshold = norm2_threshold;
	p.norm2_max_iter = norm2_max_iter;
    return Faust::hierarchical(A, p, lambda, compute_2norm_on_array, on_gpu);
}
