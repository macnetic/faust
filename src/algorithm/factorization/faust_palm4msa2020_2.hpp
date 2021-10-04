template <typename FPP, FDevice DEVICE>
void Faust::palm4msa2(const Faust::MatDense<FPP,DEVICE>& A,
		std::vector<Faust::ConstraintGeneric*> & constraints,
		Faust::TransformHelper<FPP,DEVICE>& S,
		Real<FPP>& lambda, //TODO: FPP lambda ? is it useful to have a complex lamdba ?
		//const unsigned int nites,
		const StoppingCriterion<Real<FPP>>& sc,
		const bool is_update_way_R2L,
		const FactorsFormat factors_format,
		const bool packing_RL,
		const MHTPParams<Real<FPP>> mhtp_params/*=MHTPParams<FPP>()*/,
		const bool compute_2norm_on_array,
		const Real<FPP> norm2_threshold,
		const unsigned int norm2_max_iter,
		bool constant_step_size, Real<FPP> step_size,
		const bool on_gpu /*=false*/,
		const bool is_verbose/*=false*/, const int id/*=0*/)
{
	std::chrono::duration<double> norm2_duration = std::chrono::duration<double>::zero();
	std::chrono::duration<double> fgrad_duration = std::chrono::duration<double>::zero();
	double norm1, norm2;
	/* variable environment parameters (which are interesting enough for debugging/profiling but not yet for the wrapper user API */
	char* str_env_prod_mod = getenv("PROD_MOD");
	int prod_mod = DYNPROG; // GREEDY_ALL_BEST_GENMAT; DYNPROG is a bit better to factorize the MEG matrix and not slower to factorize a Hadamard matrix
	if(str_env_prod_mod)
		prod_mod = std::atoi(str_env_prod_mod);
	auto str_env_no_lambda_error = getenv("NO_LAMBDA_ERROR");
	bool no_lambda_error = false;
	if(str_env_no_lambda_error)
		no_lambda_error = (bool) std::atoi(str_env_no_lambda_error);
	bool use_grad1 = false;
	auto str_env_use_grad1 = getenv("USE_GRAD1");
	if(str_env_use_grad1)
		use_grad1 = (bool) std::atoi(str_env_use_grad1);
	/******************************************************/
//	std::cout << "palm4msa2 "<< std::endl;
	if(constraints.size() == 0)
		throw out_of_range("No constraint passed to palm4msa.");
	const Real<FPP> lipschitz_multiplicator = 1.001;
	Faust::MatGeneric<FPP,DEVICE>* cur_fac;
	Faust::MatSparse<FPP,DEVICE>* scur_fac;
	Faust::MatDense<FPP,DEVICE>* dcur_fac;
	unsigned int nfacts = constraints.size();
	Real<FPP> error = -1; //negative error is ignored
	std::vector<std::pair<faust_unsigned_int,faust_unsigned_int>> dims;
	int norm2_flag; // return val
	for(auto c: constraints)
		dims.push_back(make_pair(c->get_rows(), c->get_cols()));
	//TODO: make it possible to receive a MatSparse A
	if(is_verbose && mhtp_params.used)
	{
		std::cout << mhtp_params.constant_step_size << std::endl;
		std::cout << mhtp_params.to_string() << std::endl;
	}
	Faust::MatDense<FPP,DEVICE> A_H = A;
	A_H.adjoint();
	if(S.size() != nfacts)
		fill_of_eyes(S, nfacts, factors_format != AllDense, dims, on_gpu);
	else if(factors_format == AllSparse)
	{
		S.convertToSparse();
	}
	else if(factors_format == AllDense)
	{
		S.convertToDense();
	}
	int i = 0, f_id, j;
	std::function<void()> init_ite, next_fid;
	std::function<bool()> updating_facs;
	std::function<bool()> is_last_fac_updated;
	// packed Fausts corresponding to each factor
	std::vector<TransformHelper<FPP,DEVICE>*> pL, pR;
	pL.resize(nfacts);// pL[i] is the Faust for all factors to the left of the factor *(S.begin()+i)
	pR.resize(nfacts);// pR[i] the same for the right of S_i
	Faust::MatDense<FPP,DEVICE> D, tmp;
	Faust::MatSparse<FPP,DEVICE> spD;
	Real<FPP> c = 1/step_size;

	for(int i=0;i<nfacts;i++)
	{
		pL[i] = new TransformHelper<FPP,DEVICE>();
		pR[i] = new TransformHelper<FPP,DEVICE>();
	}
	if(is_update_way_R2L)
	{
		init_ite = [&f_id, &nfacts, &pL, &S, &packing_RL, &on_gpu, &prod_mod]()
		{
			//pre-compute left products for each Si
			if(pL[0] != nullptr) delete pL[0];
			pL[0] = new TransformHelper<FPP,DEVICE>(); // empty faust // no factors to the left of *(S.begin())
			for(int i=1;i < nfacts; i++)
			{
				auto vec_Si_minus_1 = { *(S.begin()+i-1) };
				if(pL[i] != nullptr) delete pL[i]; //TODO: maybe to replace by a TransformHelper stored in the stack to avoid deleting each time
				pL[i] = new TransformHelper<FPP,DEVICE>(*pL[i-1], vec_Si_minus_1);
				if(packing_RL) ((TransformHelperGen<FPP,DEVICE>*)pL[i])->pack_factors(prod_mod);
			}
			// all pL[i] Fausts are composed at most of one factor matrix
			f_id = nfacts-1;
		};
		next_fid = [&f_id, &pR, &S, &packing_RL, &on_gpu, &prod_mod]()
		{
			if(f_id > 0)
			{
				if(pR[f_id-1] != nullptr)
					delete pR[f_id-1];
				auto vec_Sj = { *(S.begin()+f_id) };
				pR[f_id-1] = new Faust::TransformHelper<FPP,DEVICE>(vec_Sj, *pR[f_id]);
				if(packing_RL) ((TransformHelperGen<FPP,DEVICE>*)pR[f_id-1])->pack_factors(prod_mod);
			}
			f_id--;
		};
		is_last_fac_updated = [&f_id]() {return f_id == 0;};
		updating_facs = [&f_id]() {return f_id >= 0;};
	}
	else
	{
		init_ite = [&f_id, &pR, &S, &packing_RL, &on_gpu, &prod_mod]()
		{
			if(pR[S.size()-1] != nullptr) delete pR[S.size()-1];
			pR[S.size()-1] = new TransformHelper<FPP,DEVICE>(); // empty faust // no factors to the right of *(S.begin()+S.size()]-1)
			for(int i=S.size()-2;i >= 0; i--)
			{
				auto vec_Si_plus_1 = { *(S.begin()+i+1) };
				if(pR[i] != nullptr) delete pR[i];
				pR[i] = new TransformHelper<FPP,DEVICE>(vec_Si_plus_1, *pR[i+1]);
				if(packing_RL) ((TransformHelperGen<FPP,DEVICE>*)pR[i])->pack_factors(prod_mod);
			}
			f_id = 0;
		};
		next_fid = [&f_id, &S, &pL, &nfacts, packing_RL, &on_gpu, &prod_mod]()
		{
			if(f_id < nfacts-1)
			{
				if(pL[f_id+1] != nullptr)
					delete pL[f_id+1];
				auto vec_Sj = { *(S.begin()+f_id) };
				pL[f_id+1] = new Faust::TransformHelper<FPP,DEVICE>(*pL[f_id], vec_Sj);
				if(packing_RL) ((TransformHelperGen<FPP,DEVICE>*)pL[f_id+1])->pack_factors(prod_mod);
			}
			f_id++;
		};
		updating_facs = [&f_id, &nfacts]() {return f_id < nfacts;};
		is_last_fac_updated = [&f_id, &nfacts]() {return f_id == nfacts-1;};
	}

	while(sc.do_continue(i, error))
	{
//		std::cout << "i: " <<  i << std::endl;
//		std::cout << "nfacts:" << nfacts << std::endl;

		init_ite();
		while(updating_facs())
		{
			//						std::cout << "#f_id: " << f_id << std::endl;
			cur_fac = S.get_gen_fact_nonconst(f_id);
			if(mhtp_params.used && i%mhtp_params.palm4msa_period == 0)
			{
				perform_MHTP(mhtp_params, A, A_H, S, f_id, pL, pR, packing_RL,
						is_verbose, *constraints[f_id], norm2_max_iter, norm2_threshold,
						norm2_duration,
						fgrad_duration,
						sc, error, factors_format, prod_mod, c, lambda);
			}
			else
				update_fact(cur_fac, f_id, A, S, pL, pR, packing_RL,
						is_verbose, *constraints[f_id], norm2_max_iter, norm2_threshold,
						norm2_duration,
						fgrad_duration,
						constant_step_size, step_size,
						sc, error, factors_format, prod_mod, c, lambda, use_grad1);

			next_fid(); // f_id updated to iteration factor index (pL or pR too)
		}
		//update lambda
		update_lambda(S, pL, pR, A_H, lambda, no_lambda_error);
		if(is_verbose)
		{
			set_calc_err_ite_period(); //macro setting the variable ite_period
			if((! (i%ite_period)) && ite_period > 0)
			{
				std::cout << "PALM4MSA2020 iteration: " << i;
				auto err = calc_rel_err(S, A, lambda);
				std::cout << " relative error: " << err;
				std::cout << " (call id: " << id << ")" << std::endl;
				std::cout << " lambda=" << lambda << std::endl;
			}
		}
		i++;
	}
	S.update_total_nnz();
	// free the latest pR and pL TransformHelpers
	for(int i=0;i<nfacts;i++)
	{
		if(pL[i] != nullptr)
			delete pL[i];
		if(pR[i] != nullptr)
			delete pR[i];
	}
	if(is_verbose)
	{
		std::cout << "palm4msa spectral time=" << norm2_duration.count() << std::endl;
		std::cout << "palm4msa fgrad time=" << fgrad_duration.count() << std::endl;
	}
}

template <typename FPP, FDevice DEVICE>
void Faust::compute_n_apply_grad1(const int f_id, const Faust::MatDense<FPP,DEVICE> &A, Faust::TransformHelper<FPP,DEVICE>& S, std::vector<TransformHelper<FPP,DEVICE>*> &pL, std::vector<TransformHelper<FPP,DEVICE>*>& pR, const bool packing_RL, const Real<FPP>& lambda, const Real<FPP> &c, Faust::MatDense<FPP,DEVICE> &out /* D */, const StoppingCriterion<Real<FPP>>& sc, Real<FPP> &error, const int prod_mod)
{
	Faust::MatDense<FPP,DEVICE> tmp;
	Faust::MatDense<FPP,DEVICE> & D = out;
//	Faust::MatGeneric<FPP,DEVICE> * LorR;
	Faust::MatDense<FPP,DEVICE> * LorR;
	Faust::MatDense<FPP,DEVICE> _LorR;
	auto S_j_vec = {*(S.begin()+f_id)};
	Faust::TransformHelper<FPP, DEVICE> _LSR(/*lambda_vec,*/ *pL[f_id], S_j_vec, *pR[f_id]);
	//			tmp = _LSR.get_product(prod_mod);
	_LSR.get_product(tmp, prod_mod);
//	_LSR.get_product(tmp);
	tmp *= FPP(lambda);
	tmp -= A;
	if(sc.isCriterionErr())
		error = tmp.norm();
	FPP alpha_R = 1, alpha_L = 1, beta_R = 0, beta_L = 0; //decl in parent scope
	auto pR_sz = pR[f_id]->size();
	auto pL_sz = pL[f_id]->size();
	if(pR_sz > 0)
	{
		if(pR_sz == 1 && packing_RL) // packing_RL == true
			LorR = dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(pR[f_id]->get_gen_fact_nonconst(0));
//			LorR = pR[f_id]->get_gen_fact_nonconst(0); //normally pR[f_id] is packed (hence reduced to a single MatDense)
		else
		{
			_LorR = pR[f_id]->get_product(prod_mod);
			LorR = &_LorR;
		}
		if(pL_sz == 0)
		{ //no L factor for factor f_id
			alpha_R = - lambda/c;
			beta_R = 1;
			gemm(tmp, *LorR, D, alpha_R, beta_R, 'N', 'H');
//			gemm_gen(tmp, *LorR, D, alpha_R, beta_R, 'N', 'H');
		}
		else
			gemm(tmp, *LorR, tmp, alpha_R, beta_R, 'N', 'H');
//			gemm_gen(tmp, *LorR, tmp, alpha_R, beta_R, 'N', 'H');
	}
	if(pL_sz > 0)
	{
		if(pL_sz == 1 && packing_RL) // packing_RL == true
			LorR = dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(pL[f_id]->get_gen_fact_nonconst(0));
//			LorR = pL[f_id]->get_gen_fact_nonconst(0);
		else
		{
			_LorR = pL[f_id]->get_product(prod_mod);
			LorR = &_LorR;
		}
		alpha_L = -lambda/c;
		beta_L = 1;
		gemm(*LorR, tmp, D, alpha_L, beta_L, 'H', 'N');
//		gemm_gen(*LorR, tmp, D, alpha_L, beta_L, 'H', 'N');
	}
}

template <typename FPP, FDevice DEVICE>
void Faust::compute_n_apply_grad2(const int f_id, const Faust::MatDense<FPP,DEVICE> &A, Faust::TransformHelper<FPP,DEVICE>& S, std::vector<TransformHelper<FPP,DEVICE>*> &pL, std::vector<TransformHelper<FPP,DEVICE>*>& pR, const bool packing_RL, const Real<FPP>& lambda, const Real<FPP> &c, Faust::MatDense<FPP,DEVICE> &out /* D */, const StoppingCriterion<Real<FPP>>& sc, Real<FPP> &error, const int prod_mod)
{
	Faust::MatDense<FPP,DEVICE> tmp;
	Faust::MatDense<FPP,DEVICE> grad_over_c;
	Faust::MatDense<FPP,DEVICE> & D = out;
	Faust::MatGeneric<FPP,DEVICE> *_L, *_R;
	Faust::MatDense<FPP,DEVICE> __L, __R;
	Faust::MatDense<FPP,DEVICE> * LorR;
	Faust::MatDense<FPP,DEVICE> _LorR;
	std::vector<Faust::MatGeneric<FPP,DEVICE>*> facts;
	std::vector<char> tc_flags;
//#define mul_3_facts multiply_order_opt
#define mul_3_facts multiply_order_opt_all_ends// this one only optimizes the product on factor ends but for three factors it doesn't change anything comparing to multiply_order_opt
	tmp = A;
	auto pR_sz = pR[f_id]->size();
	auto pL_sz = pL[f_id]->size();
	if(pR_sz > 0)
	{
		if(pR_sz == 1)
			_R = pR[f_id]->get_gen_fact_nonconst(0);
		else
		{
//			__R = pR[f_id]->get_product(prod_mod); // disabled because GREEDY_ALL_BEST_GENMAT is slower than DEFAULT_L2R for Hadamard factorization
			__R = pR[f_id]->get_product();
			_R = &__R;
		}

	}
	if(pL_sz > 0)
	{
		if(pL_sz == 1)
			_L = pL[f_id]->get_gen_fact_nonconst(0);
		else
		{
//			__L = pL[f_id]->get_product(prod_mod); // disabled because GREEDY_ALL_BEST_GENMAT is slower than DEFAULT_L2R for Hadamard factorization
			__L = pL[f_id]->get_product();
			_L = &__L;
		}
	}
	if(pR_sz > 0 && pL_sz > 0)
	{
		// compute error = m_lambda*L*S*R-data
		facts = { _L, &D, _R };
		mul_3_facts(facts, tmp, (FPP) lambda, (FPP) -1.0);
		if(sc.isCriterionErr())
			error = tmp.norm();
		// compute m_lambda/c * L'*error*R'
		facts = { _L, &tmp, _R };
		tc_flags = {'H', 'N', 'H'};
	}
	else if(pR_sz > 0)
	{
		// compute error = m_lambda*L*S*R-data
		facts = { &D, _R };
		mul_3_facts(facts, tmp, (FPP) lambda, (FPP) -1.0);
		if(sc.isCriterionErr())
			error = tmp.norm();
		// compute m_lambda/c * L'*error*R'
		facts = { &tmp, _R };
		tc_flags = { 'N', 'H'};
	}
	else //if(pL_sz > 0)
	{
		// compute error = m_lambda*L*S*R-data
		facts = { _L, &D};
		mul_3_facts(facts, tmp, (FPP) lambda, (FPP) -1.0);
		if(sc.isCriterionErr())
			error = tmp.norm();
		// compute m_lambda/c * L'*error*R'
		facts = { _L, &tmp};
		tc_flags = {'H', 'N'};
	}
	// mul_3_facts(facts, D, (FPP) - lambda/c, (FPP)1, tc_flags); // this one has showed error in calculation when using complex matrices (DFT factorization)
	// do it two steps: 1) compute the gradient, then 2) apply it separately to D (the current factor)
	mul_3_facts(facts, grad_over_c, (FPP) lambda/c, (FPP)0, tc_flags);
	D -= grad_over_c;
}

	template<typename FPP, FDevice DEVICE>
void Faust::update_lambda(Faust::TransformHelper<FPP,DEVICE>& S, std::vector<TransformHelper<FPP,DEVICE>*> &pL, std::vector<TransformHelper<FPP,DEVICE>*>& pR, const MatDense<FPP, DEVICE> &A_H, Real<FPP>& lambda, bool no_lambda_error/*= false*/)
{
	Faust::MatDense<FPP,DEVICE> A_H_S;
	MatDense<FPP, DEVICE> S_mat;
	FPP tr; // A_H_S trace
	Real<FPP> nS; // S fro. norm
	auto n = S.size();
	bool packing_RL = (pR[0] == nullptr || pR[0]->size() == 1) && (pL[n-1] == nullptr || pL[n-1]->size() == 1);
	// compute S full matrix
	if(packing_RL)
	{
		if(pR[0] == nullptr || pL[n-1] == nullptr)
			throw std::logic_error("update_lambda: pR and pL weren't properly initialized.");
		// optimize S product by re-using R[0] or L[0] (S == S[0]*R[0] or S == L[n-1]*S[n-1])
		if(S.get_gen_fact(0)->getNbCol() * pR[0]->getNbRow() < pL[n-1]->getNbCol() * S.get_gen_fact(n-1)->getNbRow())
		{
			auto S0 = { S.get_gen_fact(0) };
			TransformHelper<FPP, DEVICE> _S(S0, *pR[0]);
			_S.get_product(S_mat);
		}
		else
		{
			auto Sn = { S.get_gen_fact(n-1) };
			TransformHelper<FPP, DEVICE> _S(*pL[n-1], Sn);
			_S.get_product(S_mat);
		}
	}
	else
		S.get_product(S_mat);
	gemm(A_H, S_mat, A_H_S, (FPP) 1.0, (FPP) 0.0, 'N', 'N');
	tr = A_H_S.trace();
	nS = S_mat.norm();
	if(Real<FPP>(0) == nS)
		if(no_lambda_error)
		{
			// don't change lambda
			std::cout << "WARNING: lambda didn't change because S Fro. norm is zero." << std::endl;
			return;
		}
		else
			throw std::runtime_error("Error in update_lambda: S Frobenius norm is zero, can't compute lambda.");
	if(std::isnan(std::real(tr)) || std::isnan(nS))
		if(no_lambda_error)
		{
			// don't change lambda
			std::cout << "WARNING: lambda didn't change because S contains NaN." << std::endl;
			return;
		}
		else
			throw std::runtime_error("Error in update_lambda: S (the Faust) contains nan elements in at least one of its matrices, can't compute lambda.");
	lambda = std::real(tr)/(nS*nS);
}

	template<typename FPP, FDevice DEVICE>
void Faust::update_fact(
		Faust::MatGeneric<FPP,DEVICE>* cur_fac,
		int f_id,
		const Faust::MatDense<FPP,DEVICE>& A,
		Faust::TransformHelper<FPP,DEVICE>& S,
		std::vector<TransformHelper<FPP,DEVICE>*> &pL,
		std::vector<TransformHelper<FPP,DEVICE>*> &pR,
		const bool packing_RL,
		const bool is_verbose,
		const Faust::ConstraintGeneric &constraint,
		const int norm2_max_iter,
		const Real<FPP>& norm2_threshold,
		std::chrono::duration<double>& norm2_duration,
		std::chrono::duration<double>& fgrad_duration,
		const bool constant_step_size,
		const Real<FPP> step_size,
		const StoppingCriterion<Real<FPP>>& sc,
		Real<FPP> &error,
		const FactorsFormat factors_format,
		const int prod_mod,
		Real<FPP> &c,
		const Real<FPP>& lambda,
		bool use_grad1/*= false*/)
{
	int norm2_flag;
	std::chrono::time_point<std::chrono::high_resolution_clock> spectral_stop, spectral_start;
	std::chrono::time_point<std::chrono::high_resolution_clock> fgrad_stop, fgrad_start;
	Faust::MatSparse<FPP,DEVICE>* scur_fac = nullptr;
	Faust::MatDense<FPP,DEVICE>* dcur_fac = nullptr;
	Faust::MatDense<FPP,DEVICE> D;
	Faust::MatSparse<FPP,DEVICE> spD;

	Real<FPP> nR=1,nL=1;
	if(constant_step_size)
		c = 1 / step_size;
	else
	{
		if(is_verbose)
			spectral_start = std::chrono::high_resolution_clock::now();
		if(pR[f_id]->size() > 0)
			nR = pR[f_id]->spectralNorm(norm2_max_iter, norm2_threshold, norm2_flag);
		if(pL[f_id]->size() > 0)
			nL = pL[f_id]->spectralNorm(norm2_max_iter, norm2_threshold, norm2_flag);
		if(is_verbose)
		{
			spectral_stop = std::chrono::high_resolution_clock::now();
			norm2_duration += spectral_stop-spectral_start;
		}
		c = LIPSCHITZ_MULTIPLICATOR*lambda*lambda*nR*nR*nL*nL;
	}
	if(S.is_fact_sparse(f_id))
	{
		scur_fac = dynamic_cast<Faust::MatSparse<FPP,DEVICE>*>(cur_fac);
		D = *scur_fac;
		dcur_fac = nullptr;
	}
	else
	{
		dcur_fac = dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(cur_fac); // TOFIX: possible it is not Dense... but MatDiag (sanity check on function start)
		D = *dcur_fac;
		scur_fac = nullptr;
	}

	if(is_verbose)
		fgrad_start = std::chrono::high_resolution_clock::now();
	if(typeid(D) != typeid(Faust::MatDense<FPP,Cpu>) || use_grad1)
		// compute_n_apply_grad2 is not yet supported with GPU2
		compute_n_apply_grad1(f_id, A, S, pL, pR, packing_RL, lambda, c, D, sc, error, prod_mod);
	else
		compute_n_apply_grad2(f_id, A, S, pL, pR, packing_RL, lambda, c, D, sc, error, prod_mod);
	if(is_verbose)
	{
		fgrad_stop = std::chrono::high_resolution_clock::now();
		fgrad_duration += fgrad_stop-fgrad_start;
	}

	// Update the S factor
	if(factors_format == AllDynamic)
	{
		// new_fac is a MatGeneric but underliying concrete object can be a MatSparse or a MatDense
		// new_fac is allocated in the heap
		// replace the former fact by the new one
		auto new_fac = constraint.project_gen<FPP,DEVICE,Real<FPP>>(D);
		S.replace(new_fac, f_id);
	}
	else // factors_format == AllDense or AllSparse
	{
		constraint.project<FPP,DEVICE,Real<FPP>>(D);
		// D is the prox image (always a MatDense
		// convert D to the proper format (MatSparse or MatDense)

		if(factors_format == AllSparse && dcur_fac != nullptr || factors_format == AllDense && scur_fac != nullptr)
			throw std::runtime_error("Current factor is inconsistent with the configured factors_format.");

		if(factors_format == AllSparse)
		{
			// convert to sparse then update
			spD = D;
			S.update(spD, f_id); // update is at higher level than a simple assignment
		}
		else
		{ // factors_format == AllDense
			// directly update (D is a MatDense)
			S.update(D, f_id);
		}
	}


}


