template <typename FPP, FDevice DEVICE>
void Faust::palm4msa2(const Faust::MatDense<FPP,DEVICE>& A,
		std::vector<Faust::ConstraintGeneric*> & constraints,
		Faust::TransformHelper<FPP,DEVICE>& S,
		Real<FPP>& lambda, //TODO: FPP lambda ? is it useful to have a complex lamdba ?
		//const unsigned int nites,
		const StoppingCriterion<Real<FPP>>& sc,
		const bool is_update_way_R2L,
		const bool use_csr,
		const bool packing_RL,
		const MHTPParams<FPP> mhtp_params/*=MHTPParams<FPP>()*/,
		const bool compute_2norm_on_array,
		const Real<FPP> norm2_threshold,
		const unsigned int norm2_max_iter,
		bool constant_step_size, Real<FPP> step_size,
		const bool on_gpu /*=false*/,
		const bool is_verbose/*=false*/, const int id/*=0*/)
{
	std::chrono::time_point<std::chrono::high_resolution_clock> spectral_stop, spectral_start;
	std::chrono::duration<double> spectral_duration = std::chrono::duration<double>::zero();
	std::chrono::time_point<std::chrono::high_resolution_clock> fgrad_stop, fgrad_start;
	std::chrono::duration<double> fgrad_duration = std::chrono::duration<double>::zero();
	int prod_mod = ORDER_ALL_BEST_MIXED;
	double norm1, norm2;
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
	if(is_verbose)
	{
		std::cout << "use_MHTP: " << mhtp_params.used << std::endl;
		std::cout<<"MHTP stop crit.: "<< std::endl << mhtp_params.sc.to_string() <<std::endl;
	}
	Faust::MatDense<FPP,DEVICE> A_H = A;
	A_H.adjoint();
	if(S.size() != nfacts)
		fill_of_eyes(S, nfacts, use_csr, dims, on_gpu);
	int i = 0, f_id, j;
	std::function<void()> init_ite, next_fid;
	std::function<bool()> updating_facs;
	std::function<bool()> is_last_fac_updated;
	std::function<void(Faust::MatGeneric<FPP,DEVICE>*, int f_id)> update_fac;
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
				// if the ctor args are GPU-enabled so is pL[i]
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
	update_fac = [&A, &D, &spD, &constant_step_size, &is_verbose, &spectral_start, &spectral_stop, &spectral_duration,
			   &fgrad_start, &fgrad_stop, &fgrad_duration,
			   &c, &lipschitz_multiplicator, &sc, &error, &prod_mod, &packing_RL, &constraints,
			   &pR, &pL, &norm2_max_iter, &norm2_threshold, &norm2_flag, &lambda, &S, &scur_fac, &dcur_fac,
			   &use_csr]
				   (Faust::MatGeneric<FPP,DEVICE>* cur_fac, int f_id)
	{
			Real<FPP> nR=1,nL=1;
			if(! constant_step_size)
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
					spectral_duration += spectral_stop-spectral_start;
				}
				c = lipschitz_multiplicator*lambda*lambda*nR*nR*nL*nL;
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
			if(typeid(D) != typeid(Faust::MatDense<FPP,Cpu>))
				// compute_n_apply_grad2 is not yet supported with GPU2
				compute_n_apply_grad1(f_id, A, S, pL, pR, lambda, c, D, sc, error, prod_mod, packing_RL);
			else
				compute_n_apply_grad2(f_id, A, S, pL, pR, lambda, c, D, sc, error, prod_mod, packing_RL);
			if(is_verbose)
			{
				fgrad_stop = std::chrono::high_resolution_clock::now();
				fgrad_duration += fgrad_stop-fgrad_start;
			}
			// really update now
			constraints[f_id]->project<FPP,DEVICE,Real<FPP>>(D);


			if(use_csr && dcur_fac != nullptr || !use_csr && scur_fac != nullptr)
				throw std::runtime_error("Current factor is inconsistent with use_csr.");

			if(use_csr)
			{
				spD = D;
				S.update(spD, f_id); // update is at higher level than a simple assignment
			}
			else
			{
				S.update(D, f_id);
			}

	};
	while(sc.do_continue(i, error))
	{
//		std::cout << "i: " <<  i << std::endl;
//		std::cout << "nfacts:" << nfacts << std::endl;

		init_ite();
		while(updating_facs())
		{
			//						std::cout << "#f_id: " << f_id << std::endl;
			cur_fac = S.get_gen_fact_nonconst(f_id);
			if(i%1000 == 0 && mhtp_params.used)
			{
				if(is_verbose)
					std::cout << "MHTP" << std::endl;
				j = 0;
				// set the factor to zero
				cur_fac->setZeros();
				while(mhtp_params.sc.do_continue(j)) // TODO: what about the error stop criterion?
				{
					if(mhtp_params.constant_step_size)
					{
						//TODO: add arguments to update_fac in order to avoid shunting PALM4MSA parameters with MHTP's
						constant_step_size = true;
						step_size = mhtp_params.step_size;
						c = 1 / step_size;
					}
					update_fac(cur_fac, f_id);
					j++;
					Faust::MatDense<FPP,DEVICE> A_H_S = S.multiply(A_H);
					Real<FPP> trr = std::real(A_H_S.trace());
					Real<FPP> n = S.normFro();
					if(std::numeric_limits<Real<FPP>>::epsilon() >= n)
						throw std::runtime_error("Faust Frobenius norm is zero, can't compute lambda.");
					lambda = trr/(n*n); //TODO: raise exception if n == 0
				}
				if(mhtp_params.constant_step_size) //TODO: cf above
					constant_step_size = false;
				if(is_verbose)
					std::cout << "end MHTP" << std::endl;
			}
			else
				update_fac(cur_fac, f_id);

			next_fid(); // f_id updated to iteration factor index (pL or pR too)
		}
		//update lambda
		//TODO: variable decl in parent scope
		Faust::MatDense<FPP,DEVICE> A_H_S = S.multiply(A_H);
		//		auto last_Sfac_vec = { *(S.begin()+nfacts-1), dynamic_cast<Faust::MatGeneric<FPP,DEVICE>*>(&A_H)};
		//		Faust::TransformHelper<FPP,DEVICE> A_H_S_(*pL[nfacts-1], last_Sfac_vec);
		//		A_H_S_.disable_dtor();
		//		Faust::MatDense<FPP,DEVICE> A_H_S = A_H_S_.get_product();
		Real<FPP> trr = std::real(A_H_S.trace());
		Real<FPP> n = S.normFro();
		if(std::numeric_limits<Real<FPP>>::epsilon() >= n)
			throw std::runtime_error("Faust Frobenius norm is zero, can't compute lambda.");
		lambda = trr/(n*n); //TODO: raise exception if n == 0
//		std::cout << "debug lambda: " << lambda << std::endl;
		if(is_verbose)
		{
			set_calc_err_ite_period(); //macro setting the variable ite_period
			if(! (i%ite_period))
			{
				std::cout << "PALM4MSA2020 iteration: " << i;
				auto err = calc_rel_err(S, A, lambda);
				std::cout << " relative error: " << err;
				std::cout << " (call id: " << id << ")" << std::endl;
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
		std::cout << "palm4msa spectral time=" << spectral_duration.count() << std::endl;
		std::cout << "palm4msa fgrad time=" << fgrad_duration.count() << std::endl;
	}
}

template <typename FPP, FDevice DEVICE>
void Faust::compute_n_apply_grad1(const int f_id, const Faust::MatDense<FPP,DEVICE> &A, Faust::TransformHelper<FPP,DEVICE>& S, std::vector<TransformHelper<FPP,DEVICE>*> &pL, std::vector<TransformHelper<FPP,DEVICE>*>& pR, const FPP& lambda, const Real<FPP> &c, Faust::MatDense<FPP,DEVICE> &out /* D */, const StoppingCriterion<Real<FPP>>& sc, Real<FPP> &error, const int prod_mod, const bool packing_RL)
{
	Faust::MatDense<FPP,DEVICE> tmp;
	Faust::MatDense<FPP,DEVICE> & D = out;
	Faust::MatDense<FPP,DEVICE> * LorR;
	Faust::MatDense<FPP,DEVICE> _LorR;
	auto S_j_vec = {*(S.begin()+f_id)};
	Faust::TransformHelper<FPP, DEVICE> _LSR(/*lambda_vec,*/ *pL[f_id], S_j_vec, *pR[f_id]);
	//			tmp = _LSR.get_product(prod_mod);
	_LSR.get_product(tmp, prod_mod);
	tmp *= FPP(lambda);
	tmp -= A;
	if(sc.isCriterionErr())
		error = tmp.norm();
	FPP alpha_R = 1, alpha_L = 1, beta_R = 0, beta_L = 0; //decl in parent scope
	if(pR[f_id]->size() > 0)
	{
		if(packing_RL)
			LorR = dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(pR[f_id]->get_gen_fact_nonconst(0)); //normally pR[f_id] is packed (hence reduced to a single MatDense)
		else
		{
			_LorR = pR[f_id]->get_product(prod_mod);
			LorR = &_LorR;
		}
		if(pL[f_id]->size() == 0)
		{ //no L factor for factor f_id
			alpha_R = - lambda/c;
			beta_R = 1;
			gemm(tmp, *LorR, D, alpha_R, beta_R, 'N', 'T');
		}
		else
			gemm(tmp, *LorR, tmp, alpha_R, beta_R, 'N', 'T');
	}
	if(pL[f_id]->size() > 0)
	{
		if(packing_RL)
			LorR = dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(pL[f_id]->get_gen_fact_nonconst(0));
		else
		{
			_LorR = pL[f_id]->get_product(prod_mod);
			LorR = &_LorR;
		}
		alpha_L = -lambda/c;
		beta_L = 1;
		gemm(*LorR, tmp, D, alpha_L, beta_L, 'T', 'N');
	}
}

template <typename FPP, FDevice DEVICE>
void Faust::compute_n_apply_grad2(const int f_id, const Faust::MatDense<FPP,DEVICE> &A, Faust::TransformHelper<FPP,DEVICE>& S, std::vector<TransformHelper<FPP,DEVICE>*> &pL, std::vector<TransformHelper<FPP,DEVICE>*>& pR, const FPP& lambda, const Real<FPP> &c, Faust::MatDense<FPP,DEVICE> &out /* D */, const StoppingCriterion<Real<FPP>>& sc, Real<FPP> &error, const int prod_mod, const bool packing_RL)
{
	Faust::MatDense<FPP,DEVICE> tmp;
	Faust::MatDense<FPP,DEVICE> & D = out;
	Faust::MatDense<FPP,DEVICE> *_L, *_R, __L, __R;
	Faust::MatDense<FPP,DEVICE> * LorR;
	Faust::MatDense<FPP,DEVICE> _LorR;
	std::vector<MatDense<FPP,DEVICE>*> facts;
	std::vector<char> tc_flags;
//#define mul_3_facts multiply_order_opt
#define mul_3_facts multiply_order_opt_all_ends// this one only optimizes the product on factor ends but for three factors it doesn't change anything comparing to multiply_order_opt
	tmp = A;
	if(pR[f_id]->size() > 0)
	{
		if(packing_RL)
			_R = dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(pR[f_id]->get_gen_fact_nonconst(0));
		else
		{
			__R = pR[f_id]->get_product(prod_mod);
			_R = &__R;
		}

	}
	if(pL[f_id]->size() > 0)
	{
		if(packing_RL)
			_L = dynamic_cast<Faust::MatDense<FPP,DEVICE>*>(pL[f_id]->get_gen_fact_nonconst(0));
		else
		{
			__L = pL[f_id]->get_product(prod_mod);
			_L = &__L;
		}
	}
	if(pR[f_id]->size() > 0 && pL[f_id]->size() > 0)
	{
		// compute error = m_lambda*L*S*R-data
		facts = { _L, &D, _R };
		mul_3_facts(facts, tmp, (FPP) lambda, (FPP) -1.0);
		if(sc.isCriterionErr())
			error = tmp.norm();
		// compute m_lambda/c * L'*error*R'
		facts = { _L, &tmp, _R };
		tc_flags = {'T', 'N', 'T'};
		mul_3_facts(facts, D, (FPP) - lambda/c, (FPP)1, tc_flags);
	}
	else if(pR[f_id]->size() > 0)
	{
		// compute error = m_lambda*L*S*R-data
		facts = { &D, _R };
		mul_3_facts(facts, tmp, (FPP) lambda, (FPP) -1.0);
		if(sc.isCriterionErr())
			error = tmp.norm();
		// compute m_lambda/c * L'*error*R'
		facts = { &tmp, _R };
		tc_flags = { 'N', 'T'};
		mul_3_facts(facts, D, (FPP) - lambda/c, (FPP)1, tc_flags);
	}
	else //if(pL[f_id]->size() > 0)
	{
		// compute error = m_lambda*L*S*R-data
		facts = { _L, &D};
		mul_3_facts(facts, tmp, (FPP) lambda, (FPP) -1.0);
		if(sc.isCriterionErr())
			error = tmp.norm();
		// compute m_lambda/c * L'*error*R'
		facts = { _L, &tmp};
		tc_flags = {'T', 'N'};
		mul_3_facts(facts, D, (FPP) - lambda/c, (FPP)1, tc_flags);
	}
}


